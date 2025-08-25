package com.leo.histomagnify;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.io.FileInfo;
import ij.plugin.PlugIn;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.io.*;
import java.net.URL;
import java.net.URLDecoder;
import java.nio.file.*;
import java.util.*;
import java.util.List;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

import ij.gui.Overlay;
import ij.gui.Roi;


import ij.gui.ImageRoi;
import ij.gui.TextRoi;
import ij.measure.Calibration;

public class HistoMagnify_Plugin implements PlugIn {
    private static final String SCRIPTS_ROOT = "scripts/";
    private static final String SRC_ROOT     = "src/";
    private static final String REQ_RESOURCE = "/scripts/requirements.txt";

    private static Path lastTmpDir;
    private static Map<String, File> lastSegFiles;
    private static Set<String> lastModelsRan = new HashSet<>();
    private static ImagePlus lastResult;

    //visibility toggles
    private static boolean showLegend = true;
    private static boolean showAxisDots = true;
    private static boolean showAxisColorBar = true;
    private static boolean showProcessLines = true;
    private static boolean showProcessColorBar = true;
    private static boolean showProcessContours = true;
    private static boolean showStatsBanner = true;



    private static final double DEFAULT_MIN_DIST_UM  = 0.08;
    private static final double DEFAULT_MAX_PAIR_UM  = 1.5;
    private static boolean thresholdsBootstrapped = false;


    private static Double userPixelSizeUnitsPerPixel = null;
    private static String userPixelUnit = "µm"; // display label only

    private static Double autoPixelSizeUnitsPerPixel = null;
    private static String autoPixelUnit = "µm";

    private static double EF_val = 3.5;
    private static boolean expanded = true;


    private static List<Roi> legendItems = new ArrayList<>();
    private static List<Roi> axisDotItems = new ArrayList<>();
    private static List<Roi> axisBarItems = new ArrayList<>();
    private static List<Roi> procLineItems = new ArrayList<>();
    private static List<Roi> procBarItems = new ArrayList<>();
    private static List<Roi> procContourItems = new ArrayList<>();
    private static List<Roi> rbcContourItems = new ArrayList<>();
    private static List<Roi> nucContourItems = new ArrayList<>();
    private static List<Roi> statsBannerItems = new ArrayList<>();

    private static boolean showRbcContours = true;
    private static boolean showNucContours = true;

    //remember which magnification the last result was built with
    private static String lastEffMag = null;

    private static boolean isWin = System.getProperty("os.name").toLowerCase(Locale.ROOT).contains("win");

    private static double wsProcMinDist   = 3.28;
    private static double wsProcThreshRel = 0.26;
    private static double wsProcSigma     = 2.0;
    // ize filtering
    private static double procKeepLow     = 0.20;
    private static double procKeepHigh    = 0.02;
    private static double procMaxPairPx   = 20.0;

    //20X RBC
    private static double wsRbcMinDist    = 15.0;
    private static double wsRbcThreshRel  = 0.30;
    private static double wsRbcSigma      = 2.0;
    private static double rbcKeepLow      = 0.40;
    private static double rbcKeepHigh     = 0.00;

    //20X Nuclei
    private static double wsNucMinDist    = 40.0;
    private static double wsNucThreshRel  = 0.30;
    private static double wsNucSigma      = 2.0;
    private static double nucKeepLow      = 0.40;
    private static double nucKeepHigh     = 0.00;

    private static Integer lastNucCount = null;
    private static Integer lastRbcCount = null;
    private static Double  lastGbmThicknessUm = null;
    private static Double  lastProcNndUm = null;


    private static final Map<String,String> MODEL_RESOURCE;
    static {
        Map<String,String> m = new HashMap<>();
        m.put("ACTN4",      "/models/ACTN4.hdf5");
        m.put("DAPI",       "/models/DAPI.hdf5");
        m.put("NHS_SINGLE_CHANNEL", "/models/NHS_ester_single.hdf5");
        m.put("NHS_COMBINED_ACTN4",    "/models/NHS_ester_com.hdf5");
        MODEL_RESOURCE = Collections.unmodifiableMap(m);
    }

    private static final Map<String, String> BF_MODELS;
    static {
        Map<String,String> m = new HashMap<>();
        m.put("20X", "/models/20x.hdf5");
        m.put("40X", "/models/40x.hdf5");
        BF_MODELS = Collections.unmodifiableMap(m);
    }

    //state for BF pipeline
    private static boolean standardProcessed = false;
    private static boolean showEnhancedBackground = true;
    private static Integer lastPageSelected = 0;
    private static String autoMagnification = null;
    private static String userMagnification = null;           // if user overrides, store here

    private static File lastSeg20x;
    private static File lastSeg40x;
    private static File lastEnhanced;


    @Override
    public void run(String arg) {
        if ("settings".equalsIgnoreCase(arg)) {
            openSettingsDialog();
            return;
        }

//        if ("thickness".equalsIgnoreCase(arg)) {
//            runThicknessMenu();
//            return;
//        }
//        if ("proc".equalsIgnoreCase(arg) || "process".equalsIgnoreCase(arg)) {
//            runProcessMenu();
//            return;
//        }

        if ("thickness_full".equalsIgnoreCase(arg) || "thickness/full".equalsIgnoreCase(arg)) {
            runThicknessMenu(false);
            return;
        }
        if ("thickness_roi".equalsIgnoreCase(arg) || "thickness/roi".equalsIgnoreCase(arg)) {
            runThicknessMenu(true);
            return;
        }

        if ("proc_full".equalsIgnoreCase(arg) || "process_full".equalsIgnoreCase(arg) || "proc/full".equalsIgnoreCase(arg)) {
            runProcessMenu(false);
            return;
        }
        if ("proc_roi".equalsIgnoreCase(arg) || "process_roi".equalsIgnoreCase(arg) || "proc/roi".equalsIgnoreCase(arg)) {
            runProcessMenu(true);
            return;
        }



        if ("nuc_full".equalsIgnoreCase(arg)) {
            runCountMenu("nuc", false);
            return;
        }
        if ("nuc_roi".equalsIgnoreCase(arg))  {
            runCountMenu("nuc", true);
            return;
        }
        if ("rbc_full".equalsIgnoreCase(arg)) {
            runCountMenu("rbc", false);
            return;
        }
        if ("rbc_roi".equalsIgnoreCase(arg))  {
            runCountMenu("rbc", true);
            return;
        }



        //list open images
        int count = WindowManager.getImageCount();
        if (count == 0) {
            IJ.showMessage("MagnifySeg", "No images open.");
            return;
        }
        String[] titles = new String[count];
        for (int i = 0; i < count; i++) {
            titles[i] = WindowManager.getImage(i + 1).getTitle();
        }

        //detect magnification
        ImagePlus first = WindowManager.getImage(titles[0]);
        autoMagnification = detectMagnification(first); // may be null

        GenericDialog gd = new GenericDialog("HistoMagnify – Apply Segmentation");
        gd.addChoice("Source image:", titles, titles[0]);
        gd.addNumericField("TIFF page (0-based):", lastPageSelected != null ? lastPageSelected : 0, 0);


        String[] mags = new String[]{"20X", "40X"};
        String autoMag = autoMagnification;
        if ("60X".equals(autoMag) || autoMag == null) autoMag = "40X";
        int defaultMagIdx = "20X".equals(autoMag) ? 0 : 1;
        gd.addChoice("Magnification:", mags, mags[defaultMagIdx]);

        gd.addCheckbox("Image is standard processed (skip enhancement)", standardProcessed);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        String sourceTitle = gd.getNextChoice();
        int pageIndex0 = (int) gd.getNextNumber();
        String magChoice   = gd.getNextChoice();

        boolean skipEnhancement = gd.getNextBoolean();
        standardProcessed = skipEnhancement;

        String effMag = magChoice;
        userMagnification = effMag;
        lastPageSelected = pageIndex0;


        ImagePlus original = WindowManager.getImage(sourceTitle);
        if (original == null) {
            IJ.showMessage("HistoMagnify","Can't find image: "+sourceTitle);
            return;
        }
        refreshAutoPixelSizeFrom(original);
        int zPages = Math.max(1, original.getNSlices());
        pageIndex0 = Math.min(Math.max(0, pageIndex0), zPages - 1);
        // Force to set a real world pixel size
        Calibration detCal = original.getCalibration();
        String detectedUnit = "";

        if (detCal != null && detCal.getUnit() != null) {
            detectedUnit = detCal.getUnit().trim();
        } else if (autoPixelUnit != null) {
            detectedUnit = autoPixelUnit;
        }


        //block + open Settings
        boolean missingSize = (userPixelSizeUnitsPerPixel == null || Double.isNaN(userPixelSizeUnitsPerPixel)
                || userPixelSizeUnitsPerPixel <= 0);

        if (unitLooksLikePixel(detectedUnit) || unitLooksLikePixel(userPixelUnit) || missingSize) {
            IJ.showMessage("HistoMagnify",
                    "This image appears to be uncalibrated (unit = 'pixel').\n" +
                            "Please set a real-world pixel size and unit in Settings (e.g., \u00B5m/pixel) before running.");
            openSettingsDialog(true);

//            // Re-check
//            missingSize = (userPixelSizeUnitsPerPixel == null
//                    || Double.isNaN(userPixelSizeUnitsPerPixel)
//                    || userPixelSizeUnitsPerPixel <= 0);
//            if (unitLooksLikePixel(userPixelUnit) || missingSize) {
//                IJ.showMessage("HistoMagnify",
//                        "Segmentation aborted: pixel size and unit must be set to real units (e.g., \u00B5m/pixel).");
//                return;
//            }
        }



        try {
            // ensure venv
            Path home = Paths.get(System.getProperty("user.home"));
            Path baseCache = home.resolve(".histomagnify");
            Path venvDir = baseCache.resolve("venv");
            String pyInVenv = venvDir.resolve(isWin ? "Scripts\\python.exe" : "bin/python").toString();

            boolean needCreate = true;
            if (Files.exists(venvDir)) {
                if (Files.exists(Paths.get(pyInVenv))) {
                    needCreate = false;
                    IJ.log("[HistoMagnify] Using venv at " + venvDir);
                } else {
                    IJ.log("[HistoMagnify] Detected incomplete venv, recreating…");
                    try {
                        Files.walk(venvDir).sorted(Comparator.reverseOrder())
                                .forEach(p -> { try { Files.deleteIfExists(p); } catch (IOException ignored) {} });
                    } catch (IOException ignored) {}
                }
            }

            if (needCreate) {
                Files.createDirectories(baseCache);

                String[] pyLauncher = choosePythonLauncher();
                if (pyLauncher == null) {
                    IJ.showMessage("HistoMagnify",
                            "Python 3.8+ wasn’t found on this system.\n\n" +
                                    "Install Python 3.8+ (recommended via python.org),\n" +
                                    "then re-run this command.\n\n" +
                                    "Recommended source https://www.python.org/downloads/");
                    return;
                }

                int vcode = runAndGetExit(cat(pyLauncher, "-m", "venv", venvDir.toString()),
                        "[venv] ", null);
                if (vcode != 0) {
                    IJ.showMessage("HistoMagnify",
                            "Failed to create a Python virtual environment.\n" +
                                    "Please ensure Python 3.8+ is installed and try again.");
                    return;
                }

                // write requirements and install
                Path req = baseCache.resolve("requirements.txt");
                try (InputStream in = getClass().getResourceAsStream(REQ_RESOURCE);
                     OutputStream os = Files.newOutputStream(req)) {
                    if (in == null) throw new FileNotFoundException("Missing " + REQ_RESOURCE + " in plugin JAR");
                    byte[] buf = new byte[8192]; int r;
                    while ((r = in.read(buf)) > 0) os.write(buf, 0, r);
                }

                runAndGetExit(new String[]{pyInVenv, "-m", "pip", "install", "--upgrade", "pip"}, "[pip] ", null);
                int pcode = runAndGetExit(new String[]{pyInVenv, "-m", "pip", "install", "-r", req.toString()},
                        "[pip] ", null);
                if (pcode != 0) {
                    IJ.showMessage("HistoMagnify",
                            "Python was found, but psackage installation failed.\n" +
                                    "Open the Log window for details and check internet access.");
                    return;
                }
            }


            Path tmpDir = Files.createTempDirectory("histomagnify_");

            File rawPageSnap = exportSinglePageRAW(original, pageIndex0, tmpDir);


            String origPath = guessTifPath(original);
            String tifPath = (origPath != null) ? origPath : rawPageSnap.getAbsolutePath();

            String modelRes = BF_MODELS.get(effMag.equals("20X") ? "20X" : "40X");
            if (modelRes == null) modelRes = BF_MODELS.get("40X");
            String modelName = Paths.get(modelRes).getFileName().toString();
            File modelFile = extractResource(tmpDir, modelRes, modelName);

            extractFolder(tmpDir, SCRIPTS_ROOT);

            // run segmentation
            String py  = venvDir.resolve(isWin ? "Scripts\\python.exe" : "bin/python").toString();
            File segOut = tmpDir.resolve(effMag.equals("20X") ? "seg_20x.tif" : "seg_40x.tif").toFile();
            File enhOut = tmpDir.resolve("enhanced.tif").toFile();

            String enhanced = "1";
            if (standardProcessed){
                enhanced = "0";
            }


            List<String> cmd = new ArrayList<>(Arrays.asList(
                    py, tmpDir.resolve("segment.py").toString(),
                    "--tif", tifPath,
                    "--page", String.valueOf(pageIndex0),
                    "--model", effMag.equals("20X") ? "20x" : "40x",
                    "--modeldir", tmpDir.toString(),
                    "--output", segOut.getAbsolutePath(),
                    "--save_enhanced", enhOut.getAbsolutePath(),
                    "--enhance", enhanced
            ));
            IJ.log("[HistoMagnify] Running segmentation ("+effMag+")...");
            runWithLogging(cmd.toArray(new String[0]), "[HistoMagnify] ", tmpDir.toFile());
            if (!segOut.isFile()) {
                IJ.showMessage("HistoMagnify","Segmentation failed.");
                return;
            }

            //for secondary tools
            if ("20X".equals(effMag)) lastSeg20x = segOut; else lastSeg40x = segOut;
            lastEnhanced = enhOut.isFile() ? enhOut : null;
            lastTmpDir = tmpDir;

            // 2-slice result: [Background, Overlay]
            ColorProcessor bgRGB;
            if (showEnhancedBackground && lastEnhanced != null && lastEnhanced.isFile()) {
                ImagePlus enh = IJ.openImage(lastEnhanced.getAbsolutePath());
                bgRGB = toRGB(enh.getProcessor());
                enh.close();
            } else {
                File rawFile = tmpDir.resolve("input_page_raw.tif").toFile();
                ImagePlus rawImp = IJ.openImage(rawFile.getAbsolutePath());
                bgRGB = toRGB(rawImp.getProcessor());
                rawImp.close();
            }


            ImagePlus seg = IJ.openImage(segOut.getAbsolutePath());
            ImageProcessor lbl = seg.getProcessor().convertToByteProcessor();
            seg.close();

            int w = bgRGB.getWidth();
            int h = bgRGB.getHeight();
            ColorProcessor out = (ColorProcessor) bgRGB.duplicate();


            int[][] lut20 = new int[][]{
                    {  0,  0,130},
                    {  91, 24,199},
                    {  242,91,96},
                    {240,203, 73},
                    {89,195,  71},
                    {76,  98,  246}
            };

            int[][] lut40 = new int[][] {
                    {105,105,105},
                    {255,  0,255},
                    {  0,255,255}
            };

            int[][] lut = "20X".equals(effMag) ? lut20 : lut40;
            float alpha = 0.45f;

            for (int y=0; y<h; y++) {
                for (int x=0; x<w; x++) {
                    int lab = lbl.get(x,y) & 0xff;
                    if (lab <= 0 || lab >= lut.length) continue;
                    int[] c = lut[lab];
                    int p  = out.get(x,y);
                    int r=(p>>16)&255, g=(p>>8)&255, b=p&255;
                    int nr = (int)(r*(1-alpha) + c[0]*alpha);
                    int ng = (int)(g*(1-alpha) + c[1]*alpha);
                    int nb = (int)(b*(1-alpha) + c[2]*alpha);
                    out.set(x,y, ((nr&255)<<16)|((ng&255)<<8)|(nb&255));
                }
            }

            ij.ImageStack st = new ij.ImageStack(w,h);
            st.addSlice("Background", bgRGB);
            st.addSlice("Overlay", out);
            ImagePlus result = new ImagePlus("HistoMagnify – Result ("+effMag+")", st);

            //Build legend
            Overlay ov = new Overlay();
            legendItems.clear();
            int x0 = 25, y0 = 25, box = 36;
            Font f = new Font("SansSerif", Font.PLAIN, 28);
            String[] lbls20 = {"Background","Nucleus","RBC","Tube","Glomerulus","GBM"};
            String[] lbls40 = {"Background","GBM","Podocyte foot processes"};
            String[] labels = "20X".equals(effMag) ? lbls20 : lbls40;


            int nItems = Math.min(lut.length, labels.length) - 1;
            int maxText = 0;
            for (int i = 1; i <= nItems; i++) {
                TextRoi tmp = new TextRoi(0, 0, labels[i], f);
                maxText = Math.max(maxText, (int)Math.round(tmp.getFloatWidth()));
            }
            int pad = 12;
            int totalH = nItems * box + (nItems - 1) * 14 + pad * 2;
            int panelW = box + 16 + maxText + pad * 2;
            int panelX = x0 - pad;
            int panelY = y0 - pad;

            ImageRoi legendBg = makeBackdrop(panelX, panelY, panelW, totalH, 0.35f);
            ov.add(legendBg);
            legendItems.add(legendBg);

            //egend swatches/labels
            int yy = y0;
            for (int i = 1; i <= nItems; i++) {
                int[] c = lut[i];
                Roi r = new Roi(x0, yy, box, box);
                Color cc = new Color(c[0], c[1], c[2]);
                r.setFillColor(cc);
                r.setStrokeColor(cc);
                ov.add(r);
                legendItems.add(r);

                TextRoi t = new TextRoi(x0 + box + 16, yy, labels[i], f);
                t.setStrokeColor(Color.WHITE);
                ov.add(t);
                legendItems.add(t);

                yy += box + 14;
            }
            // reset, apply visibility
            axisDotItems.clear();
            axisBarItems.clear();
            procLineItems.clear();
            procBarItems.clear();
            procContourItems.clear();
            rbcContourItems.clear();
            nucContourItems.clear();

            lastNucCount = null;
            lastRbcCount = null;
            lastGbmThicknessUm = null;
            lastProcNndUm = null;

            //bottom stats banner
            rebuildStatsBanner(ov);


            applyVisibility(ov);
            result.setOverlay(ov);
            result.show();
            result.setSlice(2);


            lastResult = result;
            lastModelsRan = new HashSet<>(Collections.singletonList(effMag)); // remember which mag we ran
            lastEffMag = effMag; // for  refreshResultFromSettings() yto pivk the right LUT/seg

        } catch (Exception e) {
            IJ.handleException(e);
        }

    }

    private void extractFolder(Path tmpDir, String folderName) throws IOException {
        URL srcUrl = getClass().getProtectionDomain()
                .getCodeSource().getLocation();
        String jarPath = URLDecoder.decode(srcUrl.getPath(), "UTF-8");
        try (JarFile jar = new JarFile(jarPath)) {
            for (JarEntry entry : Collections.list(jar.entries())) {
                String name = entry.getName();
                if (!name.startsWith(folderName)) continue;
                String rel = name.substring(folderName.length());
                if (rel.isEmpty()) continue;
                File out = tmpDir.resolve(rel).toFile();
                if (entry.isDirectory()) out.mkdirs();
                else {
                    out.getParentFile().mkdirs();
                    try (InputStream in = jar.getInputStream(entry);
                         OutputStream os = new FileOutputStream(out)) {
                        byte[] buf = new byte[8192]; int r;
                        while ((r=in.read(buf))>0) os.write(buf,0,r);
                    }
                    if (rel.endsWith(".py")) out.setExecutable(true,false);
                }
            }
        }
    }

    private File extractResource(Path tmpDir, String resPath, String name) throws IOException {
        try (InputStream in = getClass().getResourceAsStream(resPath)) {
            if (in==null) throw new FileNotFoundException(resPath);
            File out = tmpDir.resolve(name).toFile();
            out.getParentFile().mkdirs();
            try (OutputStream os = new FileOutputStream(out)) {
                byte[] buf = new byte[8192]; int len;
                while ((len=in.read(buf))>0) os.write(buf,0,len);
            }
            return out;
        }
    }


    private void runWithLogging(String[] cmd, String prefix, File workDir)
            throws IOException, InterruptedException {
        Process p = Runtime.getRuntime().exec(cmd, null, workDir);
        Thread outG = new StreamGobbler(p.getInputStream(),prefix);
        Thread errG = new StreamGobbler(p.getErrorStream(), prefix);
        outG.start();
        errG.start();
        p.waitFor();
        outG.join();
        errG.join();
    }

    private static class StreamGobbler extends Thread {
        private final InputStream is;
        private final String prefix;
        StreamGobbler(InputStream is, String prefix) {
            this.is = is; this.prefix = prefix;
        }
        @Override
        public void run() {
            try (BufferedReader r = new BufferedReader(new InputStreamReader(is))) {
                String line;
                while ((line = r.readLine())!=null) {
                    IJ.log(prefix + line);
                }
            } catch (IOException ignored) {}
        }
    }




    private void makeBinaryMaskFromLabel(File labelTif, Roi roi, int labelVal, File outTif) throws IOException {
        ImagePlus imp = IJ.openImage(labelTif.getAbsolutePath());
        if (imp == null) throw new IOException("Open failed: "+labelTif);
        if (roi != null) imp.setRoi(roi);
        ImageProcessor ip = (roi != null ? imp.getProcessor().crop() : imp.getProcessor()).convertToByteProcessor();
        int w = ip.getWidth(), h = ip.getHeight();
        byte[] pix = (byte[]) ip.getPixels();
        for (int i=0;i<w*h;i++) {
            int v = pix[i] & 0xff;
            pix[i] = (byte)((v==labelVal)?255:0);
        }
        new ij.io.FileSaver(new ImagePlus("mask", ip)).saveAsTiff(outTif.getAbsolutePath());
        imp.close();
    }



    private void runThicknessMenu(boolean useROI) {
        //basic guard
        if (lastSeg40x == null || !lastSeg40x.isFile()) {
            IJ.showMessage("MagnifySeg",
                    "GBM thickness can only be run on 40X images.\n" +
                            "Please run the 40X brightfield segmentation first.");
            return;
        }
        if (lastTmpDir == null) lastTmpDir = lastSeg40x.getParentFile().toPath();
        if (useROI && (lastResult == null)) {
            IJ.showMessage("MagnifySeg",
                    "No result window found to draw an ROI.\n" +
                            "Re-run segmentation (or use the full-frame command).");
            return;
        }


        if (userPixelSizeUnitsPerPixel == null || userPixelSizeUnitsPerPixel <= 0) {
            IJ.showMessage("MagnifySeg",
                    "Pixel size is not set.\nOpen Settings and enter a pixel size (or reset to auto) before running thickness.");
            openSettingsDialog(true);
            if (userPixelSizeUnitsPerPixel == null || userPixelSizeUnitsPerPixel <= 0) return;
        }


        Roi roi = null;
        if (useROI) {
            if (lastResult.getWindow() == null) lastResult.show();
            lastResult.getWindow().toFront();
            IJ.selectWindow(lastResult.getTitle());
            IJ.setTool("rectangle");
            new ij.gui.WaitForUserDialog("ROI",
                    "Draw a rectangle on the result and click OK.").show();
            roi = lastResult.getRoi();
            if (roi == null || roi.getBounds().width == 0 || roi.getBounds().height == 0) {
                IJ.showMessage("MagnifySeg", "No ROI drawn.");
                return;
            }
        } else {
            if (lastResult != null) lastResult.killRoi();
        }

        // 4) Build GBM mask from 40X labels (GBM = 1)
        File gbmSeg  = lastSeg40x;
        int  gbmLabel = 1;

        File gbmMask  = lastTmpDir.resolve("gbm_mask.tif").toFile();
        File thickTxt = lastTmpDir.resolve("thickness.txt").toFile();
        File thickCsv = lastTmpDir.resolve("thickness_points.csv").toFile();
        try {
            makeBinaryMaskFromLabel(gbmSeg, roi, gbmLabel, gbmMask);
        } catch (IOException ex) {
            IJ.handleException(ex); return;
        }

        Path venvDir = Paths.get(System.getProperty("user.home")).resolve(".histomagnify").resolve("venv");
        String pyExe = venvDir.resolve(isWin ? "Scripts\\python.exe" : "bin/python").toString();

        File metricsPy = lastTmpDir.resolve("metrics.py").toFile();

        String[] cmdT = new String[]{
                pyExe, metricsPy.getAbsolutePath(),
                "--task","thickness",
                "--mask", gbmMask.getAbsolutePath(),
                "--out_txt", thickTxt.getAbsolutePath(),
                "--out_csv", thickCsv.getAbsolutePath()
        };
        try {
            IJ.log("[MagnifySeg] Computing GBM thickness...");
            runWithLogging(cmdT, "[metrics:thick] ", lastTmpDir.toFile());

        } catch (Exception ex) {
            IJ.handleException(ex); return;
        }

        //overlay points and show value
        Overlay ov2 = lastResult.getOverlay();
        if (ov2 == null) {
            ov2 = new Overlay();
            lastResult.setOverlay(ov2);
        }
        for (Roi r: axisDotItems) ov2.remove(r);
        for (Roi r: axisBarItems) ov2.remove(r);
        axisDotItems.clear();
        axisBarItems.clear();

        double val = Double.NaN;
        List<double[]> pts = new ArrayList<>();
        try {
            try {
                String s = new String(java.nio.file.Files.readAllBytes(thickTxt.toPath())).trim();
                double unitsPerPixel = unitsPerPixel(userPixelSizeUnitsPerPixel, expanded, EF_val);
                val = Double.parseDouble(s) * unitsPerPixel;   // px -> units
            } catch (Exception ignore) {}

            try (BufferedReader br = new BufferedReader(new FileReader(thickCsv))) {
                String line;
                java.awt.Rectangle b = (roi != null) ? roi.getBounds() : new java.awt.Rectangle(0,0,0,0);
                while ((line = br.readLine()) != null) {
                    String[] t = line.split(",");
                    if (t.length >= 2) {
                        double x = Double.parseDouble(t[0]) + b.x;
                        double y = Double.parseDouble(t[1]) + b.y;
                        double unitsPerPixel = unitsPerPixel(userPixelSizeUnitsPerPixel, expanded, EF_val);
                        double tvUnits = (t.length >= 3) ? (Double.parseDouble(t[2]) * unitsPerPixel) : Double.NaN;
                        pts.add(new double[]{x, y, tvUnits});
                    }
                }
            }
        } catch (Exception ignore) {}

        if (!pts.isEmpty()) {
            double vmin = Double.POSITIVE_INFINITY, vmax = Double.NEGATIVE_INFINITY;
            boolean haveVals = false;
            for (double[] p : pts) {
                if (!Double.isNaN(p[2])) {
                    haveVals = true;
                    vmin = Math.min(vmin, p[2]);
                    vmax = Math.max(vmax, p[2]);
                }
            }
            if (!haveVals) { vmin = 0; vmax = 1; }

            for (double[] p : pts) {
                double v = haveVals ? p[2] : 0.5;
                double t = (vmax > vmin) ? (v - vmin) / (vmax - vmin) : 0.0;
                int[] rgb = viridis(t);
                int cx = (int) Math.round(p[0]);
                int cy = (int) Math.round(p[1]);
                ij.gui.OvalRoi dot = new ij.gui.OvalRoi(cx - 1, cy - 1, 2, 2);
                Color col = new Color(rgb[0], rgb[1], rgb[2]);
                dot.setFillColor(col);
                dot.setStrokeColor(col);
                dot.setPosition(2);
                ov2.add(dot);
                axisDotItems.add(dot);
            }
            List<Roi> bar = makeColorBar(ov2,
                    (int)(lastResult.getWidth()-110), 40, 20, Math.max(120, lastResult.getHeight()/4),
                    vmin, vmax, "Thickness (" + userPixelUnit + ")", "viridis");
            axisBarItems.addAll(bar);
        }

        lastGbmThicknessUm = Double.isNaN(val) ? null : val;


        rebuildStatsBanner(ov2);
        applyVisibility(ov2);
        lastResult.updateAndDraw();
        IJ.showMessage("MagnifySeg", "Avg GBM thickness: " + val + " " + userPixelUnit);
    }


    private void runProcessMenu(boolean useROI)  {
        if (lastSeg40x == null || !lastSeg40x.isFile()) {
            IJ.showMessage("MagnifySeg",
                    "Process distance requires a 40X brightfield segmentation.\n" +
                            "Please run the 40X segmentation first.");
            return;
        }

        if (lastTmpDir == null) lastTmpDir = lastSeg40x.getParentFile().toPath();

        if (useROI && (lastResult == null)) {
            IJ.showMessage("MagnifySeg",
                    "No result window found to draw an ROI.\n" +
                            "Re-run segmentation (or use the full-frame command).");
            return;
        }


        if (userPixelSizeUnitsPerPixel == null || userPixelSizeUnitsPerPixel <= 0) {
            IJ.showMessage("MagnifySeg",
                    "Pixel size is not set.\nOpen Settings and enter a pixel size (or reset to auto) before running process distance.");
            openSettingsDialog(true);
            if (userPixelSizeUnitsPerPixel == null || userPixelSizeUnitsPerPixel <= 0) return;
        }
        double pxSize = userPixelSizeUnitsPerPixel;

        double maxPairPx = procMaxPairPx;

        Roi roi = null;
        if (useROI) {
            if (lastResult.getWindow() == null) lastResult.show();
            lastResult.getWindow().toFront();
            IJ.selectWindow(lastResult.getTitle());
            IJ.setTool("rectangle");
            new ij.gui.WaitForUserDialog("ROI",
                    "Draw a rectangle on the result and click OK.").show();
            roi = lastResult.getRoi();
            if (roi == null || roi.getBounds().width == 0 || roi.getBounds().height == 0) {
                IJ.showMessage("MagnifySeg", "No ROI drawn.");
                return;
            }
        } else {
            if (lastResult != null) lastResult.killRoi();
        }


        File procSrc = lastSeg40x;
        File procMask = lastTmpDir.resolve("proc_mask.tif").toFile();
        File procTxt  = lastTmpDir.resolve("proc.txt").toFile();
        File procCsv  = lastTmpDir.resolve("proc_pairs.csv").toFile();
        File procLabels = lastTmpDir.resolve("proc_labels.tif").toFile();
        File procEdges  = lastTmpDir.resolve("proc_contours.tif").toFile();
        File procOuter  = lastTmpDir.resolve("proc_outer_contours.tif").toFile();

        try {
            makeBinaryMaskFromLabel(procSrc, roi, 2, procMask);
        } catch (IOException ex) {
            IJ.handleException(ex); return;
        }

        //metrics.py
        Path venvDir = Paths.get(System.getProperty("user.home"))
                .resolve(".histomagnify").resolve("venv");
        String pyExe = venvDir.resolve(isWin ? "Scripts\\python.exe" : "bin/python").toString();

        File metricsPy = lastTmpDir.resolve("metrics.py").toFile();

        double px = (expanded && EF_val > 0) ? (pxSize / EF_val) : pxSize;

        String[] cmdP = new String[]{
                pyExe, metricsPy.getAbsolutePath(),
                "--task","proc",
                "--mask", procMask.getAbsolutePath(),
                "--max_pair_px", String.valueOf(maxPairPx),

                "--ws_min_dist",   String.valueOf(wsProcMinDist),
                "--ws_thresh_rel", String.valueOf(wsProcThreshRel),
                "--ws_sigma",      String.valueOf(wsProcSigma),

                "--keep_low",  String.valueOf(procKeepLow),
                "--keep_high", String.valueOf(procKeepHigh),

                "--out_labels",   procLabels.getAbsolutePath(),
                "--out_contours", procEdges.getAbsolutePath(),
                "--out_outer_contours", procOuter.getAbsolutePath(),

                "--out_txt", procTxt.getAbsolutePath(),
                "--out_csv", procCsv.getAbsolutePath()
        };

        try {
            IJ.log("[MagnifySeg] Computing process NND...");
            runWithLogging(cmdP, "[metrics:proc] ", lastTmpDir.toFile());
        } catch (Exception ex) {
            IJ.handleException(ex); return;
        }

        //verlay + colorbar
        Overlay ov2 = lastResult.getOverlay();
        if (ov2 == null) { ov2 = new Overlay(); lastResult.setOverlay(ov2); }

        for (Roi r: procLineItems) ov2.remove(r);
        for (Roi r: procBarItems) ov2.remove(r);
        for (Roi r: procContourItems) ov2.remove(r);
        procLineItems.clear();
        procBarItems.clear();
        procContourItems.clear();

        try {
            java.awt.Rectangle b = (roi != null) ? roi.getBounds() : new java.awt.Rectangle(0,0,0,0);

            Roi edgesOuter = makeColoredMaskRoi(procOuter, new Color(255, 240, 6), 1.0f, 1);
            edgesOuter.setPosition(2);
            if (roi != null) edgesOuter.setLocation(b.x, b.y);
            ov2.add(edgesOuter);
            procContourItems.add(edgesOuter);

            Roi edgesSplit = makeColoredMaskRoi(procEdges, new Color(255, 240, 6), 1.0f, 1);
            edgesSplit.setPosition(2);
            if (roi != null) edgesSplit.setLocation(b.x, b.y);
            ov2.add(edgesSplit);
            procContourItems.add(edgesSplit);

        } catch (IOException ignore) {
            IJ.log("[MagnifySeg] No process contours produced.");
        }

        double val = Double.NaN;
        class Pair { double x0,y0,x1,y1,distPx; }
        List<Pair> pairs = new ArrayList<>();
        try {
            try {
                String s = new String(java.nio.file.Files.readAllBytes(procTxt.toPath())).trim();
                double unitsPerPixel = unitsPerPixel(userPixelSizeUnitsPerPixel, expanded, EF_val);
                val = Double.parseDouble(s) * unitsPerPixel;   // px -> units

            } catch (Exception ignore) {}
            try (BufferedReader br = new BufferedReader(new FileReader(procCsv))) {
                String line;
                java.awt.Rectangle b = (roi != null) ? roi.getBounds() : new java.awt.Rectangle(0,0,0,0);
                while ((line = br.readLine()) != null) {
                    String[] t = line.split(",");
                    if (t.length >= 4) {
                        double px0 = Double.parseDouble(t[0]) + b.x;
                        double py0 = Double.parseDouble(t[1]) + b.y;
                        double px1 = Double.parseDouble(t[2]) + b.x;
                        double py1 = Double.parseDouble(t[3]) + b.y;
                        Pair p = new Pair();
                        p.x0=px0; p.y0=py0; p.x1=px1; p.y1=py1;
                        p.distPx = Math.hypot(px1-px0, py1-py0);
                        pairs.add(p);
                    }
                }
            }
        } catch (Exception ignore) {}

        if (!pairs.isEmpty()) {
            double dmin = Double.POSITIVE_INFINITY, dmax = Double.NEGATIVE_INFINITY;
            for (Pair p : pairs) {
                double d = p.distPx * px;
                dmin = Math.min(dmin, d);
                dmax = Math.max(dmax, d);
            }
            if (!(dmax > dmin)) { dmin = 0; dmax = Math.max(1e-9, dmax); }

            int cap = 0;
            for (Pair p : pairs) {
                if (cap++ > 2000) break;
                double d = p.distPx * px;
                double t = (dmax > dmin) ? (d - dmin) / (dmax - dmin) : 0.0;
                int[] rgb = hot(t);
                ij.gui.Line ln = new ij.gui.Line(p.x0, p.y0, p.x1, p.y1);
                ln.setStrokeColor(new Color(rgb[0], rgb[1], rgb[2]));
                ln.setStrokeWidth(1.0);
                ln.setPosition(2);
                ov2.add(ln);
                procLineItems.add(ln);
            }

            List<Roi> bar = makeColorBar(ov2,
                    (int)(lastResult.getWidth()-110), Math.max(120, lastResult.getHeight()/4)+100, 20, Math.max(120, lastResult.getHeight()/4),
                    dmin, dmax, "Proc dist (" + userPixelUnit + ")", "hot");
            procBarItems.addAll(bar);
        }

        lastProcNndUm = Double.isNaN(val) ? null : val;
        rebuildStatsBanner(ov2);
        applyVisibility(ov2);
        lastResult.updateAndDraw();

        IJ.showMessage("MagnifySeg", "Process mean nearest-neighbor distance: " + val + " " + userPixelUnit);
    }

    private void openSettingsDialog() { openSettingsDialog(false); }


    private void openSettingsDialog(boolean forcePixelFocus) {
        GenericDialog gd = new GenericDialog("HistoMagnify Settings");


        double shownPx = (userPixelSizeUnitsPerPixel != null) ? userPixelSizeUnitsPerPixel : Double.NaN;
        String unitShown = (userPixelUnit != null && !userPixelUnit.trim().isEmpty()) ? userPixelUnit : "units";
        final boolean fieldsInMicrons = !Double.isNaN(shownPx) && !"units".equals(unitShown);
        final double efDisp = expanded ? EF_val : 1.0;

        double dispProcMin = wsProcMinDist;
        double dispProcMax = procMaxPairPx;
        double dispRbcMin  = wsRbcMinDist;
        double dispNucMin  = wsNucMinDist;
        if (fieldsInMicrons) {
            dispProcMin = wsProcMinDist * (userPixelSizeUnitsPerPixel / efDisp);
            dispProcMax = procMaxPairPx * (userPixelSizeUnitsPerPixel / efDisp);
            dispRbcMin  = wsRbcMinDist * (userPixelSizeUnitsPerPixel / efDisp);
            dispNucMin  = wsNucMinDist * (userPixelSizeUnitsPerPixel / efDisp);
        }


        Panel twoCol = new Panel(new GridLayout(1, 2, 10, 0));

        //Left column
        BorderedPanel left = new BorderedPanel(new GridBagLayout());
        GridBagConstraints L = gbc(0, 0, 1.0, 0.0);
        left.add(headingLabel("General"), L);
        L.gridy++;
        Checkbox cbShowLegend = new Checkbox("Show legend", showLegend);
        left.add(indent(cbShowLegend, 14), L);
        L.gridy++;
        Checkbox cbShowBanner = new Checkbox("Show stats banner", showStatsBanner);
        left.add(indent(cbShowBanner, 14), L);

        L.gridy++; left.add(separator(), L);

        // 40X visibility
        L.gridy++;
        left.add(headingLabel("40X – Display / Overlay"), L);
        L.gridy++;
        Checkbox cbAxisDots = new Checkbox("Show GBM axis", showAxisDots);
        left.add(indent(cbAxisDots, 14), L);
        L.gridy++;
        Checkbox cbAxisBar = new Checkbox("Show GBM color bar", showAxisColorBar);
        left.add(indent(cbAxisBar, 14), L);
        L.gridy++;
        Checkbox cbProcLines = new Checkbox("Show process lines", showProcessLines);
        left.add(indent(cbProcLines, 14), L);
        L.gridy++;
        Checkbox cbProcBar = new Checkbox("Show process color bar", showProcessColorBar);
        left.add(indent(cbProcBar, 14), L);
        L.gridy++;
        Checkbox cbProcContours = new Checkbox("Show process contours", showProcessContours);
        left.add(indent(cbProcContours, 14), L);

        L.gridy++; left.add(separator(), L);

        // 20X visibility
        L.gridy++;
        left.add(headingLabel("20X – Display / Overlay"), L);
        L.gridy++;
        Checkbox cbRbcContours = new Checkbox("Show RBC contours", showRbcContours);
        left.add(indent(cbRbcContours, 14), L);
        L.gridy++;
        Checkbox cbNucContours = new Checkbox("Show nuclei contours", showNucContours);
        left.add(indent(cbNucContours, 14), L);

        L.gridy++; left.add(separator(), L);

        // Background/expansion toggles
        L.gridy++;
        left.add(headingLabel("Background & Expansion"), L);
        L.gridy++;
        Checkbox cbShowEnh = new Checkbox("Show enhanced/standard image as background", showEnhancedBackground);
        left.add(indent(cbShowEnh, 14), L);
        L.gridy++;
        Checkbox cbExpanded = new Checkbox("Image is expanded (apply EF)", expanded);
        left.add(indent(cbExpanded, 14), L);

        if (autoPixelSizeUnitsPerPixel != null) {
            L.gridy++;
            Checkbox cbResetAuto    = new Checkbox("Reset pixel size to auto (" +
                    String.format("%.6g", autoPixelSizeUnitsPerPixel) + " " + autoPixelUnit + "/pixel)", false);
            cbResetAuto.setName("cbResetAuto");
            left.add(indent(cbResetAuto, 14), L);
        }

        //Right column
        BorderedPanel right = new BorderedPanel(new GridBagLayout());
        GridBagConstraints R = gbc(0, 0, 1.0, 0.0);

        //general inputs
        right.add(headingLabel("General Inputs"), R);
        R.gridy++;
        Panel generalInputs = new Panel(new GridBagLayout());
        GridBagConstraints G = gbc(0, 0, 1.0, 0.0);

        //default
        String[] mags = new String[]{"20X", "40X"};
        String autoMag = autoMagnification;
        if ("60X".equals(autoMag) || autoMag == null) autoMag = "40X";
        int idx = ("20X".equals(userMagnification) ? 0 :
                ("40X".equals(userMagnification) ? 1 :
                        ("20X".equals(autoMag) ? 0 : 1)));
        Choice chDefaultMag = new Choice();
        chDefaultMag.add("20X"); chDefaultMag.add("40X");
        chDefaultMag.select(idx);
        addRow(generalInputs, "Default magnification:", chDefaultMag, G);

        //pixel size unit + value
        TextField tfUnit = new TextField(unitShown, 10);
        G.gridy++;
        addRow(generalInputs, "Pixel size unit:", tfUnit, G);

        TextField tfPxSize = new TextField(Double.isNaN(shownPx) ? "" : String.valueOf(shownPx), 10);
        G.gridy++;
        addRow(generalInputs, "Pixel size (" + unitShown + "/pixel):", tfPxSize, G);

        //EF
        TextField tfEF = new TextField(String.valueOf(EF_val), 10);
        G.gridy++;
        addRow(generalInputs, "Expansion factor (EF):", tfEF, G);

        right.add(generalInputs, R);

        R.gridy++;
        right.add(separator(), R);

        // 40X inputs
        R.gridy++;
        right.add(headingLabel("40X - Watershed (Process)"), R);
        R.gridy++;
        Panel p40w = new Panel(new GridBagLayout());
        GridBagConstraints P40 = gbc(0, 0, 1.0, 0.0);

        TextField tfProcMin = new TextField(String.format("%.4f", fieldsInMicrons ? dispProcMin : wsProcMinDist), 8);
        addRow(p40w, fieldsInMicrons ? "Min distance (" + userPixelUnit + "):" : "Min distance (px):", tfProcMin, P40);
        P40.gridy++;
        TextField tfProcThr = new TextField(String.valueOf(wsProcThreshRel), 8);
        addRow(p40w, "Peak threshold (0-1):", tfProcThr, P40);
        P40.gridy++;
        TextField tfProcSig = new TextField(String.valueOf(wsProcSigma), 8);
        addRow(p40w, "Gaussian sigma:", tfProcSig, P40);
        P40.gridy++;
        TextField tfProcMax = new TextField(String.format("%.4f", fieldsInMicrons ? dispProcMax : procMaxPairPx), 8);
        addRow(p40w, fieldsInMicrons ? "Max pair distance (" + userPixelUnit + "):" : "Max pair distance (px):", tfProcMax, P40);

        right.add(p40w, R);

        R.gridy++;
        Panel p40s = new Panel(new GridBagLayout());
        GridBagConstraints P40s = gbc(0, 0, 1.0, 0.0);
        right.add(headingLabel("40X - Size filter (Process)"), R);
        R.gridy++;
        TextField tfKeepLow  = new TextField(String.valueOf(procKeepLow), 6);
        addRow(p40s, "Exclude the lowest (0-1):", tfKeepLow, P40s);
        P40s.gridy++;
        TextField tfKeepHigh = new TextField(String.valueOf(procKeepHigh), 6);
        addRow(p40s, "Exclude the highest (0-1):", tfKeepHigh, P40s);
        right.add(p40s, R);

        R.gridy++;
        right.add(separator(), R);

        // 20X inputs
        R.gridy++;
        right.add(headingLabel("20X - Watershed (RBC)"), R);
        R.gridy++;
        Panel p20r = new Panel(new GridBagLayout());
        GridBagConstraints P20r = gbc(0, 0, 1.0, 0.0);

        TextField tfRbcMin  = new TextField(String.format("%.4f", fieldsInMicrons ? dispRbcMin : wsRbcMinDist), 8);
        addRow(p20r, fieldsInMicrons ? "Min distance (" + userPixelUnit + "):" : "Min distance (px):", tfRbcMin, P20r);
        P20r.gridy++;
        TextField tfRbcThr  = new TextField(String.valueOf(wsRbcThreshRel), 8);
        addRow(p20r, "Peak threshold (0-1):", tfRbcThr, P20r);
        P20r.gridy++;
        TextField tfRbcSig  = new TextField(String.valueOf(wsRbcSigma), 8);
        addRow(p20r, "Gaussian sigma:", tfRbcSig, P20r);
        P20r.gridy++;
        TextField tfRbcLow  = new TextField(String.valueOf(rbcKeepLow), 6);
        addRow(p20r, "Exclude the lowest (0-1):", tfRbcLow, P20r);
        P20r.gridy++;
        TextField tfRbcHigh = new TextField(String.valueOf(rbcKeepHigh), 6);
        addRow(p20r, "Exclude the highest (0-1):", tfRbcHigh, P20r);

        right.add(p20r, R);

        R.gridy++;
        right.add(separator(), R);

        R.gridy++;
        right.add(headingLabel("20X - Watershed (Nuclei)"), R);
        R.gridy++;
        Panel p20n = new Panel(new GridBagLayout());
        GridBagConstraints P20n = gbc(0, 0, 1.0, 0.0);

        TextField tfNucMin  = new TextField(String.format("%.4f", fieldsInMicrons ? dispNucMin : wsNucMinDist), 8);
        addRow(p20n, fieldsInMicrons ? "Min distance (" + userPixelUnit + "):" : "Min distance (px):", tfNucMin, P20n);
        P20n.gridy++;
        TextField tfNucThr  = new TextField(String.valueOf(wsNucThreshRel), 8);
        addRow(p20n, "Peak threshold (0-1):", tfNucThr, P20n);
        P20n.gridy++;
        TextField tfNucSig  = new TextField(String.valueOf(wsNucSigma), 8);
        addRow(p20n, "Gaussian sigma:", tfNucSig, P20n);
        P20n.gridy++;
        TextField tfNucLow  = new TextField(String.valueOf(nucKeepLow), 6);
        addRow(p20n, "Exclude the lowest (0-1):", tfNucLow, P20n);
        P20n.gridy++;
        TextField tfNucHigh = new TextField(String.valueOf(nucKeepHigh), 6);
        addRow(p20n, "Exclude the highest (0-1):", tfNucHigh, P20n);

        right.add(p20n, R);

        //compose two columns
        twoCol.add(left);
        twoCol.add(right);
        gd.addPanel(twoCol);

        // Focus pixel size field if forced
        if (forcePixelFocus) {
            java.awt.EventQueue.invokeLater(() -> {
                tfPxSize.requestFocus();
                tfPxSize.selectAll();
            });
        }

        // Show dialog
        gd.showDialog();
        if (gd.wasCanceled()) return;

        // READBACK
        showLegend           = cbShowLegend.getState();
        showStatsBanner      = cbShowBanner.getState();
        showEnhancedBackground = cbShowEnh.getState();
        showAxisDots         = cbAxisDots.getState();
        showAxisColorBar     = cbAxisBar.getState();
        showProcessLines     = cbProcLines.getState();
        showProcessColorBar  = cbProcBar.getState();
        showProcessContours  = cbProcContours.getState();
        showRbcContours      = cbRbcContours.getState();
        showNucContours      = cbNucContours.getState();
        boolean expandedNew  = cbExpanded.getState();

        boolean doResetToAuto = false;
        for (Component c : left.getComponents()) {
            if (c instanceof Checkbox && "cbResetAuto".equals(c.getName())) {
                doResetToAuto = ((Checkbox)c).getState();
                break;
            }
        }

        userMagnification = chDefaultMag.getSelectedItem();

        String unitIn = tfUnit.getText().trim();
        double pxIn   = parseDoubleSafe(tfPxSize.getText(), Double.NaN);
        double EF_new = parseDoubleSafe(tfEF.getText(), EF_val);
        double efUsed = expandedNew ? EF_new : 1.0;

        double inProcMin = parseDoubleSafe(tfProcMin.getText(), wsProcMinDist);
        wsProcThreshRel  = parseDoubleSafe(tfProcThr.getText(), wsProcThreshRel);
        wsProcSigma      = parseDoubleSafe(tfProcSig.getText(), wsProcSigma);
        double inProcMax = parseDoubleSafe(tfProcMax.getText(), procMaxPairPx);

        procKeepLow  = parseDoubleSafe(tfKeepLow.getText(),  procKeepLow);
        procKeepHigh = parseDoubleSafe(tfKeepHigh.getText(), procKeepHigh);

        double inRbcMin = parseDoubleSafe(tfRbcMin.getText(), wsRbcMinDist);
        wsRbcThreshRel  = parseDoubleSafe(tfRbcThr.getText(), wsRbcThreshRel);
        wsRbcSigma      = parseDoubleSafe(tfRbcSig.getText(), wsRbcSigma);
        rbcKeepLow      = parseDoubleSafe(tfRbcLow.getText(), rbcKeepLow);
        rbcKeepHigh     = parseDoubleSafe(tfRbcHigh.getText(), rbcKeepHigh);

        double inNucMin = parseDoubleSafe(tfNucMin.getText(), wsNucMinDist);
        wsNucThreshRel  = parseDoubleSafe(tfNucThr.getText(), wsNucThreshRel);
        wsNucSigma      = parseDoubleSafe(tfNucSig.getText(), wsNucSigma);
        nucKeepLow      = parseDoubleSafe(tfNucLow.getText(), nucKeepLow);
        nucKeepHigh     = parseDoubleSafe(tfNucHigh.getText(), nucKeepHigh);



        Double pxForConv = userPixelSizeUnitsPerPixel;
        if (doResetToAuto && autoPixelSizeUnitsPerPixel != null) {
            pxForConv = autoPixelSizeUnitsPerPixel;
        } else if (!Double.isNaN(pxIn) && pxIn > 0) {
            pxForConv = pxIn;
        }

        // Convert distances back to pixels
        if (fieldsInMicrons && pxForConv != null && pxForConv > 0) {
            wsProcMinDist = (inProcMin * efUsed) / pxForConv;
            procMaxPairPx = (inProcMax * efUsed) / pxForConv;
            wsRbcMinDist = (inRbcMin * efUsed) / pxForConv;
            wsNucMinDist = (inNucMin * efUsed) / pxForConv;
        } else {
            wsProcMinDist = inProcMin;
            procMaxPairPx = inProcMax;
            wsRbcMinDist = inRbcMin;
            wsNucMinDist = inNucMin;
        }

        expanded = expandedNew;
        EF_val = EF_new;
        if (doResetToAuto && autoPixelSizeUnitsPerPixel != null) {
            userPixelSizeUnitsPerPixel = autoPixelSizeUnitsPerPixel;
            userPixelUnit = autoPixelUnit;
        } else {
            userPixelUnit = unitIn.isEmpty() ? (autoPixelUnit != null ? autoPixelUnit : "units") : unitIn;
            if (!Double.isNaN(pxIn) && pxIn > 0) userPixelSizeUnitsPerPixel = pxIn;
        }



        if (lastResult != null) {
            applyVisibility(lastResult.getOverlay());
            refreshResultFromSettings();
            lastResult.updateAndDraw();
        } else {
            IJ.showStatus("[MagnifySeg] Settings updated.");
        }
    }



    private void applyVisibility(Overlay ov) {
        if (ov == null) return;

        setGroupVisible(ov, legendItems, showLegend);
        setGroupVisible(ov, axisDotItems, showAxisDots);
        setGroupVisible(ov, axisBarItems, showAxisColorBar);
        setGroupVisible(ov, procLineItems, showProcessLines);
        setGroupVisible(ov, procBarItems, showProcessColorBar);
        setGroupVisible(ov, procContourItems,showProcessContours);
        setGroupVisible(ov, rbcContourItems, showRbcContours);
        setGroupVisible(ov, nucContourItems, showNucContours);
        setGroupVisible(ov, statsBannerItems, showStatsBanner);
    }

    private void setGroupVisible(Overlay ov, List<Roi> group, boolean visible) {
        if (ov == null || group == null || group.isEmpty()) return;

        if (!visible) {
            for (int i = ov.size() - 1; i >= 0; i--) {
                Roi r = ov.get(i);
                if (group.contains(r)) ov.remove(i);
            }
            return;
        }

        for (Roi r : group) {
            if (!overlayContains(ov, r)) ov.add(r);
        }
    }

    private static int[] turbo(double t) {
        t = Math.max(0, Math.min(1, t));
        double r = Math.min(1, Math.max(0, 1.5*t - 0.1));
        double g = Math.min(1, Math.max(0, 1.5 - Math.abs(2*t - 1.0)*1.5));
        double b = Math.min(1, Math.max(0, 1.2*(1.0 - t)));
        return new int[]{ (int)(255*r), (int)(255*g), (int)(255*b) };
    }

    private List<Roi> makeColorBar(Overlay ov, int x, int y, int w, int h,
                                   double vmin, double vmax, String label, String cmapName) {
        List<Roi> items = new ArrayList<>();

        Font small = new Font("SansSerif", Font.PLAIN, 25);
        Font labFont = new Font("SansSerif", Font.PLAIN, 30);

        int maxNumW = 0;
        {
            TextRoi m1 = new TextRoi(0, 0, formatVal(vmax), small);
            TextRoi m2 = new TextRoi(0, 0, formatVal(vmin), small);
            maxNumW = (int)Math.round(Math.max(m1.getFloatWidth(), m2.getFloatWidth()));
        }
        int pad = 10;
        int gapNums = 12;
        int panelW = w + gapNums + maxNumW + pad * 2;
        int panelH = h + pad * 2;
        int panelX = x - pad;
        int panelY = y - pad;

        ImageRoi bg = makeBackdrop(panelX, panelY, panelW, panelH, 0.35f);
        ov.add(bg); items.add(bg);

        ColorProcessor cp = new ColorProcessor(w, h);
        for (int yy = 0; yy < h; yy++) {
            double t = 1.0 - (yy / (double)(h - 1));
            int[] rgb = mapColor(cmapName, t);
            int packed = ((rgb[0] & 255) << 16) | ((rgb[1] & 255) << 8) | (rgb[2] & 255);
            for (int xx = 0; xx < w; xx++) cp.set(xx, yy, packed);
        }
        ImageRoi bar = new ImageRoi(x, y, cp);
        ov.add(bar); items.add(bar);

        TextRoi tMax = new TextRoi(x + w + 12, y - 5, formatVal(vmax), small);
        tMax.setStrokeColor(Color.WHITE);
        TextRoi tMin = new TextRoi(x + w + 12, y + h - 16, formatVal(vmin), small);
        tMin.setStrokeColor(Color.WHITE);
        ov.add(tMax); items.add(tMax);
        ov.add(tMin); items.add(tMin);

        // vertical label
        if (label != null && !label.isEmpty()) {
            TextRoi tLab = new TextRoi(0, 0, label, labFont);
            int textLen = (int)Math.round(tLab.getFloatWidth());
            int xLab = x + 62;
            int yLab = y + (h - textLen) / 2;
            tLab.setLocation(xLab, yLab);
            tLab.setAngle(-90.0);
            tLab.setStrokeColor(Color.WHITE);
            ov.add(tLab); items.add(tLab);
        }
        return items;
    }


    private void rebuildStatsBanner(Overlay ov) {
        if (lastResult == null) return;
        if (ov == null) { ov = new Overlay(); lastResult.setOverlay(ov); }

        for (Roi r : statsBannerItems) ov.remove(r);
        statsBannerItems.clear();

        // Compose the text line by magnification
        String sNuc = (lastNucCount != null) ? String.valueOf(lastNucCount) : "N/A";
        String sRbc = (lastRbcCount != null) ? String.valueOf(lastRbcCount) : "N/A";
        String sThk = (lastGbmThicknessUm != null && !Double.isNaN(lastGbmThicknessUm))
                ? (formatVal(lastGbmThicknessUm) + " " + userPixelUnit) : "N/A";
        String sNnd = (lastProcNndUm != null && !Double.isNaN(lastProcNndUm))
                ? (formatVal(lastProcNndUm) + " " + userPixelUnit) : "N/A";

        String line;
        if ("40X".equalsIgnoreCase(String.valueOf(lastEffMag))) {
            line = " GBM thickness: " + sThk + "   |   Process NND: " + sNnd;
        } else if ("20X".equalsIgnoreCase(String.valueOf(lastEffMag))) {
            line = " Nuclei count: " + sNuc + "   |   RBC count: " + sRbc;
        } else {
            line = " Nuclei count: " + sNuc + "   |   RBC count: " + sRbc + "   |   GBM thickness: " + sThk +
                    "   |   Process NND: " + sNnd;
        }

        int w = lastResult.getWidth();
        int h = lastResult.getHeight();

        Font f = new Font("SansSerif", Font.PLAIN, 15);
        TextRoi text = new TextRoi(0, 0, line, f);
        int textH = (int) Math.round(text.getBounds().getHeight());

        int panelW = w;
        int panelH = Math.max(1, textH);
        int x = 0;
        int y = Math.max(0, h - panelH);

        //black background bar
        ImageRoi bg = makeBackdrop(x, y, panelW, panelH, 0.35f);
        bg.setPosition(2);
        ov.add(bg); statsBannerItems.add(bg);


        int tx = x;
        int ty = y;
        text.setLocation(tx, ty);
        text.setStrokeColor(Color.WHITE);
        text.setPosition(2);
        ov.add(text); statsBannerItems.add(text);
    }





    private String formatVal(double v) {
        if (Math.abs(v) >= 100 || Math.abs(v) < 0.01) return String.format("%.2e", v);
        return String.format("%.3f", v);
    }

    // Simple linear interpolation between two RGB colors
    private static int[] lerp(int[] a, int[] b, double t) {
        t = Math.max(0, Math.min(1, t));
        return new int[] {
                (int)Math.round(a[0] + (b[0]-a[0]) * t),
                (int)Math.round(a[1] + (b[1]-a[1]) * t),
                (int)Math.round(a[2] + (b[2]-a[2]) * t)
        };
    }

    private static final int[][] VIRIDIS_ANCHORS = new int[][]{
            {68, 1, 84},
            {59, 82, 139},
            {33, 145, 140},
            {94, 201, 98},
            {253, 231, 37}
    };
    private static int[] viridis(double t) {
        t = Math.max(0, Math.min(1, t));
        double pos = t * (VIRIDIS_ANCHORS.length - 1);
        int i = (int)Math.floor(pos);
        int j = Math.min(i + 1, VIRIDIS_ANCHORS.length - 1);
        double f = pos - i;
        return lerp(VIRIDIS_ANCHORS[i], VIRIDIS_ANCHORS[j], f);
    }

    //black → red → yellow → white
    private static int[] hot(double t) {
        t = Math.max(0, Math.min(1, t));
        double r, g, b;
        if (t < 1.0/3.0) {             // [0, 1/3): ramp up red
            r = 3*t; g = 0; b = 0;
        } else if (t < 2.0/3.0) {      // [1/3, 2/3): red=1, ramp up green
            r = 1; g = 3*t - 1; b = 0;
        } else {                        // [2/3, 1]: red=1, green=1, ramp up blue
            r = 1; g = 1; b = 3*t - 2;
        }
        r = Math.max(0, Math.min(1, r));
        g = Math.max(0, Math.min(1, g));
        b = Math.max(0, Math.min(1, b));
        return new int[]{ (int)(255*r), (int)(255*g), (int)(255*b) };
    }

    private static int[] mapColor(String cmap, double t) {
        switch (cmap == null ? "" : cmap.toLowerCase()) {
            case "viridis": return viridis(t);
            case "hot": return hot(t);
            case "turbo":
            default: return turbo(t);
        }
    }

    //read pixel size /unit from image
    private void refreshAutoPixelSizeFrom(ImagePlus imp) {
        autoPixelSizeUnitsPerPixel = null;
        autoPixelUnit = "µm";

        if (imp == null) return;


        Calibration cal = imp.getCalibration();
        if (cal != null && cal.pixelWidth > 0) {
            autoPixelSizeUnitsPerPixel = cal.pixelWidth;
            if (cal.getUnit() != null && !cal.getUnit().trim().isEmpty()) {
                autoPixelUnit = cal.getUnit().trim();
            }
        } else {
            FileInfo fi = imp.getOriginalFileInfo();
            if (fi != null && fi.pixelWidth > 0) {
                autoPixelSizeUnitsPerPixel = fi.pixelWidth;
                if (fi.unit != null && !fi.unit.trim().isEmpty()) {
                    autoPixelUnit = fi.unit.trim();
                }
            }
        }

        if (userPixelSizeUnitsPerPixel == null && autoPixelSizeUnitsPerPixel != null) {
            userPixelSizeUnitsPerPixel = autoPixelSizeUnitsPerPixel;
            userPixelUnit = autoPixelUnit;
        }
        if (!thresholdsBootstrapped && userPixelSizeUnitsPerPixel != null && userPixelSizeUnitsPerPixel > 0) {
            double ef = expanded ? EF_val : 1.0;
            //to px
            wsProcMinDist = (DEFAULT_MIN_DIST_UM * ef) / userPixelSizeUnitsPerPixel;
            procMaxPairPx = (DEFAULT_MAX_PAIR_UM * ef) / userPixelSizeUnitsPerPixel;
            thresholdsBootstrapped = true;
        }
    }


    private ColorProcessor toRGB(ImageProcessor ip) {
        if (ip == null) {
            return null;
        }
        if (ip instanceof ColorProcessor) {
            return (ColorProcessor) ip.duplicate();
        } else {
            ImageProcessor g = ip.convertToByteProcessor();
            int w = g.getWidth(), h = g.getHeight();
            byte[] p = (byte[]) g.getPixelsCopy();
            ColorProcessor cp = new ColorProcessor(w, h);
            cp.setRGB(p, p, p);
            return cp;
        }
    }

    private Roi makeColoredMaskRoi(File tif, Color color, float alpha, int thicknessPx) throws IOException {
        ImagePlus imp = IJ.openImage(tif.getAbsolutePath());
        if (imp == null) throw new IOException("Open failed: " + tif);
        ImageProcessor ip = imp.getProcessor().convertToByteProcessor();
        int w = ip.getWidth(), h = ip.getHeight();
        byte[] src = (byte[]) ip.getPixels();
        boolean[] mask = new boolean[w*h];
        for (int i = 0; i < mask.length; i++) mask[i] = (src[i] & 0xff) != 0;

        int dilations = Math.max(0, thicknessPx - 1);
        for (int iter = 0; iter < dilations; iter++) {
            boolean[] dst = Arrays.copyOf(mask, mask.length);
            for (int y = 0; y < h; y++) {
                int row = y * w;
                for (int x = 0; x < w; x++) {
                    if (!mask[row + x]) continue;
                    for (int dy = -1; dy <= 1; dy++) {
                        int yy = y + dy;
                        if (yy < 0 || yy >= h) continue;
                        int base = yy * w;
                        for (int dx = -1; dx <= 1; dx++) {
                            int xx = x + dx; if (xx < 0 || xx >= w) continue;
                            dst[base + xx] = true;
                        }
                    }
                }
            }
            mask = dst;
        }

        ColorProcessor cp = new ColorProcessor(w, h);
        int packed = ((color.getRed() & 255) << 16) | ((color.getGreen() & 255) << 8) | (color.getBlue() & 255);
        for (int i = 0; i < w*h; i++) if (mask[i]) cp.set(i % w, i / w, packed);

        ImageRoi roi = new ImageRoi(0, 0, cp);
        roi.setZeroTransparent(true);
        roi.setOpacity(alpha);
        imp.close();
        return roi;
    }


    private String detectMagnification(ImagePlus imp) {
        Object infoObj = imp.getProperty("Info");
        if (infoObj instanceof String) {
            String info = ((String) infoObj).toLowerCase(Locale.ROOT);
            if (info.contains(" 60x") || info.contains(" x60") || info.contains("60x")) return "60X";
            if (info.contains(" 40x") || info.contains(" x40") || info.contains("40x")) return "40X";
            if (info.contains(" 20x") || info.contains(" x20") || info.contains("20x")) return "20X";
        }
        String t = imp.getTitle().toLowerCase(Locale.ROOT);
        if (t.contains("60x")) return "60X";
        if (t.contains("40x")) return "40X";
        if (t.contains("20x")) return "20X";
        return null;
    }

    private String guessTifPath(ImagePlus imp) {
        FileInfo fi = imp.getOriginalFileInfo();
        if (fi != null && fi.directory != null && fi.fileName != null) {
            String path = fi.directory + fi.fileName;
            if (path.toLowerCase(Locale.ROOT).endsWith(".tif") || path.toLowerCase(Locale.ROOT).endsWith(".tiff")) {
                File f = new File(path);
                if (f.isFile()) return f.getAbsolutePath();
            }
        }
        return null;
    }


    private void runCountMenu(String mode, boolean useROI) {
        if (lastSeg20x == null || !lastSeg20x.isFile() || lastTmpDir == null) {
            IJ.showMessage("HistoMagnify",
                    "Run 20X segmentation first (Apply segmentation model, magnification 20X).");
            return;
        }
        if (useROI && (lastResult == null)) {
            IJ.showMessage("HistoMagnify",
                    "No result window found to draw an ROI.\n" +
                            "Re-run segmentation (or use the full-frame command).");
            return;
        }
        int classId = "nuc".equals(mode) ? 1 : 2;

        // ROI
        Roi roi = null;
        if (useROI) {
            if (lastResult.getWindow() == null) lastResult.show();
            lastResult.getWindow().toFront();
            IJ.selectWindow(lastResult.getTitle());
            IJ.setTool("rectangle");
            new ij.gui.WaitForUserDialog("ROI", "Draw a rectangle on the result and click OK.").show();
            roi = lastResult.getRoi();
            if (roi == null || roi.getBounds().width == 0 || roi.getBounds().height == 0) {
                IJ.showMessage("HistoMagnify", "No ROI drawn.");
                return;
            }
        } else {
            if (lastResult != null) lastResult.killRoi();
        }

        File mask = lastTmpDir.resolve(mode+"_mask.tif").toFile();
        try {
            makeBinaryMaskFromLabel(lastSeg20x, roi, classId, mask);
        } catch (IOException e) {
            IJ.handleException(e); return;
        }

        Path venvDir = Paths.get(System.getProperty("user.home")).resolve(".histomagnify").resolve("venv");
        String py = venvDir.resolve(isWin ? "Scripts\\python.exe" : "bin/python").toString();
        File metrics = lastTmpDir.resolve("metrics.py").toFile();
        File outTxt = lastTmpDir.resolve(mode+"_count.txt").toFile();
        File outLabels = lastTmpDir.resolve(mode+"_labels.tif").toFile();
        File outSplitContours = lastTmpDir.resolve(mode+"_contours.tif").toFile();
        File outUnsplitContours = lastTmpDir.resolve(mode+"_outer_contours.tif").toFile();


        try {
            double px = (expanded && EF_val > 0) ? (userPixelSizeUnitsPerPixel / EF_val) : userPixelSizeUnitsPerPixel;
            List<String> args = new ArrayList<>(Arrays.asList(
                    py, metrics.getAbsolutePath(),
                    "--task", mode,
                    "--mask", mask.getAbsolutePath(),
                    "--ws_min_dist", String.valueOf("rbc".equals(mode) ? wsRbcMinDist : wsNucMinDist),
                    "--ws_thresh_rel", String.valueOf("rbc".equals(mode) ? wsRbcThreshRel : wsNucThreshRel),
                    "--ws_sigma", String.valueOf("rbc".equals(mode) ? wsRbcSigma : wsNucSigma),
                    "--keep_low", String.valueOf("rbc".equals(mode) ? rbcKeepLow : nucKeepLow),
                    "--keep_high", String.valueOf("rbc".equals(mode) ? rbcKeepHigh : nucKeepHigh),
                    "--out_txt", outTxt.getAbsolutePath(),
                    "--out_labels", outLabels.getAbsolutePath(),
                    "--out_contours", outSplitContours.getAbsolutePath(), // split
                    "--out_outer_contours", outUnsplitContours.getAbsolutePath()// unsplit
            ));
            runWithLogging(args.toArray(new String[0]), "[metrics:"+mode+"] ", lastTmpDir.toFile());
        } catch (Exception e) {
            IJ.handleException(e);
            return;
        }
        File lbls = lastTmpDir.resolve(mode+"_labels.tif").toFile();
        File splitContours = lastTmpDir.resolve(mode+"_contours.tif").toFile();
        File unsplitContours = lastTmpDir.resolve(mode+"_outer_contours.tif").toFile();

        Overlay ov = lastResult.getOverlay();
        if (ov == null) { ov = new Overlay(); lastResult.setOverlay(ov); }

        List<Roi> bucket = "rbc".equals(mode) ? rbcContourItems : nucContourItems;

        // clear old
        for (Roi r : bucket) ov.remove(r);
        bucket.clear();

        Color col = "rbc".equals(mode) ? Color.WHITE : new Color(255,255,0); // RBC=white, Nuc=yellow
        float alpha = 1.0f;

        java.awt.Rectangle b = (roi != null) ? roi.getBounds() : new java.awt.Rectangle(0,0,0,0);

        try {
            Roi uns = makeColoredMaskRoi(unsplitContours, col, alpha, 2);
            uns.setPosition(2);
            if (roi != null) uns.setLocation(b.x, b.y);
            ov.add(uns); bucket.add(uns);
        } catch (IOException ignore) {
            IJ.log("[MagnifySeg] No unsplit contours for " + mode);
        }
        try {
            Roi spl = makeColoredMaskRoi(splitContours, col, alpha, 2);
            spl.setPosition(2);
            if (roi != null) spl.setLocation(b.x, b.y);
            ov.add(spl); bucket.add(spl);
        } catch (IOException ignore) {
            IJ.log("[MagnifySeg] No split contours for " + mode);
        }

        applyVisibility(ov);
        lastResult.updateAndDraw();

        // read and show
        try {
            String s = new String(java.nio.file.Files.readAllBytes(outTxt.toPath())).trim();
            int cnt = Integer.parseInt(s);
            if ("nuc".equals(mode)) lastNucCount = cnt; else lastRbcCount = cnt;
            rebuildStatsBanner(ov);
            applyVisibility(ov);
            lastResult.updateAndDraw();

            IJ.showMessage("HistoMagnify", ("nuc".equals(mode) ? "Total nuclei: " : "Total RBCs: ") + cnt);
        } catch (Exception ignore) {
            IJ.showMessage("HistoMagnify", "Counting failed.");
        }
    }

    private File exportSinglePageRAW(ImagePlus imp, int page0, Path tmpDir) throws IOException {
        int zSaved = imp.getZ();
        int  cSaved = imp.getC();
        int tSaved = imp.getT();
        try {
            int zPages = Math.max(1, imp.getNSlices());
            int z = Math.min(Math.max(0, page0), zPages - 1);
            imp.setPosition(1, z + 1, tSaved);
            ImageProcessor ip = imp.getProcessor().duplicate();
            File out = tmpDir.resolve("input_page_raw.tif").toFile();
            new ij.io.FileSaver(new ImagePlus("page", ip)).saveAsTiff(out.getAbsolutePath());
            return out;
        } finally {
            imp.setPosition(cSaved, zSaved, tSaved);
        }
    }

    private void refreshResultFromSettings() {
        if (lastResult == null) return;

        try {
            // Choose background source and keep it in RGB
            ColorProcessor bgRGB;
            if (showEnhancedBackground && lastEnhanced != null && lastEnhanced.isFile()) {
                ImagePlus enh = IJ.openImage(lastEnhanced.getAbsolutePath());
                bgRGB = toRGB(enh.getProcessor());
                enh.close();
            } else {
                File raw = (lastTmpDir != null) ? lastTmpDir.resolve("input_page_raw.tif").toFile() : null;
                if (raw == null || !raw.isFile()) {
                    IJ.log("[MagnifySeg] Raw page snapshot not found; re-run segmentation to toggle background.");
                    return;
                }
                ImagePlus rawImp = IJ.openImage(raw.getAbsolutePath());
                bgRGB = toRGB(rawImp.getProcessor());
                rawImp.close();
            }

            // Rebuild the two-slice stack using current LUT and last segmentation
            File segFile = ("20X".equals(lastEffMag) ? lastSeg20x : lastSeg40x);
            if (segFile == null || !segFile.isFile()) return;

            ImagePlus seg = IJ.openImage(segFile.getAbsolutePath());
            ImageProcessor lbl = seg.getProcessor().convertToByteProcessor();
            seg.close();

            int w = bgRGB.getWidth(), h = bgRGB.getHeight();
            ColorProcessor out = (ColorProcessor) bgRGB.duplicate();

            int[][] lut20 = new int[][]{
                    {  0,  0,130},
                    { 91, 24,199},
                    {242, 91, 96},
                    {240,203, 73},
                    { 89,195, 71},
                    { 76, 98,246}
            };
            int[][] lut40 = new int[][]{
                    {105,105,105},
                    {255,  0,255},
                    {  0,255,255}
            };
            int[][] lut = "20X".equals(lastEffMag) ? lut20 : lut40;

            float alpha = 0.45f;
            for (int y=0; y<h; y++) for (int x=0; x<w; x++) {
                int lab = lbl.get(x,y) & 0xff;
                if (lab <= 0 || lab >= lut.length) continue;
                int[] c = lut[lab];
                int p = out.get(x,y);
                int r=(p>>16)&255, g=(p>>8)&255, b0=p&255;
                int nr = (int)(r*(1-alpha) + c[0]*alpha);
                int ng = (int)(g*(1-alpha) + c[1]*alpha);
                int nb = (int)(b0*(1-alpha) + c[2]*alpha);
                out.set(x,y, ((nr&255)<<16)|((ng&255)<<8)|(nb&255));
            }

            ij.ImageStack st = new ij.ImageStack(w,h);
            st.addSlice("Background", bgRGB);
            st.addSlice("Overlay", out);

            Overlay ov = lastResult.getOverlay();
            lastResult.setStack(st);
            lastResult.setOverlay(ov);

            rebuildStatsBanner(ov);
            applyVisibility(ov);

            lastResult.setSlice(2);
            lastResult.updateAndDraw();
        } catch (Exception ex) {
            IJ.handleException(ex);
        }
    }


    private static boolean unitLooksLikePixel(String unit) {
        if (unit == null) return true;
        String u = unit.trim().toLowerCase(Locale.ROOT);
        return u.isEmpty() || "pixel".equals(u) || "pixels".equals(u) || "px".equals(u);
    }

    private boolean overlayContains(Overlay ov, Roi target) {
        if (ov == null || target == null) return false;
        for (int i = 0; i < ov.size(); i++) {
            if (ov.get(i) == target) return true;
        }
        return false;
    }

    private ImageRoi makeBackdrop(int x, int y, int w, int h, float opacity) {
        ColorProcessor cp = new ColorProcessor(w, h);
        int packedBlack = 0; // 0x000000
        for (int yy = 0; yy < h; yy++) {
            for (int xx = 0; xx < w; xx++) cp.set(xx, yy, packedBlack);
        }
        ImageRoi bg = new ImageRoi(x, y, cp);
        bg.setZeroTransparent(false);
        bg.setOpacity(opacity);
        return bg;
    }

    double unitsPerPixel(double pxSize, boolean expanded, double EF) {
        return (expanded && EF > 0) ? (pxSize / EF) : pxSize;
    }

    // Thin and rounded border panel to group sections
    private static class BorderedPanel extends Panel {
        BorderedPanel(LayoutManager lm) { super(lm); }
        @Override public void paint(Graphics g) {
            super.paint(g);
            g.setColor(new Color(160,160,160));
            g.drawRoundRect(0, 0, getWidth()-1, getHeight()-1, 10, 10);
        }
    }

    private static GridBagConstraints gbc(int x, int y, double wx, double wy) {
        GridBagConstraints c = new GridBagConstraints();
        c.gridx = x; c.gridy = y;
        c.weightx = wx; c.weighty = wy;
        c.anchor = GridBagConstraints.NORTHWEST;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.insets = new Insets(4, 8, 4, 8);
        return c;
    }

    private static Label headingLabel(String text) {
        Label l = new Label(text);
        l.setFont(new Font("SansSerif", Font.BOLD, 13));
        return l;
    }

    //Thin separator line
    private static Panel separator() {
        Panel p = new Panel();
        p.setPreferredSize(new Dimension(10, 1));
        p.setBackground(new Color(120,120,120));
        return p;
    }

    private static Panel indent(Component c, int leftPad) {
        Panel w = new Panel(new GridBagLayout());
        GridBagConstraints gc = new GridBagConstraints();
        gc.gridx = 0; gc.gridy = 0; gc.weightx = 1.0; gc.fill = GridBagConstraints.HORIZONTAL;
        gc.insets = new Insets(0, leftPad, 0, 0);
        w.add(c, gc);
        return w;
    }

    private static void addRow(Panel p, String label, Component field, GridBagConstraints base) {
        GridBagConstraints lc = (GridBagConstraints) base.clone();
        lc.gridx = 0; lc.weightx = 0.0; lc.fill = GridBagConstraints.NONE; lc.anchor = GridBagConstraints.WEST;
        Label lab = new Label(label);
        p.add(lab, lc);

        GridBagConstraints fc = (GridBagConstraints) base.clone();
        fc.gridx = 1; fc.weightx = 1.0; fc.fill = GridBagConstraints.HORIZONTAL; fc.anchor = GridBagConstraints.EAST;
        p.add(field, fc);
    }


    private static double parseDoubleSafe(String s, double defVal) {
        try {
            return Double.parseDouble(s.trim());
        } catch (Exception e) {
            return defVal;
        }
    }

    private int runAndGetExit(String[] cmd, String prefix, File workDir)
            throws IOException, InterruptedException {
        Process p = Runtime.getRuntime().exec(cmd, null, workDir);
        Thread outG = (prefix == null) ? null : new StreamGobbler(p.getInputStream(), prefix);
        Thread errG = (prefix == null) ? null : new StreamGobbler(p.getErrorStream(), prefix);
        if (outG != null) outG.start();
        if (errG != null) errG.start();
        int code = p.waitFor();
        if (outG != null) outG.join();
        if (errG != null) errG.join();
        return code;
    }


    private boolean commandWorks(String... probeCmd) {
        try {
            return runAndGetExit(probeCmd, null, null) == 0;
        } catch (Exception e) {
            return false;
        }
    }


    private boolean pythonOK(String... launcher) {
        //must be 3.8+ and have the venv module
        String versionCheck = "import sys,importlib.util; " +
                "sys.exit(0 if sys.version_info[:2]>=(3,8) and importlib.util.find_spec('venv') else 1)";
        try {
            return runAndGetExit(cat(launcher, "-c", versionCheck), null, null) == 0;
        } catch (Exception e) {
            return false;
        }
    }

    private String[] choosePythonLauncher() {
        if (isWin) {
            if (pythonOK("py","-3"))          return new String[]{"py","-3"};
            if (pythonOK("py"))               return new String[]{"py"};
            if (pythonOK("python"))           return new String[]{"python"};
            if (pythonOK("python3"))          return new String[]{"python3"};
        } else {
            if (pythonOK("python3"))          return new String[]{"python3"};
            if (pythonOK("python"))           return new String[]{"python"};
        }
        return null;
    }


    private static String[] cat(String[] a, String... b) {
        String[] out = Arrays.copyOf(a, a.length + b.length);
        System.arraycopy(b, 0, out, a.length, b.length);
        return out;
    }



}
