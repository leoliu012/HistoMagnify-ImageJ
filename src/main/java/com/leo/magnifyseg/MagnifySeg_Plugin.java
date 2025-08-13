package com.leo.magnifyseg;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.io.FileInfo;
import ij.plugin.PlugIn;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

import java.io.*;
import java.net.URL;
import java.net.URLDecoder;
import java.nio.file.*;
import java.util.*;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

import ij.gui.Overlay;
import ij.gui.Roi;
import java.awt.Color;
import java.awt.Font;

import ij.gui.ImageRoi;
import ij.gui.TextRoi;
import ij.measure.Calibration;

public class MagnifySeg_Plugin implements PlugIn {
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


    private static double procMaxPairPx = 560.0;


    private static double wsProcMinDist   = 150;
    private static final double DEFAULT_MIN_DIST_UM  = 0.08;
    private static final double DEFAULT_MAX_PAIR_UM  = 1.5;
    private static boolean thresholdsBootstrapped = false;
    private static double wsProcThreshRel = 0.26;  // 0..1
    private static double wsProcSigma     = 0.0;   // Gaussian grad sigma


    private static Double userPixelSizeUnitsPerPixel = null;
    private static String userPixelUnit = "µm"; // display label only

    private static Double autoPixelSizeUnitsPerPixel = null;
    private static String autoPixelUnit = "µm";

    private static double EF_val = 7.0;
    private static boolean expanded = true;


    private static List<Roi> legendItems = new ArrayList<>();
    private static List<Roi> axisDotItems = new ArrayList<>();
    private static List<Roi> axisBarItems = new ArrayList<>();
    private static List<Roi> procLineItems = new ArrayList<>();
    private static List<Roi> procBarItems = new ArrayList<>();
    private static List<Roi> procContourItems = new ArrayList<>();


    private static final Map<String,String> MODEL_RESOURCE;
    static {
        Map<String,String> m = new HashMap<>();
        m.put("ACTN4",      "/models/ACTN4.hdf5");
        m.put("DAPI",       "/models/DAPI.hdf5");
        m.put("NHS_SINGLE_CHANNEL", "/models/NHS_ester_single.hdf5");
        m.put("NHS_COMBINED_ACTN4",    "/models/NHS_ester_com.hdf5");
        MODEL_RESOURCE = Collections.unmodifiableMap(m);
    }

    @Override
    public void run(String arg) {
        if ("settings".equalsIgnoreCase(arg)) {
            openSettingsDialog();
            return;
        }

        // Secondary menu entries
        if ("thickness".equalsIgnoreCase(arg)) {
            runThicknessMenu();  // will guard NHS was run
            return;
        }
        if ("proc".equalsIgnoreCase(arg) || "process".equalsIgnoreCase(arg)) {
            runProcessMenu();    // will guard ACTN4 was run
            return;
        }

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


        if ("thickness".equalsIgnoreCase(arg)) { runThicknessMenu(true); return; }
        if ("proc".equalsIgnoreCase(arg) || "process".equalsIgnoreCase(arg)) { runProcessMenu(true); return; }

        //list open images
        int count = WindowManager.getImageCount();
        if (count == 0) {
            IJ.showMessage("MagnifySeg", "No images open.");
            return;
        }
        String[] titles = new String[count];
        for (int i = 0; i < count; i++)
            titles[i] = WindowManager.getImage(i+1).getTitle();

        //channel choices for the dialog from all open images
        int maxChannels = 1;
        for (int i = 0; i < count; i++) {
            ImagePlus imp = WindowManager.getImage(i+1);
            if (imp != null) {
                maxChannels = Math.max(maxChannels, Math.max(1, imp.getNChannels()));
            }
        }
        String[] chanOpts = new String[maxChannels + 1];
        chanOpts[0] = "(None)";
        for (int i = 1; i <= maxChannels; i++) chanOpts[i] = "Channel " + i;


        //stain→channel assignment above model selection
        GenericDialog gd = new GenericDialog("MagnifySeg Parameters");
        gd.addChoice("Source image:", titles, titles[0]);
        gd.addStringField("ND2 file path (override):", "");
        gd.addNumericField("Z-slice (0-based):", 0, 0);

        gd.addChoice("ACTN4 channel:",     chanOpts, chanOpts[3]);
        gd.addChoice("DAPI channel:",      chanOpts, chanOpts[1]);
        gd.addChoice("NHS Ester channel:", chanOpts, chanOpts[2]);


        gd.addCheckbox("Use single-channel model for NHS", false);

        gd.addMessage("Models to run:");
        gd.addCheckbox("Run ACTN4", false);
        gd.addCheckbox("Run DAPI",  false);
        gd.addCheckbox("Run NHS",   false);

        gd.addMessage("Background to show:");
        gd.addCheckbox("Use ACTN4", false);
        gd.addCheckbox("Use DAPI",  false);
        gd.addCheckbox("Use NHS Ester (default)", true);

        gd.showDialog();
        if (gd.wasCanceled()) return;


        String sourceTitle = gd.getNextChoice();
        String nd2Override = gd.getNextString().trim();
        int z = (int) gd.getNextNumber();


        String actn4Choice = gd.getNextChoice();
        String dapiChoice  = gd.getNextChoice();
        String nhsChoice   = gd.getNextChoice();

        boolean nhsUseSingle = gd.getNextBoolean();

        // models to run
        boolean runACTN4 = gd.getNextBoolean();
        boolean runDAPI  = gd.getNextBoolean();
        boolean runNHS   = gd.getNextBoolean();

        boolean bgACTN4 = gd.getNextBoolean();
        boolean bgDAPI  = gd.getNextBoolean();
        boolean bgNHS   = gd.getNextBoolean();

        ImagePlus original = WindowManager.getImage(sourceTitle);
        if (original == null) {
            IJ.showMessage("MagnifySeg","Can't find image: "+sourceTitle);
            return;
        }
        int nChannels = Math.max(1, original.getNChannels());
        refreshAutoPixelSizeFrom(original);

        int actn4Chan = parseChannelChoice(actn4Choice);
        int dapiChan  = parseChannelChoice(dapiChoice);
        int nhsChan   = parseChannelChoice(nhsChoice);

        // validation
        if (runACTN4 && actn4Chan < 0) { IJ.showMessage("MagnifySeg","Select a channel for ACTN4 or untick it."); return; }
        if (runDAPI && dapiChan < 0) { IJ.showMessage("MagnifySeg","Select a channel for DAPI or untick it.");  return; }
        if (runNHS && nhsChan < 0) { IJ.showMessage("MagnifySeg","Select a channel for NHS or untick it.");   return; }
        // Background validations: if a box is checked, its channel must be assigned
        if (bgACTN4 && actn4Chan < 0) { IJ.showMessage("MagnifySeg","'Use ACTN4' background is checked, but ACTN4 channel isn't assigned."); return; }
        if (bgDAPI && dapiChan < 0) { IJ.showMessage("MagnifySeg","'Use DAPI' background is checked, but DAPI channel isn't assigned.");   return; }
        if (bgNHS && nhsChan < 0) { IJ.showMessage("MagnifySeg","'Use NHS Ester' background is checked, but NHS channel isn't assigned."); return; }


        // Ensure chosen channel indices are in range
        if ((actn4Chan >= nChannels) || (dapiChan >= nChannels) || (nhsChan >= nChannels)) {
            IJ.showMessage("MagnifySeg","Chosen channel exceeds the channel count of the selected source image.");
            return;
        }

        // Prevent duplicate assignment of the same channel to multiple stains
        Set<Integer> used = new HashSet<>();
        if (actn4Chan >= 0 && !used.add(actn4Chan)) { IJ.showMessage("MagnifySeg","Duplicate channel assignment detected."); return; }
        if (dapiChan  >= 0 && !used.add(dapiChan))  { IJ.showMessage("MagnifySeg","Duplicate channel assignment detected."); return; }
        if (nhsChan   >= 0 && !used.add(nhsChan))   { IJ.showMessage("MagnifySeg","Duplicate channel assignment detected."); return; }



        //determine ND2 path
        String nd2Path = null;
        FileInfo fi = original.getOriginalFileInfo();
        if (fi!=null && fi.directory!=null && fi.fileName!=null)
            nd2Path = fi.directory + fi.fileName;
        if ((nd2Path==null || !new File(nd2Path).isFile()) && !nd2Override.isEmpty())
            nd2Path = nd2Override;
        if (nd2Path==null || !new File(nd2Path).isFile()) {
            IJ.showMessage("MagnifySeg","Cannot locate ND2 file.");
            return;
        }
        File nd2File = new File(nd2Path);

        try {
            //persistent venv (~/.magnifyseg/venv)
            Path home = Paths.get(System.getProperty("user.home"));
            Path baseCache = home.resolve(".magnifyseg");
            Path venvDir = baseCache.resolve("venv");
            if (!Files.exists(venvDir)) {
                IJ.log("[MagnifySeg] Creating venv...");
                Files.createDirectories(baseCache);
                runWithLogging(
                        new String[]{"python3","-m","venv",venvDir.toString()},
                        "[venv] ", null);
                IJ.log("[MagnifySeg] Installing deps...");
                String pipExe = venvDir.resolve("bin/pip").toString();
                Path reqTemp  = baseCache.resolve("requirements.txt");
                try (InputStream in = getClass().getResourceAsStream(REQ_RESOURCE);
                     OutputStream os = Files.newOutputStream(reqTemp)) {
                    byte[] buf = new byte[8192]; int r;
                    while ((r=in.read(buf))>0) os.write(buf,0,r);
                }
                runWithLogging(
                        new String[]{pipExe,"install","-r",reqTemp.toString()},
                        "[pip] ", null
                );
                IJ.log("[MagnifySeg] Venv ready at "+venvDir);
            } else {
                IJ.log("[MagnifySeg] Reusing venv at "+venvDir);
            }

            //per run temp dir
            Path tmpDir = Files.createTempDirectory("magnifyseg_run_");

            //extract chosen plane
            extractFolder(tmpDir, SCRIPTS_ROOT);
            extractFolder(tmpDir, SRC_ROOT);

            String pyExe = venvDir.resolve("bin/python").toString();
            File script  = tmpDir.resolve("segment.py").toFile();

            // Build model list to run
            java.util.List<String> modelsToRun = new ArrayList<>();
            if (runACTN4) modelsToRun.add("ACTN4");
            if (runDAPI) modelsToRun.add("DAPI");
            if (runNHS) {
                if (actn4Chan >= 0 && !nhsUseSingle) {
                    modelsToRun.add("NHS_COMBINED_ACTN4");
                } else {
                    modelsToRun.add("NHS_SINGLE_CHANNEL");
                }
            }

            // Re-check
            if (modelsToRun.contains("ACTN4") && actn4Chan < 0) { IJ.showMessage("MagnifySeg","Select a channel for ACTN4 or untick it."); return; }
            if (modelsToRun.contains("DAPI")  && dapiChan  < 0) { IJ.showMessage("MagnifySeg","Select a channel for DAPI or untick it.");  return; }
            if (modelsToRun.contains("NHS_SINGLE_CHANNEL") && nhsChan < 0) { IJ.showMessage("MagnifySeg","Select a channel for NHS or untick it."); return; }
            if (modelsToRun.contains("NHS_COMBINED_ACTN4") && (nhsChan < 0 || actn4Chan < 0)) {
                IJ.showMessage("MagnifySeg","Combined NHS+ACTN4 model needs both NHS and ACTN4 channels assigned.");
                return;
            }

            if (modelsToRun.isEmpty()) { IJ.showMessage("MagnifySeg","No models selected to run."); return; }

            // Ensure weights are present
            for (String mk : modelsToRun) {
                String resPath = MODEL_RESOURCE.get(mk);
                String fname = Paths.get(resPath).getFileName().toString();
                File weights = extractResource(tmpDir, resPath, fname);
                if (!weights.isFile()) {
                    IJ.showMessage("MagnifySeg","Weights not found for "+mk);
                    return;
                }
            }

            //Run each on its assigned channel
            Map<String, File> segFiles = new LinkedHashMap<>();
            for (String mk : modelsToRun) {
                File out = tmpDir.resolve("seg_" + mk + ".tif").toFile();

                List<String> cmd = new ArrayList<>();
                cmd.add(pyExe);
                cmd.add(script.getAbsolutePath());
                cmd.add("--nd2"); cmd.add(nd2File.getAbsolutePath());
                cmd.add("--z"); cmd.add(String.valueOf(z));
                cmd.add("--modeldir"); cmd.add(tmpDir.toString());
                cmd.add("--model"); cmd.add(mk);
                cmd.add("--output"); cmd.add(out.getAbsolutePath());

                if (mk.equals("ACTN4")) {
                    cmd.add("--channel"); cmd.add(String.valueOf(actn4Chan));
                } else if (mk.equals("DAPI")) {
                    cmd.add("--channel"); cmd.add(String.valueOf(dapiChan));
                } else if (mk.equals("NHS_SINGLE_CHANNEL")) {
                    cmd.add("--channel"); cmd.add(String.valueOf(nhsChan));
                } else if (mk.equals("NHS_COMBINED_ACTN4")) {
                    // Pass NHS as channel A and ACTN4 as channel B
                    cmd.add("--channel"); cmd.add(String.valueOf(nhsChan));
                    cmd.add("--channel2"); cmd.add(String.valueOf(actn4Chan));
                }

                IJ.log("[MagnifySeg] Running segmentation ("+mk+")...");
                runWithLogging(cmd.toArray(new String[0]), "[seg:"+mk+"] ", tmpDir.toFile());
                if (!out.isFile()) { IJ.showMessage("MagnifySeg","Segmentation failed for "+mk); return; }
                segFiles.put(mk, out);
            }


            //RGB background from assigned stains: R=NHS, G=ACTN4, B=DAPI
            int w = original.getWidth(), h = original.getHeight();
            byte[] R = new byte[w*h], G = new byte[w*h], B = new byte[w*h];


            int bgCount = 0;
            bgCount += (bgNHS   && nhsChan   >= 0) ? 1 : 0;
            bgCount += (bgACTN4 && actn4Chan >= 0) ? 1 : 0;
            bgCount += (bgDAPI  && dapiChan  >= 0) ? 1 : 0;

            //none selected -> default to NHS  else first available stain
            if (bgCount == 0) {
                if (nhsChan >= 0) { bgNHS = true; bgCount = 1; }
                else if (actn4Chan >= 0) { bgACTN4 = true; bgCount = 1; }
                else if (dapiChan  >= 0) { bgDAPI  = true; bgCount = 1; }
            }


            // Helper to load once
            final double P_LOW  = 1.0;
            final double P_HIGH = 99.7;
            byte[] pACT = (actn4Chan >= 0) ? enhanced8bitPlane(original, actn4Chan, z, P_LOW, P_HIGH) : null;
            byte[] pDAP = (dapiChan  >= 0) ? enhanced8bitPlane(original, dapiChan,  z, P_LOW, P_HIGH) : null;
            byte[] pNHS = (nhsChan   >= 0) ? enhanced8bitPlane(original, nhsChan,   z, P_LOW, P_HIGH) : null;


            if (bgCount == 1) {
                // Single background -> grayscale
                byte[] src = bgNHS ? pNHS : bgACTN4 ? pACT : pDAP;
                if (src == null) src = new byte[w*h]; // fallback to black if somehow null
                R = src; G = src; B = src;
            } else {

                if (bgNHS && pNHS != null) R = pNHS;               // NHS -> Red
                if (bgACTN4 && pACT != null) G = pACT;                 // ACTN4 -> Green
                if (bgDAPI && pDAP != null) B = pDAP;          // DAPI -> Blue
            }



            ColorProcessor bgRGB = new ColorProcessor(w, h);
            bgRGB.setRGB(R, G, B);
            ColorProcessor bgOverlay = (ColorProcessor) bgRGB.duplicate();

            // Overlay colors
            Map<String,int[]> overlayColors = new HashMap<>();
            overlayColors.put("ACTN4", new int[]{0,255,0});              // green
            overlayColors.put("DAPI",  new int[]{0,0,255});              // blue
            overlayColors.put("NHS_SINGLE_CHANNEL", new int[]{255,0,255}); // magenta
            overlayColors.put("NHS_COMBINED_ACTN4", new int[]{255,0,255}); // magenta


            float alpha = 0.40f;
            for (Map.Entry<String,File> e : segFiles.entrySet()) {
                String mk = e.getKey();
                ImagePlus segImp = IJ.openImage(e.getValue().getAbsolutePath());
                if (segImp == null) { IJ.showMessage("MagnifySeg", "Cannot open "+e.getValue()); return; }
                ImageProcessor lbl = segImp.getProcessor().convertToByteProcessor();

                if (mk.equals("NHS_SINGLE_CHANNEL") || mk.equals("NHS_COMBINED_ACTN4")) {
                    // NHS special rules:
                    // - label 1 => magenta
                    // - label 2 => red (unless Run DAPI is checked, then skip to avoid duplicate nuclei)
                    final int[] MAGENTA = new int[]{255,0,255};
                    final int[] RED     = new int[]{255,0,0};

                    for (int ypx = 0; ypx < h; ypx++) {
                        for (int xpx = 0; xpx < w; xpx++) {
                            int lab = lbl.getPixel(xpx, ypx);
                            int[] col;
                            if (lab == 1) {
                                col = MAGENTA;
                            } else if (lab == 2) {
                                if (runDAPI) continue; // hide NHS nuclei when DAPI model is used
                                col = RED;
                            } else {
                                continue; // ignore background / other labels
                            }

                            int pix0 = bgOverlay.getPixel(xpx, ypx);
                            int br = (pix0>>16)&0xff, bgc=(pix0>>8)&0xff, bb=pix0&0xff;
                            int nr = (int)(br * (1-alpha) + col[0] * alpha);
                            int ng = (int)(bgc * (1-alpha) + col[1] * alpha);
                            int nb = (int)(bb * (1-alpha) + col[2] * alpha);
                            int blended = ((nr&0xff)<<16)|((ng&0xff)<<8)|(nb&0xff);
                            bgOverlay.set(xpx, ypx, blended);
                        }
                    }
                } else {
                    int[] col = overlayColors.getOrDefault(mk, new int[]{255,255,0});
                    for (int ypx = 0; ypx < h; ypx++) {
                        for (int xpx = 0; xpx < w; xpx++) {
                            int lab = lbl.getPixel(xpx, ypx);
                            if (lab > 0) {
                                int pix0 = bgOverlay.getPixel(xpx, ypx);
                                int br  = (pix0>>16)&0xff, bgc=(pix0>>8)&0xff, bb=pix0&0xff;
                                int nr  = (int)(br  * (1-alpha) + col[0] * alpha);
                                int ng  = (int)(bgc * (1-alpha) + col[1] * alpha);
                                int nb  = (int)(bb  * (1-alpha) + col[2] * alpha);
                                int blended = ((nr&0xff)<<16)|((ng&0xff)<<8)|(nb&0xff);
                                bgOverlay.set(xpx, ypx, blended);
                            }
                        }
                    }
                }
            }

            // 2-slice:Background, Overlay]
            ij.ImageStack stack = new ij.ImageStack(w, h);
            stack.addSlice("Background", bgRGB);
            stack.addSlice("Overlay", bgOverlay);
            ImagePlus result = new ImagePlus("Result", stack);

            // Legend for only the models we ran
            Overlay ov = new Overlay();
            int x0 = 25, y0 = 25, box = 40;
            Font f = new Font("SansSerif", Font.PLAIN, 40);
            int yy = y0;
            for (String mk : modelsToRun) {
                int[] col = overlayColors.get(mk);
                Roi r = new Roi(x0, yy, box, box);
                Color cc = new Color(col[0], col[1], col[2]);
                r.setFillColor(cc);
                r.setStrokeColor(cc);
                ov.add(r);


                String label = mk.equals("NHS_SINGLE_CHANNEL") ? "GBM"
                        : mk.equals("NHS_COMBINED_ACTN4") ? "GBM (Multi Channel)"
                        : mk.equals("DAPI") ? "Nuclei"
                        : mk.equals("ACTN4") ? "Process"
                        : mk;


                ij.gui.TextRoi t = new ij.gui.TextRoi(x0+box+20, yy, label, f);
                t.setStrokeColor(Color.WHITE);
                ov.add(t);
                yy += box + 20;
            }
            legendItems.clear();
            for (int i = 0; i < ov.size(); i++) legendItems.add(ov.get(i));

            axisDotItems.clear();
            axisBarItems.clear();
            procLineItems.clear();
            procBarItems.clear();
            procContourItems.clear();


            applyVisibility(ov); // respect current settings

            result.setOverlay(ov);
            result.show();

            MagnifySeg_Plugin.lastTmpDir = tmpDir;
            MagnifySeg_Plugin.lastSegFiles = new LinkedHashMap<>(segFiles);
            MagnifySeg_Plugin.lastModelsRan = new HashSet<>(modelsToRun);
            MagnifySeg_Plugin.lastResult = result;
            return;


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

    private File extractResource(Path tmpDir, String resPath, String name)
            throws IOException {
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


    private int parseChannelChoice(String choice) {
        if (choice == null || choice.startsWith("(None)")) return -1;
        try {
            String[] parts = choice.trim().split("\\s+"); // "Channel N"
            return Integer.parseInt(parts[1]) - 1;        // 0-based
        } catch (Exception e) {
            return -1;
        }
    }

    // Normalize by max -> scale to [0,255] -> percentile stretch (pLow..pHigh) -> 8-bit
    private byte[] enhanced8bitPlane(ImagePlus imp, int c0 , int z0,
                                     double pLow, double pHigh) {
        int cSaved = imp.getC(), zSaved = imp.getZ(), tSaved = imp.getT();
        try {
            imp.setPosition(c0+1, z0+1, tSaved);
            ImageProcessor ip = imp.getProcessor().convertToFloatProcessor(); // raw float copy
            int w = ip.getWidth(), h = ip.getHeight(), n = w*h;
            float[] f = (float[]) ip.getPixels();

            // normalize by max
            float max = 0f;
            for (int i = 0; i < n; i++) if (f[i] > max) max = f[i];
            // build 8-bit pre-enhance
            int[] hist = new int[256];
            byte[] src8 = new byte[n];
            if (max > 0f) {
                for (int i = 0; i < n; i++) {
                    int v = (int) Math.round((f[i] / max) * 255.0);
                    if (v < 0) v = 0; else if (v > 255) v = 255;
                    src8[i] = (byte) v;
                    hist[v]++;
                }
            } else {
                // all zeros
                return new byte[n];
            }

            // percentile thresholds on 8-bit histogram
            int lo = percentileFromHist(hist, pLow, n);
            int hi = percentileFromHist(hist, pHigh, n);
            if (hi <= lo) hi = lo + 1; // avoid div by zero

            // linear stretch to [0,255]
            byte[] out = new byte[n];
            for (int i = 0; i < n; i++) {
                int v = src8[i] & 0xff;
                if (v <= lo) out[i] = 0;
                else if (v >= hi) out[i] = (byte) 255;
                else {
                    int vv = (int) Math.round((v - lo) * 255.0 / (hi - lo));
                    out[i] = (byte) (vv < 0 ? 0 : vv > 255 ? 255 : vv);
                }
            }
            return out;
        } finally {
            imp.setPosition(cSaved, zSaved, tSaved);
        }
    }

    private int percentileFromHist(int[] hist, double p, int total) {
        if (total <= 0) return 0;
        long target = Math.round((p / 100.0) * (total - 1));
        long cum = 0;
        for (int i = 0; i < 256; i++) {
            cum += hist[i];
            if (cum > target) return i;
        }
        return 255;
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

    private void makeBinaryMaskFromPositive(File labelTif, Roi roi, File outTif) throws IOException {
        ImagePlus imp = IJ.openImage(labelTif.getAbsolutePath());
        if (imp == null) throw new IOException("Open failed: "+labelTif);
        if (roi != null) imp.setRoi(roi);
        ImageProcessor ip = (roi != null ? imp.getProcessor().crop() : imp.getProcessor()).convertToByteProcessor();
        int w = ip.getWidth(), h = ip.getHeight();
        byte[] pix = (byte[]) ip.getPixels();
        for (int i=0;i<w*h;i++) {
            int v = pix[i] & 0xff;
            pix[i] = (byte)((v>0)?255:0);
        }
        new ij.io.FileSaver(new ImagePlus("mask", ip)).saveAsTiff(outTif.getAbsolutePath());
        imp.close();
    }
    private void runThicknessMenu() { runThicknessMenu(true); }
    private void runProcessMenu()   { runProcessMenu(true);  }

    private void runThicknessMenu(boolean useROI) {
        // guards
        if (lastSegFiles == null || lastResult == null || lastTmpDir == null) {
            IJ.showMessage("MagnifySeg", "Run segmentation first.");
            return;
        }
        boolean haveNHS = lastModelsRan != null && (
                lastModelsRan.contains("NHS_SINGLE_CHANNEL") ||
                        lastModelsRan.contains("NHS_COMBINED_ACTN4")
        );
        if (!haveNHS) {
            IJ.showMessage("MagnifySeg", "GBM thickness requires the NHS model. Re-run segmentation with \"Run NHS\" checked.");
            return;
        }

        //Must have a pixel size. If missing, force user to go to Settings.
        if (userPixelSizeUnitsPerPixel == null || userPixelSizeUnitsPerPixel <= 0) {
            IJ.showMessage("MagnifySeg",
                    "Pixel size is not set.\nOpen Settings and enter a pixel size (or reset to auto) before running thickness.");
            openSettingsDialog();
            if (userPixelSizeUnitsPerPixel == null || userPixelSizeUnitsPerPixel <= 0) return;
        }
        double pxSize = userPixelSizeUnitsPerPixel;

        //ROI- a rectangle.
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
            // full image-ignore any current selection
            if (lastResult != null) lastResult.killRoi();
        }

        // choose NHS seg
        File nhsSeg = lastSegFiles.containsKey("NHS_COMBINED_ACTN4")
                ? lastSegFiles.get("NHS_COMBINED_ACTN4")
                : lastSegFiles.get("NHS_SINGLE_CHANNEL");
        if (nhsSeg == null) {
            IJ.showMessage("MagnifySeg", "No NHS (GBM) segmentation available.");
            return;
        }

        // prepare masks & outputs
        File gbmMask  = lastTmpDir.resolve("gbm_mask.tif").toFile();
        File thickTxt = lastTmpDir.resolve("thickness.txt").toFile();
        File thickCsv = lastTmpDir.resolve("thickness_points.csv").toFile();
        try {
            makeBinaryMaskFromLabel(nhsSeg, roi, 1, gbmMask); // label==1 => GBM
        } catch (IOException ex) {
            IJ.handleException(ex); return;
        }

        Path venvDir = Paths.get(System.getProperty("user.home"))
                .resolve(".magnifyseg").resolve("venv");
        String pyExe = venvDir.resolve("bin/python").toString();
        File metricsPy = lastTmpDir.resolve("metrics.py").toFile();

        // run metrics
        double px = expanded && EF_val > 0 ? (pxSize / EF_val) : pxSize;
        String[] cmdT = new String[]{
                pyExe, metricsPy.getAbsolutePath(),
                "--task","thickness",
                "--mask", gbmMask.getAbsolutePath(),
                "--px", String.valueOf(px),
                "--out_txt", thickTxt.getAbsolutePath(),
                "--out_csv", thickCsv.getAbsolutePath()
        };
        try {
            IJ.log("[MagnifySeg] Computing GBM thickness...");
            runWithLogging(cmdT, "[metrics:thick] ", lastTmpDir.toFile());
            IJ.log("Finished computing GBM thickness");
        } catch (Exception ex) {
            IJ.handleException(ex); return;
        }

        // overlay points + show value
        Overlay ov2 = lastResult.getOverlay();
        if (ov2 == null) {
            ov2 = new Overlay();
            lastResult.setOverlay(ov2);
        }

        // clear previous
        for (Roi r: axisDotItems) ov2.remove(r);
        for (Roi r: axisBarItems) ov2.remove(r);
        axisDotItems.clear();
        axisBarItems.clear();

        double val = Double.NaN;
        List<double[]> pts = new ArrayList<>();
        try {
            try {
                String s = new String(java.nio.file.Files.readAllBytes(thickTxt.toPath())).trim();
                val = Double.parseDouble(s);
            } catch (Exception ignore) {}

            try (BufferedReader br = new BufferedReader(new FileReader(thickCsv))) {
                String line;
                java.awt.Rectangle b = (roi != null) ? roi.getBounds() : new java.awt.Rectangle(0,0,0,0);
                while ((line = br.readLine()) != null) {
                    String[] t = line.split(",");
                    if (t.length >= 2) {
                        double x = Double.parseDouble(t[0]) + b.x;
                        double y = Double.parseDouble(t[1]) + b.y;
                        double tv = (t.length >= 3) ? Double.parseDouble(t[2]) : Double.NaN;
                        pts.add(new double[]{x,y,tv});
                    }
                }
            }
        } catch (Exception ignore) {}

        if (!pts.isEmpty()) {
            // compute color scale
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



            // a filled dot at each skeleton pixel
            for (double[] p : pts) {
                double v = haveVals ? p[2] : 0.5;

                double t = (vmax > vmin) ? (v - vmin) / (vmax - vmin) : 0.0;
                int[] rgb = viridis(t);

                int cx = (int) Math.round(p[0]);
                int cy = (int) Math.round(p[1]);
                ij.gui.OvalRoi dot = new ij.gui.OvalRoi(cx - 1, cy - 1, 4, 4);
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

        applyVisibility(ov2);
        lastResult.updateAndDraw();

        IJ.showMessage("MagnifySeg", "Avg GBM thickness: " + val + " " + userPixelUnit);
    }


    private void runProcessMenu(boolean useROI)  {

        // guards
        if (lastSegFiles == null || lastResult == null || lastTmpDir == null) {
            IJ.showMessage("MagnifySeg", "Run segmentation first.");
            return;
        }
        if (lastModelsRan == null || !lastModelsRan.contains("ACTN4")) {
            IJ.showMessage("MagnifySeg", "Process distance requires the ACTN4 model. Re-run segmentation with \"Run ACTN4\" checked.");
            return;
        }

        // We need a pixel size. If missing, force user to go to Settings.
        if (userPixelSizeUnitsPerPixel == null || userPixelSizeUnitsPerPixel <= 0) {
            IJ.showMessage("MagnifySeg",
                    "Pixel size is not set.\nOpen Settings and enter a pixel size (or reset to auto) before running process distance.");
            openSettingsDialog();
            if (userPixelSizeUnitsPerPixel == null || userPixelSizeUnitsPerPixel <= 0) return;
        }
        double pxSize = userPixelSizeUnitsPerPixel;

        double maxPairPx = procMaxPairPx;

        // ROI per Settings
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


        // ACTN4 seg -> positive mask
        File actn4Seg = lastSegFiles.get("ACTN4");
        if (actn4Seg == null) {
            IJ.showMessage("MagnifySeg", "No ACTN4 segmentation available.");
            return;
        }
        File procMask = lastTmpDir.resolve("proc_mask.tif").toFile();
        File procTxt  = lastTmpDir.resolve("proc.txt").toFile();
        File procCsv  = lastTmpDir.resolve("proc_pairs.csv").toFile();
        File procLabels = lastTmpDir.resolve("proc_labels.tif").toFile();
        File procEdges  = lastTmpDir.resolve("proc_contours.tif").toFile();
        File procOuter  = lastTmpDir.resolve("proc_outer_contours.tif").toFile();


        try {
            makeBinaryMaskFromPositive(actn4Seg, roi, procMask);
        } catch (IOException ex) {
            IJ.handleException(ex); return;
        }

        // venv + metrics.py
        Path venvDir = Paths.get(System.getProperty("user.home"))
                .resolve(".magnifyseg").resolve("venv");
        String pyExe = venvDir.resolve("bin/python").toString();
        File metricsPy = lastTmpDir.resolve("metrics.py").toFile();

        // run metrics
        double px = expanded && EF_val > 0 ? (pxSize / EF_val) : pxSize;

        String[] cmdP = new String[]{
                pyExe, metricsPy.getAbsolutePath(),
                "--task","proc",
                "--mask", procMask.getAbsolutePath(),
                "--px", String.valueOf(px),
                "--max_pair_px", String.valueOf(maxPairPx),

                "--ws_min_dist",   String.valueOf(wsProcMinDist),
                "--ws_thresh_rel", String.valueOf(wsProcThreshRel),
                "--ws_sigma",      String.valueOf(wsProcSigma),

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

        Overlay ov2 = lastResult.getOverlay();
        if (ov2 == null) { ov2 = new Overlay(); lastResult.setOverlay(ov2); }

        for (Roi r: procLineItems)    ov2.remove(r);
        for (Roi r: procBarItems)     ov2.remove(r);
        for (Roi r: procContourItems) ov2.remove(r);
        procLineItems.clear();
        procBarItems.clear();
        procContourItems.clear();

        try {
            java.awt.Rectangle b = (roi != null) ? roi.getBounds() : new java.awt.Rectangle(0,0,0,0);

            //unsplit outer boundary
            Roi edgesOuter = makeColoredMaskRoi(procOuter, new Color(255, 240, 6), 1.0f);
            edgesOuter.setPosition(2);
            if (roi != null) edgesOuter.setLocation(b.x, b.y);
            ov2.add(edgesOuter);
            procContourItems.add(edgesOuter);

            //split watershed boundary
            Roi edgesSplit = makeColoredMaskRoi(procEdges, new Color(255, 240, 6), 1.0f);
            edgesSplit.setPosition(2);
            if (roi != null) edgesSplit.setLocation(b.x, b.y);
            ov2.add(edgesSplit);
            procContourItems.add(edgesSplit);

        } catch (IOException ignore) {
            IJ.log("[MagnifySeg] No process contours produced.");
        }

        //draw NND lines + colorbar
        double val = Double.NaN;
        class Pair { double x0,y0,x1,y1,distPx; }
        List<Pair> pairs = new ArrayList<>();
        try {
            try {
                String s = new String(java.nio.file.Files.readAllBytes(procTxt.toPath())).trim();
                val = Double.parseDouble(s);
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
                ln.setStrokeWidth(2.0);
                ln.setPosition(2);
                ov2.add(ln);
                procLineItems.add(ln);
            }

            List<Roi> bar = makeColorBar(ov2,
                    (int)(lastResult.getWidth()-110), Math.max(120, lastResult.getHeight()/4)+100, 20, Math.max(120, lastResult.getHeight()/4),
                    dmin, dmax, "Proc dist (" + userPixelUnit + ")", "hot");
            procBarItems.addAll(bar);
        }

        applyVisibility(ov2);
        lastResult.updateAndDraw();

        IJ.showMessage("MagnifySeg", "Process mean nearest-neighbor distance: " + val + " " + userPixelUnit);
    }


    private void openSettingsDialog() {
        GenericDialog gd = new GenericDialog("MagnifySeg Settings");
        // existing visibility toggles
        gd.addCheckbox("Show legend", showLegend);
        gd.addCheckbox("Show GBM axis", showAxisDots);
        gd.addCheckbox("Show GBM color bar", showAxisColorBar);
        gd.addCheckbox("Show process lines", showProcessLines);
        gd.addCheckbox("Show process color bar", showProcessColorBar);
        gd.addCheckbox("Show process contours", showProcessContours);

        gd.addMessage("");



        double shownPx = (userPixelSizeUnitsPerPixel != null) ? userPixelSizeUnitsPerPixel : Double.NaN;
        String unitShown = (userPixelUnit != null && !userPixelUnit.trim().isEmpty()) ? userPixelUnit : "units";

        gd.addStringField("Pixel size unit:", unitShown);
        gd.addNumericField("Pixel size (" + unitShown + "/pixel):", shownPx, 6);

        if (autoPixelSizeUnitsPerPixel != null) {
            gd.addCheckbox("Reset pixel size to auto (" +
                    String.format("%.6g", autoPixelSizeUnitsPerPixel) + " " + autoPixelUnit + "/pixel)", false);
        } else {
            gd.addMessage("No pixel size found in metadata; please enter a value.");

        }
        gd.addCheckbox("Image is expanded (apply EF)", expanded);
        gd.addNumericField("Expansion factor (EF):", EF_val, 3);

        gd.addMessage(" ");
        gd.addMessage("Watershed (Process)");

        boolean fieldsInMicrons = !Double.isNaN(shownPx) && !"units".equals(unitShown);

        if (fieldsInMicrons) {
            double ef = expanded ? EF_val : 1.0;
            double min_um = wsProcMinDist * (userPixelSizeUnitsPerPixel / ef);
            double max_um = procMaxPairPx * (userPixelSizeUnitsPerPixel / ef);

            gd.addNumericField("Min distance - If applicable: After EF correction ("+userPixelUnit+")", min_um, 2);
            gd.addNumericField("Peak threshold (0–1):", wsProcThreshRel, 2);
            gd.addNumericField("Gaussian sigma:", wsProcSigma, 2);
            gd.addNumericField("Max pair distance - If applicable: After EF correction ("+userPixelUnit+")", max_um, 1);
        } else {
            gd.addNumericField("Min distance (px):", wsProcMinDist, 2);
            gd.addNumericField("Peak threshold (0–1):", wsProcThreshRel, 2);
            gd.addNumericField("Gaussian sigma:", wsProcSigma, 2);
            gd.addNumericField("Max pair distance (px):", procMaxPairPx, 1);
        }


        gd.showDialog();
        if (gd.wasCanceled()) return;

        // read back toggles
        showLegend = gd.getNextBoolean();
        showAxisDots = gd.getNextBoolean();
        showAxisColorBar = gd.getNextBoolean();
        showProcessLines = gd.getNextBoolean();
        showProcessColorBar= gd.getNextBoolean();
        showProcessContours = gd.getNextBoolean();



        // read pixel size unit + value
        String unitIn = gd.getNextString().trim();
        double pxIn   = gd.getNextNumber();

        boolean doResetToAuto = false;
        if (autoPixelSizeUnitsPerPixel != null) {
            doResetToAuto = gd.getNextBoolean();
        }

        expanded = gd.getNextBoolean();
        EF_val = gd.getNextNumber();

        double ef = expanded ? EF_val : 1.0;

        double inMin = gd.getNextNumber();
        wsProcThreshRel = gd.getNextNumber();
        wsProcSigma  = gd.getNextNumber();
        double inMax = gd.getNextNumber();

        if (fieldsInMicrons) {
            wsProcMinDist = (inMin * ef) / userPixelSizeUnitsPerPixel;
            procMaxPairPx = (inMax * ef) / userPixelSizeUnitsPerPixel;
        } else {
            wsProcMinDist = inMin;
            procMaxPairPx = inMax;
        }


        if (doResetToAuto) {
            userPixelSizeUnitsPerPixel = autoPixelSizeUnitsPerPixel;
            userPixelUnit = autoPixelUnit;
        } else {
            userPixelUnit = unitIn.isEmpty() ? (autoPixelUnit != null ? autoPixelUnit : "units") : unitIn;
            if (Double.isNaN(pxIn) || pxIn <= 0) {
                // leave as null; metrics will block until user sets it
                userPixelSizeUnitsPerPixel = null;
            } else {
                userPixelSizeUnitsPerPixel = pxIn;
            }
        }

        if (lastResult != null) {
            applyVisibility(lastResult.getOverlay());
            lastResult.updateAndDraw();
        } else {
            IJ.showStatus("[MagnifySeg] Settings updated.");
        }
    }

    private void applyVisibility(Overlay ov) {
        if (ov == null) return;

        setGroupVisible(ov,legendItems, showLegend);
        setGroupVisible(ov, axisDotItems, showAxisDots);
        setGroupVisible(ov, axisBarItems, showAxisColorBar);
        setGroupVisible(ov, procLineItems, showProcessLines);
        setGroupVisible(ov, procBarItems, showProcessColorBar);
        setGroupVisible(ov, procContourItems, showProcessContours);

    }

    private void setGroupVisible(Overlay ov, List<Roi> group, boolean visible) {
        if (ov == null || group == null || group.isEmpty()) return;

        for (int i = ov.size() - 1; i >= 0; i--) {
            Roi r = ov.get(i);
            if (group.contains(r)) ov.remove(i);
        }
        if (visible) {
            for (Roi r : group) ov.add(r);
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
        ColorProcessor cp = new ColorProcessor(w, h);
        for (int yy = 0; yy < h; yy++) {
            double t = 1.0 - (yy / (double)(h - 1));
            int[] rgb = mapColor(cmapName, t);
            int packed = ((rgb[0] & 255) << 16) | ((rgb[1] & 255) << 8) | (rgb[2] & 255);
            for (int xx = 0; xx < w; xx++) cp.set(xx, yy, packed);
        }
        ImageRoi bar = new ImageRoi(x, y, cp);
        ov.add(bar); items.add(bar);

        Font small = new Font("SansSerif", Font.PLAIN, 25);
        TextRoi tMax = new TextRoi(x + w + 12, y - 5, formatVal(vmax), small);
        tMax.setStrokeColor(Color.WHITE);
        TextRoi tMin = new TextRoi(x + w + 12, y + h - 16, formatVal(vmin), small);
        tMin.setStrokeColor(Color.WHITE);
        ov.add(tMax); ov.add(tMin);
        items.add(tMax); items.add(tMin);

        if (label != null && !label.isEmpty()) {
            Font labFont = new Font("SansSerif", Font.BOLD, 30);

            TextRoi tLab = new TextRoi(0, 0, label, labFont);

            int textLen = (int) Math.round(tLab.getFloatWidth());

            int xLab = x + 62;
            int yLab = y + (h - textLen) / 2;

            tLab.setLocation(xLab, yLab);
            tLab.setAngle(-90.0);              // vertical (portrait) label
            tLab.setStrokeColor(Color.WHITE);
            ov.add(tLab); items.add(tLab);
        }
        return items;
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

    // Viridis ia 5 anchors from Matplotlib's palette
    private static final int[][] VIRIDIS_ANCHORS = new int[][]{
            {68, 1, 84},    // 0.00
            {59, 82, 139},  // 0.25
            {33, 145, 140}, // 0.50
            {94, 201, 98},  // 0.75
            {253, 231, 37}  // 1.00
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
            case "hot":     return hot(t);
            case "turbo":
            default:        return turbo(t);
        }
    }

    //read pixel size and unit from image's calibration
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
            // Fallback: fileInfo sometimes also carries calibration
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
            wsProcMinDist = (DEFAULT_MIN_DIST_UM * ef) / userPixelSizeUnitsPerPixel; // -> px
            procMaxPairPx = (DEFAULT_MAX_PAIR_UM * ef) / userPixelSizeUnitsPerPixel; // -> px
            thresholdsBootstrapped = true;
        }
    }

    private Roi makeColoredMaskRoi(File tif, Color color, float alpha) throws IOException {
        ImagePlus imp = IJ.openImage(tif.getAbsolutePath());
        if (imp == null) throw new IOException("Open failed: " + tif);
        ImageProcessor ip = imp.getProcessor().convertToByteProcessor();
        int w = ip.getWidth(), h = ip.getHeight();

        ColorProcessor cp = new ColorProcessor(w, h);
        int packed = ((color.getRed() & 255) << 16) | ((color.getGreen() & 255) << 8) | (color.getBlue() & 255);
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if ((ip.get(x,y) & 0xff) != 0) cp.set(x, y, packed);
            }
        }
        ImageRoi roi = new ImageRoi(0, 0, cp);
        roi.setZeroTransparent(true);
        roi.setOpacity(alpha);
        imp.close();
        return roi;
    }


}
