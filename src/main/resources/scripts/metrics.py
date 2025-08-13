#!/usr/bin/env python3
import argparse
import numpy as np
import tifffile as tiff
from scipy import ndimage as ndi
from skimage.morphology import medial_axis
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from skimage.measure import regionprops, label as sklabel


def gbm_thickness(mask_path, px_size, out_txt, out_csv):
    m = tiff.imread(mask_path).astype(bool)
    if m.ndim == 3:  # sometimes as a 1slice stack
        m = m[0]
    if not np.any(m):
        with open(out_txt, "w") as f:
            f.write("nan")
        open(out_csv, "w").close()
        return

    skel, dist = medial_axis(m, return_distance=True)
    diam = 2.0 * dist[skel] * float(px_size)  # units
    mean_val = float(np.mean(diam)) if diam.size else float("nan")

    ys, xs = np.nonzero(skel)
    with open(out_csv, "w") as f:
        for x, y, d in zip(xs, ys, diam):
            f.write(f"{int(x)},{int(y)},{d:.6f}\n")
    with open(out_txt, "w") as f:
        f.write(f"{mean_val:.6f}")


#Process
def ws_split(mask_u8, min_dist, thresh_rel, sigma):
    mask = mask_u8.astype(bool)
    if not np.any(mask):
        return np.zeros_like(mask, dtype=np.uint16)

    dist = ndi.distance_transform_edt(mask)
    peaks = peak_local_max(
        dist,
        labels=mask,
        footprint=np.ones((3, 3)),
        min_distance=max(1, int(round(min_dist))),
        threshold_rel=float(thresh_rel),
        exclude_border=False,
    )
    markers = np.zeros_like(dist, dtype=np.int32)
    for i, (r, c) in enumerate(peaks, start=1):
        markers[r, c] = i

    grad = ndi.gaussian_gradient_magnitude(dist.astype(np.float32), sigma=float(sigma))
    labels = watershed(grad, markers=markers, mask=mask)
    return labels.astype(np.uint16)


def labels_to_contours(labels_u16):
    lab = labels_u16.astype(np.int32)
    h, w = lab.shape
    edges = np.zeros((h, w), np.uint8)

    # mark label boundaries in 4-neighborhood
    edges[1:, :]  |= ((lab[1:, :]  != 0) & (lab[1:, :]  != lab[:-1, :])).astype(np.uint8)
    edges[:-1, :] |= ((lab[:-1, :] != 0) & (lab[:-1, :] != lab[1:, :])).astype(np.uint8)
    edges[:, 1:]  |= ((lab[:, 1:]  != 0) & (lab[:, 1:]  != lab[:, :-1])).astype(np.uint8)
    edges[:, :-1] |= ((lab[:, :-1] != 0) & (lab[:, :-1] != lab[:, 1:])).astype(np.uint8)

    return (edges * 255).astype(np.uint8)


def proc_ws_nnd(
    mask_path,
    px_size,
    max_pair_px,
    ws_min_dist,
    ws_thresh_rel,
    ws_sigma,
    out_txt,
    out_csv,
    out_labels=None,
    out_contours=None,
    out_outer_contours=None,
):
    mask = tiff.imread(mask_path).astype(np.uint8)
    if mask.ndim == 3:
        mask = mask[0]

    labels = ws_split(mask, ws_min_dist, ws_thresh_rel, ws_sigma)

    #if watershed produced no regions
    if labels.max() == 0:
        labels = sklabel(mask.astype(bool), connectivity=2).astype(np.uint16)

    if out_labels:
        tiff.imwrite(out_labels, labels, dtype=np.uint16)
    if out_contours:
        tiff.imwrite(out_contours, labels_to_contours(labels), dtype=np.uint8)

    #unsplit contours from the pre-watershed mask
    if out_outer_contours:
        outer = labels_to_contours((mask > 0).astype(np.uint16))
        tiff.imwrite(out_outer_contours, outer, dtype=np.uint8)

    if labels.max() == 0:
        with open(out_txt, "w") as f:
            f.write("nan\n")
        open(out_csv, "w").close()
        return

    cents = [(p.centroid[1], p.centroid[0]) for p in regionprops(labels)]

    pairs = []
    dists_px = []
    thr = float(max_pair_px)
    for i, (x0, y0) in enumerate(cents):
        best_d, best_j = None, None
        for j, (x1, y1) in enumerate(cents):
            if i == j:
                continue
            d = float(np.hypot(x1 - x0, y1 - y0))
            if d < thr and (best_d is None or d < best_d):
                best_d, best_j = d, j
        if best_j is not None:
            pairs.append((x0, y0, cents[best_j][0], cents[best_j][1]))
            dists_px.append(best_d)

    with open(out_csv, "w") as f:
        for x0, y0, x1, y1 in pairs:
            f.write(f"{x0:.3f},{y0:.3f},{x1:.3f},{y1:.3f}\n")

    if dists_px:
        mean_units = float(np.mean(dists_px)) * float(px_size)
        with open(out_txt, "w") as f:
            f.write(f"{mean_units:.6f}\n")
    else:
        with open(out_txt, "w") as f:
            f.write("NaN\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=["thickness", "proc"], required=True)
    ap.add_argument("--mask", required=True)
    ap.add_argument("--px", type=float, required=True)
    ap.add_argument("--max_pair_px", type=float, default=20.0)
    ap.add_argument("--out_txt", required=True)
    ap.add_argument("--out_csv", required=True)
    # watershed knobs for PROCESS ONLY
    ap.add_argument("--ws_min_dist", type=float, default=3.28)
    ap.add_argument("--ws_thresh_rel", type=float, default=0.26)
    ap.add_argument("--ws_sigma", type=float, default=0.0)
    ap.add_argument("--out_labels", type=str, default=None)
    ap.add_argument("--out_contours", type=str, default=None)
    ap.add_argument("--out_outer_contours", type=str, default=None)
    args = ap.parse_args()

    if args.task == "thickness":
        gbm_thickness(args.mask, args.px, args.out_txt, args.out_csv)
    else:
        proc_ws_nnd(
            args.mask,
            args.px,
            args.max_pair_px,
            args.ws_min_dist,
            args.ws_thresh_rel,
            args.ws_sigma,
            args.out_txt,
            args.out_csv,
            args.out_labels,
            args.out_contours,
            args.out_outer_contours,
        )


if __name__ == "__main__":
    main()
