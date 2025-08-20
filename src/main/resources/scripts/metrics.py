#!/usr/bin/env python3
import argparse
import numpy as np
import tifffile as tiff
from scipy import ndimage as ndi
from skimage.morphology import medial_axis
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from skimage.measure import regionprops, label as sklabel
from scipy.ndimage import binary_fill_holes as _fill_holes



def _read_plane_bool(path):
    a = tiff.imread(path)
    a = np.asarray(a).squeeze()
    if a.ndim > 2:
        a = a[..., 0]
    return (a > 0)


def gbm_thickness(mask_path, out_txt, out_csv):
    m = _read_plane_bool(mask_path)
    if not np.any(m):
        with open(out_txt, "w") as f: f.write("nan")
        open(out_csv, "w").close()
        return

    skel, dist = medial_axis(m, return_distance=True)
    diam_px = 2.0 * dist[skel]
    mean_px = float(np.mean(diam_px)) if diam_px.size else float("nan")

    ys, xs = np.nonzero(skel)
    with open(out_csv, "w") as f:
        for x, y, d in zip(xs, ys, diam_px):
            f.write(f"{int(x)},{int(y)},{d:.6f}\n")  #in px
    with open(out_txt, "w") as f:
        f.write(f"{mean_px:.6f}")



def ws_split(mask_u8, min_dist, thresh_rel, sigma):
    mask = mask_u8.astype(bool)
    if not np.any(mask):
        return np.zeros_like(mask, dtype=np.uint16)

    dist = ndi.distance_transform_edt(mask)

    sigma = float(sigma)
    dist_s = ndi.gaussian_filter(dist, sigma=sigma) if sigma > 0.0 else dist

    #prevents 1 px peak
    md = max(1, int(round(min_dist)))
    fp = np.ones((2 * md + 1, 2 * md + 1), dtype=bool)

    peaks = peak_local_max(
        dist_s,
        labels=mask,
        footprint=fp,
        min_distance=md,
        threshold_rel=float(thresh_rel),
        exclude_border=False,
    )
    markers = np.zeros_like(dist, dtype=np.int32)
    for i, (r, c) in enumerate(peaks, start=1):
        markers[r, c] = i

    if markers.max() == 0:
        return sklabel(mask, connectivity=2).astype(np.uint16)

    grad = ndi.gaussian_gradient_magnitude(dist_s.astype(np.float32), sigma=sigma)
    labels = watershed(grad, markers=markers, mask=mask)
    return labels.astype(np.uint16)



def labels_to_contours(labels_u16):
    lab = labels_u16.astype(np.int32)
    h, w = lab.shape
    edges = np.zeros((h, w), np.uint8)

    edges[1:, :]  |= ((lab[1:, :]  != 0) & (lab[1:, :]  != lab[:-1, :])).astype(np.uint8)
    edges[:-1, :] |= ((lab[:-1, :] != 0) & (lab[:-1, :] != lab[1:, :])).astype(np.uint8)
    edges[:, 1:]  |= ((lab[:, 1:]  != 0) & (lab[:, 1:]  != lab[:, :-1])).astype(np.uint8)
    edges[:, :-1] |= ((lab[:, :-1] != 0) & (lab[:, :-1] != lab[:, 1:])).astype(np.uint8)

    return (edges * 255).astype(np.uint8)

def _merge_holes(mask_bool):
    return _fill_holes(mask_bool).astype(np.uint8)

def count_components(mask_path, out_txt):
    m = _read_plane_bool(mask_path)
    structure = np.array([[1,1,1],
                          [1,1,1],
                          [1,1,1]], dtype=np.uint8)
    _, num = ndi.label(m.astype(np.uint8), structure=structure)
    with open(out_txt, "w") as f:
        f.write(str(int(num)) + "\n")


def _percentile_bounds(vals, keep_low, keep_high):
    if len(vals) == 0:
        return -np.inf, np.inf
    keep_low  = float(max(0.0, min(1.0, keep_low)))
    keep_high = float(max(0.0, min(1.0, keep_high)))
    lo = -np.inf if keep_low  == 0.0 else float(np.percentile(vals, keep_low * 100.0))
    hi =  np.inf if keep_high == 0.0 else float(np.percentile(vals, (1.0 - keep_high) * 100.0))
    return lo, hi



def nuc_rbc_count(
    mask_path,
    keep_low,
    keep_high,
    ws_min_dist,
    ws_thresh_rel,
    ws_sigma,
    out_txt,
    out_labels=None,
    out_contours=None,
    out_outer_contours=None,
):
    raw = _read_plane_bool(mask_path)
    mask = _merge_holes(raw)

    if not np.any(mask):
        if out_labels: tiff.imwrite(out_labels, np.zeros_like(mask, np.uint16))
        if out_contours: tiff.imwrite(out_contours, np.zeros_like(mask, np.uint8))
        if out_outer_contours: tiff.imwrite(out_outer_contours, np.zeros_like(mask, np.uint8))
        with open(out_txt, "w") as f: f.write("0\n")
        return

    #parents
    cc_map = sklabel(mask, connectivity=2)
    parent_counts = np.bincount(cc_map.ravel())[1:]
    parent_ids = list(range(1, len(parent_counts) + 1))
    areas_par = parent_counts.tolist()

    #children
    labels = ws_split(mask.astype(np.uint8), ws_min_dist, ws_thresh_rel, ws_sigma).astype(np.uint16)
    if labels.max() == 0:
        labels = sklabel(mask, connectivity=2).astype(np.uint16)
    child_counts = np.bincount(labels.ravel())[1:] if labels.max() > 0 else np.array([], dtype=int)
    child_ids = list(range(1, len(child_counts) + 1))
    areas_child = child_counts.tolist()

    all_areas = areas_par + areas_child
    lo, hi = _percentile_bounds(all_areas, keep_low, keep_high)

    kept_parents  = {pid for (pid, a) in zip(parent_ids, areas_par)  if (a > lo) and (a < hi)}
    kept_children = {cid for (cid, a) in zip(child_ids,  areas_child) if (a > lo) and (a < hi)}

    frags_by_parent = {}
    for pid in parent_ids:
        comp = (cc_map == pid)
        kids = set(np.unique(labels[comp])) - {0}
        frags_by_parent[pid] = kids

    parents_with_kept_kids = {p for p, kids in frags_by_parent.items() if len(kept_children.intersection(kids)) > 0}
    parents_for_metrics = kept_parents - parents_with_kept_kids

    kept_uns_mask     = np.isin(cc_map, list(kept_parents))
    kept_split_labels = np.where(np.isin(labels, list(kept_children)), labels, 0).astype(np.uint16)

    if out_labels:
        tiff.imwrite(out_labels, kept_split_labels)
    if out_contours:
        tiff.imwrite(out_contours, labels_to_contours(kept_split_labels))
    if out_outer_contours:
        tiff.imwrite(out_outer_contours, labels_to_contours(kept_uns_mask.astype(np.uint16)))

    total = len(kept_children) + len(parents_for_metrics)
    with open(out_txt, "w") as f:
        f.write(str(int(total)) + "\n")



def proc_ws_nnd(
    mask_path,
    max_pair_px,
    ws_min_dist,
    ws_thresh_rel,
    ws_sigma,
    out_txt,
    out_csv,
    out_labels=None,
    out_contours=None,
    out_outer_contours=None,
    keep_low=0.0,
    keep_high=0.0,
):
    raw = _read_plane_bool(mask_path)
    mask = _merge_holes(raw)
    if not np.any(mask):
        if out_labels: tiff.imwrite(out_labels, np.zeros_like(mask, np.uint16))
        if out_contours: tiff.imwrite(out_contours, np.zeros_like(mask, np.uint8))
        if out_outer_contours: tiff.imwrite(out_outer_contours, np.zeros_like(mask, np.uint8))
        if out_csv: open(out_csv, "w").close()
        with open(out_txt, "w") as f: f.write("nan\n")
        return

    #all pabels
    cc_map = sklabel(mask, connectivity=2)
    parent_counts = np.bincount(cc_map.ravel())[1:]
    parent_ids = list(range(1, len(parent_counts) + 1))
    areas_par = parent_counts.tolist()

    #child labels
    labels = ws_split(mask.astype(np.uint8), ws_min_dist, ws_thresh_rel, ws_sigma).astype(np.uint16)
    if labels.max() == 0:
        labels = sklabel(mask, connectivity=2).astype(np.uint16)
    child_counts = np.bincount(labels.ravel())[1:] if labels.max() > 0 else np.array([], dtype=int)
    child_ids = list(range(1, len(child_counts) + 1))
    areas_child = child_counts.tolist()


    all_areas = areas_par + areas_child
    lo, hi = _percentile_bounds(all_areas, keep_low, keep_high)
    kept_parents = {pid for (pid, a) in zip(parent_ids, areas_par)  if (a > lo) and (a < hi)}
    kept_children = {cid for (cid, a) in zip(child_ids,  areas_child) if (a > lo) and (a < hi)}

    # exclude parent with kept children
    frags_by_parent = {}
    for pid in parent_ids:
        comp = (cc_map == pid)
        kids = set(np.unique(labels[comp])) - {0}
        frags_by_parent[pid] = kids
    parents_with_kept_kids = {p for p, kids in frags_by_parent.items() if len(kept_children.intersection(kids)) > 0}
    parents_for_metrics = kept_parents - parents_with_kept_kids

    kept_uns_mask = np.isin(cc_map, list(kept_parents))
    kept_split_labels = np.where(np.isin(labels, list(kept_children)), labels, 0).astype(np.uint16)

    if out_labels:
        tiff.imwrite(out_labels, kept_split_labels)
    if out_contours:
        tiff.imwrite(out_contours, labels_to_contours(kept_split_labels)) # split contours
    if out_outer_contours:
        tiff.imwrite(out_outer_contours, labels_to_contours(kept_uns_mask.astype(np.uint16)))  # parent contours

    #centroids-children and parents_for_metrics
    cents = []
    for pid in parents_for_metrics:
        yy, xx = np.nonzero(cc_map == pid)
        if len(xx) > 0:
            cents.append((float(xx.mean()), float(yy.mean())))
    for rp in regionprops(kept_split_labels):
        y, x = rp.centroid
        cents.append((float(x), float(y)))

    pairs, dists_px = [], []
    thr_max = float(max_pair_px)
    thr_min = float(ws_min_dist)

    if thr_min < thr_max:
        for i, (x0, y0) in enumerate(cents):
            best_d, best_j = None, None
            for j, (x1, y1) in enumerate(cents):
                if i == j:
                    continue
                d = float(np.hypot(x1 - x0, y1 - y0))
                if (d >= thr_min) and (d < thr_max) and (best_d is None or d < best_d):
                    best_d, best_j = d, j
            if best_j is not None:
                pairs.append((x0, y0, cents[best_j][0], cents[best_j][1]))
                dists_px.append(best_d)

    # CSV + mean
    if out_csv:
        with open(out_csv, "w") as f:
            for x0, y0, x1, y1 in pairs:
                f.write(f"{x0:.3f},{y0:.3f},{x1:.3f},{y1:.3f}\n")

    if dists_px:
        mean_px = float(np.mean(dists_px)) if dists_px else float("nan")
        with open(out_txt, "w") as f:
            f.write((f"{mean_px:.6f}\n") if not np.isnan(mean_px) else "nan\n")
    else:
        with open(out_txt, "w") as f:
            f.write("nan\n")




def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=["thickness", "proc", "nuc", "rbc"], required=True)
    ap.add_argument("--mask", required=True)
    ap.add_argument("--max_pair_px", type=float, default=20.0)
    ap.add_argument("--out_txt", required=True)
    # CSV is required for thickness/proc, optional/ignored for nuc/rbc
    ap.add_argument("--out_csv", required=False, default=None)
    # watershed knobs for PROCESS ONLY
    ap.add_argument("--ws_min_dist", type=float, default=3.28)
    ap.add_argument("--ws_thresh_rel", type=float, default=0.26)
    ap.add_argument("--ws_sigma", type=float, default=0.0)
    ap.add_argument("--out_labels", type=str, default=None)
    ap.add_argument("--out_contours", type=str, default=None)
    ap.add_argument("--out_outer_contours", type=str, default=None)
    ap.add_argument("--keep_low", type=float, default=0.0)
    ap.add_argument("--keep_high", type=float, default=0.0)
    args = ap.parse_args()

    if args.task in ("nuc", "rbc"):
        nuc_rbc_count(
           mask_path=args.mask,
           keep_low=args.keep_low,
           keep_high=args.keep_high,
           ws_min_dist=args.ws_min_dist,
           ws_thresh_rel=args.ws_thresh_rel,
           ws_sigma=args.ws_sigma,
           out_txt=args.out_txt,
           out_labels=getattr(args, "out_labels", None),
           out_contours=getattr(args, "out_contours", None),
           out_outer_contours=getattr(args, "out_outer_contours", None),
       )
        return


    # tasks below require CSV
    if args.out_csv is None:
        raise SystemExit("For task '%s', --out_csv is required." % args.task)

    if args.task == "thickness":
        gbm_thickness(args.mask, args.out_txt, args.out_csv)
    else:
        proc_ws_nnd(args.mask, args.max_pair_px, args.ws_min_dist, args.ws_thresh_rel, args.ws_sigma,
                    args.out_txt, args.out_csv, args.out_labels, args.out_contours, args.out_outer_contours,
                    keep_low=args.keep_low, keep_high=args.keep_high)




if __name__ == "__main__":
    main()
