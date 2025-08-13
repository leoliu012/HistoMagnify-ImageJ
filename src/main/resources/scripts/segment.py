#!/usr/bin/env python3
#/src/main/resources/scripts
import os, argparse
import numpy as np
import copy
import tifffile
from nd2reader import ND2Reader
from src.core.model_archi import multi_unet_model_trans
from src.core.segmentation import run_patches
from src.core.image_enhance import (
    auto_enhance_single_channel,
    auto_enhance_multi_channel
)

MODEL_MAP      = {
    'ACTN4':      'ACTN4.hdf5',
    'DAPI':       'DAPI.hdf5',
    'NHS_SINGLE_CHANNEL': 'NHS_ester_single.hdf5',
    'NHS_COMBINED_ACTN4':    'NHS_ester_com.hdf5',
}
MODEL_CHANNELS = {'ACTN4':1,'DAPI':1,'NHS_SINGLE_CHANNEL':1,'NHS_COMBINED_ACTN4':2}
MODEL_CLASSES  = {'ACTN4':2,'DAPI':2,'NHS_SINGLE_CHANNEL':3,'NHS_COMBINED_ACTN4':3}
PATCH_SIZE = 576

def read_plane(nd2_path, z, c, t=0, c2=None):
    with ND2Reader(nd2_path) as rdr:
        rdr.bundle_axes = ['y','x']
        sizes = rdr.sizes
        iter_axes = [ax for ax in ('t','z','c') if sizes.get(ax,1) > 1]
        rdr.iter_axes = iter_axes

        frames = [frame for frame in rdr]
        stack = np.stack(frames, axis=0)  # flattened over iter_axes

        # multipliers for flat indexing
        lengths = [sizes[ax] for ax in iter_axes]
        multipliers = []
        for i in range(len(iter_axes)):
            prod = 1
            for L in lengths[i+1:]:
                prod *= L
            multipliers.append(prod)

        def flat_index(ti, zi, ci):
            picks = []
            for ax in iter_axes:
                if ax == 't':
                    picks.append(ti)
                elif ax == 'z':
                    picks.append(zi)
                elif ax == 'c':
                    picks.append(ci)
            return sum(pi * mul for pi, mul in zip(picks, multipliers)) if picks else 0

        if c2 is None:
            fi = flat_index(t, z, c)
            plane = stack[fi].astype(np.float32)
            return plane
        else:
            fi1 = flat_index(t, z, c)   # e.g., NHS
            fi2 = flat_index(t, z, c2)  # e.g., ACTN4
            plane1 = stack[fi1].astype(np.float32)
            plane2 = stack[fi2].astype(np.float32)
            return np.stack((plane1, plane2), axis=0)  # shape (2, H, W)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--nd2', required=True)
    p.add_argument('--z', type=int, required=True)
    p.add_argument('--channel', type=int, required=True,
                   help="0-based; use -1 for NHS_COM two-channel")
    p.add_argument('--channel2', type=int, default=-1,
                   help="Second channel (only for NHS_COMBINED_ACTN4)")
    p.add_argument('--modeldir', required=True)
    p.add_argument('--model', choices=MODEL_MAP, required=True)
    p.add_argument('--output', required=True,
                   help="segmentation label TIFF")
    p.add_argument('--plow', type=float, default=1.0)
    p.add_argument('--phigh', type=float, default=99.7)
    args = p.parse_args()

    #read & normalize
    # read & normalize (support 2-ch for combined model)
    if args.model == 'NHS_COMBINED_ACTN4':
        if args.channel2 < 0:
            raise ValueError("NHS_COMBINED_ACTN4 requires --channel2 (ACTN4 channel).")
        plane = read_plane(args.nd2, args.z, args.channel, c2=args.channel2).astype(np.float32)
    else:
        plane = read_plane(args.nd2, args.z, args.channel).astype(np.float32)

    mx = plane.max()
    plane = plane / mx if mx > 0 else plane



    #save to TIFF for enhancement & patching
    tmp = os.path.join(args.modeldir, 'plane.tif')
    tifffile.imwrite(tmp, (plane * 255).astype(np.uint8))

    # enhance in-place at tmp
    if plane.ndim == 3 and plane.shape[0] == 2:
        auto_enhance_multi_channel(tmp, tmp, p_low=args.plow, p_high=args.phigh)
    else:
        auto_enhance_single_channel(tmp, tmp, p_low=args.plow, p_high=args.phigh)

    #keep a copy for debugging
    enhanced = tifffile.imread(tmp)
    enhanced_path = os.path.join(args.modeldir, 'enhanced.tif')
    tifffile.imwrite(enhanced_path, enhanced)



    model = multi_unet_model_trans(
        n_classes=MODEL_CLASSES[args.model],
        IMG_HEIGHT=PATCH_SIZE,
        IMG_WIDTH=PATCH_SIZE,
        IMG_CHANNELS=MODEL_CHANNELS[args.model]
    )
    weights = os.path.join(args.modeldir, MODEL_MAP[args.model])
    model.load_weights(weights)

    seg_map = run_patches(
        tmp, model,
        PATCH_SIZE, PATCH_SIZE,
        MODEL_CLASSES[args.model],
        PATCH_SIZE, PATCH_SIZE
    )

    tifffile.imwrite(args.output, seg_map.astype(np.uint8))
    print(f"Segmentation saved to {args.output}")

if __name__=='__main__':
    main()
