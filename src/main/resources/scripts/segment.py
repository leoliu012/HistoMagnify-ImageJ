
import cv2
import numpy as np
import tifffile
import os
import argparse
from PIL import Image, ImageSequence

from src.core.model_archi import multi_unet_model_trans
from src.core.segmentation import run_patches

MODEL_MAP = {'20x': '20x.hdf5', '40x': '40x.hdf5'}
MODEL_CHANNELS = {'20x': 1, '40x': 1}
MODEL_CLASSES  = {'20x': 6, '40x': 3}
MODEL_PARAMS = {
    '20x': dict(P_HEIGHT=720, P_WIDTH=960, MODEL_HEIGHT=576, MODEL_WIDTH=768),
    '40x': dict(P_HEIGHT=576, P_WIDTH=768, MODEL_HEIGHT=576, MODEL_WIDTH=768),
}

def read_page(tif_path, page_index):
    with Image.open(tif_path) as im:
        frame = None
        for i, f in enumerate(ImageSequence.Iterator(im)):
            if i == page_index:
                frame = f.copy()
                break
    if frame is None:
        raise RuntimeError(f"Page {page_index} not found")
    return np.array(frame)

def to_gray8(arr):
    a = np.asarray(arr)

    if a.dtype == np.bool_:
        a = a.astype(np.uint8) * 255

    if a.dtype == np.uint16:
        g = (a >> 8).astype(np.uint8)

    elif a.dtype == np.uint8:
        g = a

    else:
        a = a.astype(np.float32)
        amin, amax = float(np.nanmin(a)), float(np.nanmax(a))
        if amax > amin:
            g = ((a - amin) * (255.0 / (amax - amin))).astype(np.uint8)
        else:
            g = np.zeros_like(a, dtype=np.uint8)

    if g.ndim == 3:
        if g.shape[2] >= 3:
            g = cv2.cvtColor(g, cv2.COLOR_RGB2GRAY)
        else:
            g = g[..., 0]

    return g


def clahe_first_page(arr):
    g8 = to_gray8(arr)
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    return clahe.apply(g8)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--tif', required=True)
    p.add_argument('--page', type=int, default=0)
    p.add_argument('--modeldir', required=True)
    p.add_argument('--model', choices=MODEL_MAP.keys(), required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--enhance', type=int, default=1)   # 1 = CLAHE, 0 = skip
    p.add_argument('--plow',  type=float, default=1.0)
    p.add_argument('--phigh', type=float, default=99.7)
    p.add_argument('--save_enhanced', default=None)

    args = p.parse_args()

    plane = read_page(args.tif, args.page)

    if int(args.enhance) == 1:

        g8_for_model = clahe_first_page(plane)
        if args.save_enhanced:
            tifffile.imwrite(args.save_enhanced, g8_for_model)
    else:
        g8_for_model = to_gray8(plane)


    tmp_plane = os.path.join(args.modeldir, 'plane.tif')
    tifffile.imwrite(tmp_plane, g8_for_model)

    n_classes = MODEL_CLASSES[args.model]
    params    = MODEL_PARAMS[args.model]
    weights   = os.path.join(args.modeldir, MODEL_MAP[args.model])

    model = multi_unet_model_trans(
        n_classes=n_classes,
        IMG_HEIGHT=params['MODEL_HEIGHT'],
        IMG_WIDTH=params['MODEL_WIDTH'],
        IMG_CHANNELS=MODEL_CHANNELS[args.model]
    )
    model.load_weights(weights)

    seg_map = run_patches(
        tmp_plane,
        model,
        params['P_HEIGHT'],
        params['P_WIDTH'],
        n_classes,
        params['MODEL_WIDTH'],
        params['MODEL_HEIGHT']
    )

    tifffile.imwrite(args.output, seg_map.astype(np.uint8))
    print(f"Segmentation saved to {args.output}")

if __name__ == '__main__':
    main()
