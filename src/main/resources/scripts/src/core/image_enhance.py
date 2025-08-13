import tifffile
from PIL import Image, ImageSequence
import numpy as np
import cv2
import os
import glob


def auto_enhance_single_channel(input_path, out_path, p_low=1, p_high=99.7):
    pil = Image.open(input_path)
    page = next(ImageSequence.Iterator(pil))
    first = page.copy()
    pil.close()

    arr = np.array(first).astype(np.float32)

    lo = np.percentile(arr, p_low)
    hi = np.percentile(arr, p_high)

    # arr * contrast + brightness
    if hi > lo:
        contrast = 255.0 / (hi - lo)
        brightness = -lo * contrast
        stretched = arr * contrast + brightness
    else:
        stretched = arr

    stretched = np.clip(stretched, 0, 255).astype(np.uint8)

    # enhanced = clahe.apply(stretched)

    out_img = Image.fromarray(stretched)
    out_img.save(out_path, format="TIFF")
    return out_img


def auto_enhance_multi_channel(input_path, out_path, p_low=1, p_high=99.7):
    pil = Image.open(input_path)
    pages = [page.copy() for page in ImageSequence.Iterator(pil)]
    pil.close()

    enhanced = []
    for page in pages:
        arr = np.array(page).astype(np.float32)
        lo = np.percentile(arr, p_low)
        hi = np.percentile(arr, p_high)

        if hi > lo:
            contrast = 255.0 / (hi - lo)
            brightness = -lo * contrast
            stretched = arr * contrast + brightness
        else:
            stretched = arr  # flat if no dynamic range
        stretched = np.clip(stretched, 0, 255).astype(np.uint8)
        enhanced.append(stretched)

    # enhanced = clahe.apply(stretched)
    combined = np.stack(enhanced, axis=0)



    tifffile.imwrite(out_path, combined, photometric='minisblack')
    print(combined.shape)
    return combined