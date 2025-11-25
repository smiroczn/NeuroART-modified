
from cellpose import models, io, plot
import numpy as np
import pandas as pd
from scipy.ndimage import center_of_mass
from scipy.io import loadmat
import tifffile
import argparse


def segment(calcium_image):
    CHANNELS = [0, 0]
    model = models.Cellpose(gpu=True, model_type='cyto3')

    # Use mean image across time for segmentation
    mean_img = calcium_image.mean(axis=0)

    # If your image has channels as the last axis (H x W x C), channel_axis=2 is correct.
    # If it's just H x W (single-channel), you can set channel_axis=None instead.
    masks, flows, styles, diams = model.eval(
        [mean_img],
        diameter=0,
        channels=CHANNELS,
        channel_axis=2
    )

    # masks is a list (one per input image); take the first
    return masks[0]


# Construct the argument parser
ap = argparse.ArgumentParser()

# Add the arguments to the parser
ap.add_argument("-p", "--fPath", required=True)         # file path
ap.add_argument("-n", "--nFrames", type=int, required=True)  # number of frames in the initial batch

args = vars(ap.parse_args())

n_frames = args['nFrames']
image_path = args['fPath']

# Read only the first n_frames from the multipage TIFF
with tifffile.TiffFile(image_path) as tif:
    calcium_image = tif.asarray(key=range(n_frames))

# --- Run Cellpose segmentation ---
masks = segment(calcium_image)

# === NEW BLOCK: save full label masks ===
# Save the 2D label image (each integer label = one cell, 0 = background)
tifffile.imwrite('masks.tif', masks.astype('uint16'))

# --- Compute centroids from masks (unchanged) ---
labels = np.unique(masks)        # Get unique labels for masks
labels = labels[labels > 0]      # Exclude background (label 0)

centroids = []
for label in labels:
    centroid = center_of_mass(masks == label)
    centroids.append(centroid)   # Each centroid is a tuple (y, x)

centroids_df = pd.DataFrame(centroids, columns=['Y', 'X'])
centroids_df.to_csv('centroids.csv', index=False)  # Save centroids to CSV
