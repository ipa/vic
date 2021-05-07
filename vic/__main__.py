# -*- coding: utf-8 -*-
"""
@author: Iwan Paolucci
"""

import argparse
import os
import sys
from datetime import date

import numpy as np
import pandas as pd

from utils.niftireader import load_image
from vic import compute_vic, summarize_vic

np.set_printoptions(suppress=True, precision=4)

def get_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("-t", "--tumor", required=True, help="path to the tumor segmentation")
    ap.add_argument("-a", "--ablation", required=True, help="path to the ablation segmentation")
    ap.add_argument("-l", "--liver", required=False, help="path to the liver segmentation")
    ap.add_argument("-o", "--output", required=True, help="output margin (xlsx)")
    ap.add_argument("-s", "--subject-id")
    ap.add_argument("-ti", "--tumor-id")
    args = vars(ap.parse_args())
    return args


if __name__ == '__main__':

    args = get_args()
    tumor_file = args['tumor']
    ablation_file = args['ablation']
    liver_file = args['liver']
    output_file = args['output']
    subject_id = args['subject_id']
    tumor_id = args['tumor_id']

    tumor, tumor_np = load_image(tumor_file)
    has_tumor_segmented = np.sum(tumor_np.astype(np.uint8)) > 0
    if has_tumor_segmented is False:
        print('No tumor segmentation mask found in the file provided...program exiting')
        sys.exit(1)

    ablation, ablation_np = load_image(ablation_file)
    has_ablation_segmented = np.sum(ablation_np.astype(np.uint8)) > 0
    if has_ablation_segmented is False:
        print('No ablation segmentation mask found in the file provided...program exiting')
        sys.exit(1)

    if liver_file is not None:
        liver, liver_np = load_image(liver_file)
        has_liver_segmented = np.sum(liver_np.astype(np.uint8)) > 0
    else:
        has_liver_segmented = False

    pixdim = ablation.header['pixdim']
    spacing = (pixdim[1], pixdim[2], pixdim[3])

    vic_ml = compute_vic(tumor_np, ablation_np, liver_np, spacing)
    print(vic_ml)
    data = summarize_vic(subject_id, tumor_id, vic_ml)
    df = pd.DataFrame([data])
    df.to_csv(output_file, index=False, float_format='%.4f')
