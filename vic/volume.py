import numpy as np
from scipy import ndimage
import pandas as pd
import nibabel as nib

def compute_vic(tumor_mask, ablation_mask, liver_mask, spacing_mm, exclusion_distance=5):
    """
    Function computing the surface distances between 2 binary segmentation Nifti images.
    :param tumor_mask: tumor file in Nifti (NiBabel) format
    :param ablation_mask: ablation file in Nifti (NiBabel) format
    :param liver_mask: liver file Nifti (NiBabel) format
    :param spacing_mm: spacing extracted from the Nifti files provided above. spacing should be the same for all.
    :param exclusion_distance: The exclusion distance to "remove" voxels from the liver capsule within this distance. Works only for subcapsular cases when liver segmentation provided.
    :return: List of Volumes of insufficient ablation with 0-10 mm margin
    """

    distmap_mask =  ndimage.morphology.distance_transform_edt(~tumor_mask, sampling=spacing_mm)

    # distmap_nib = nib.Nifti1Image(distmap_mask.astype(np.uint8), affine = np.eye(4))
    # nib.save(distmap_nib, "distmap.nii.gz")

    voxel_volume = spacing_mm[0] * spacing_mm[1] * spacing_mm[2]

    ablation_mask = ablation_mask.astype(np.int8)

    vic_ml = []

    for i in range(0, 11):
        volume_added_margin = (distmap_mask <= i).astype(np.int8)
        # volume_added_margin_nib = nib.Nifti1Image(volume_added_margin.astype(np.int8), affine = np.eye(4))
        # nib.save(volume_added_margin_nib, "added_margin_{0}.nii.gz".format(i))
        ablation_volume = np.count_nonzero(ablation_mask) * voxel_volume / 1000.0
        tumor_volume = np.count_nonzero(volume_added_margin) * voxel_volume / 1000.0
        insufficient_volume = (volume_added_margin - ablation_mask) > 0
        # insufficient_nib = nib.Nifti1Image(insufficient_volume.astype(np.int8), affine = np.eye(4))
        # nib.save(insufficient_nib, "insufficient_{0}.nii.gz".format(i))
        volume_insufficient_margin = np.count_nonzero(insufficient_volume) * voxel_volume
        vic_ml.append(volume_insufficient_margin / 1000.0)
        print("VIC {0} = {1}\t({2} - {3})".format(i, vic_ml[i], tumor_volume, ablation_volume))

    return vic_ml


def summarize_vic(subject_id, tumor_id, vic_ml):
    df = pd.DataFrame(data=[{
            'Subject': subject_id,
            'Tumor': tumor_id,
            'VIC_0': vic_ml[0],
            'VIC_1': vic_ml[1],
            'VIC_2': vic_ml[2],
            'VIC_3': vic_ml[3],
            'VIC_4': vic_ml[4],
            'VIC_5': vic_ml[5],
            'VIC_6': vic_ml[6],
            'VIC_7': vic_ml[7],
            'VIC_8': vic_ml[8],
            'VIC_9': vic_ml[9],
            'VIC_10': vic_ml[10]}])
    return df
