import numpy as np
from scipy import ndimage

def compute_vic(tumor_mask, ablation_mask, liver_mask, spacing_mm, exclusion_distance=5):
    """
    Function computing the surface distances between 2 binary segmentation Nifti images.
    :param tumor_mask: tumor file in numpy format
    :param ablation_mask: ablation file innumpy format
    :param liver_mask: liver file numpy format
    :param spacing_mm: spacing of the volumes provided above. spacing should be the same for all.
    :param exclusion_distance: The exclusion distance to "remove" voxels from the liver mask.
    :return: List of Volumes of insufficient ablation with 0-10 mm margin
    """

    distmap_mask =  ndimage.morphology.distance_transform_edt(~tumor_mask, sampling=spacing_mm)
    voxel_volume = spacing_mm[0] * spacing_mm[1] * spacing_mm[2]

    ablation_mask = ablation_mask.astype(np.int8)
    if liver_mask is not None:
        # print('has liver mask')
        liver_mask = liver_mask.astype(np.int8)
        distmap_exclusion = ndimage.morphology.distance_transform_edt(liver_mask, sampling=spacing_mm)
        distmap_mask[distmap_exclusion <= exclusion_distance] = np.inf

    vic_ml = []

    for i in range(0, 11):
        volume_added_margin = (distmap_mask <= i).astype(np.int8)
        ablation_volume = np.count_nonzero(ablation_mask) * voxel_volume / 1000.0
        tumor_volume = np.count_nonzero(volume_added_margin) * voxel_volume / 1000.0
        insufficient_volume = (volume_added_margin - ablation_mask) > 0
        volume_insufficient_margin = np.count_nonzero(insufficient_volume) * voxel_volume
        vic_ml.append(volume_insufficient_margin / 1000.0)
        # print("VIC {0} = {1}\t({2} - {3})".format(i, vic_ml[i], tumor_volume, ablation_volume))

    return vic_ml


def summarize_vic(subject_id, tumor_id, vic_ml):
    data={
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
            'VIC_10': vic_ml[10]}
    return data
