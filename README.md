# Volume of Insufficient Coverage (VIC)

***!!!This package is under active development and not save for work!!!***

## Credits

This package is an implementation of the VIC metric describe by Kaye et al. to assess the coverage of thermal ablations.

Kaye, E. A., Cornelis, F. H., Petre, E. N., Tyagi, N., Shady, W., Shi, W., Zhang, Z., Solomon, S. B., Sofocleous, C. T. & Durack, J. C. Volumetric 3D assessment of ablation zones after thermal ablation of colorectal liver metastases to improve prediction of local tumor progression. Eur. Radiol. 29, 2698–2705 (2019). https://doi.org/10.1007/s00330-018-5809-0

## Installation

    python setup.py install

## Usage

    python -m vic --tumor [tumor.nii] --ablation [ablation.nii] --liver [liver.nii] --output [output.csv] --subject-id [patient-1] --tumor-id [tumor-1]
