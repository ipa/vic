import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="volume-of-insufficient-coverage",
    version="0.1.0",
    author="Iwan Paolucci",
    author_email="iwan.paolucci@gmail.com",
    description="Calculation of the volume of insufficient margin after ablation treatment",
    long_description="Implemented according to Kaye et al. Volumetric 3D assessment of ablation zones after thermal ablation of colorectal liver metastases to improve prediction of local tumor progression. ",
    long_description_content_type="text/markdown",
    url="https://github.com/ipa/vic",
    packages=setuptools.find_packages(),
    install_requires = [
        "nibabel>=3.2",
        "numpy>=1.19",
        "pandas>=1.1",
        "scipy>=1.5"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT ",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
