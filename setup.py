import setuptools
from pathlib import Path

source_root = Path(".")

with open(source_root / "README.md", "r") as fh:
    long_description = fh.read()

version = "0.1.2"

with open(source_root / "vic" / "version.py", "w") as fh:
    fh.writelines([
        f'__version__ = "{version}" \n '
    ])


setuptools.setup(
    name="volume-of-insufficient-coverage",
    version=version,
    author="Iwan Paolucci",
    author_email="iwan.paolucci@gmail.com",
    description="Calculation of the volume of insufficient margin after ablation treatment",
    long_description=long_description,
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
