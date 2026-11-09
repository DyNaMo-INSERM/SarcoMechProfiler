# SarcoMechProfiler

## Introduction

SarcoMechProfiler is a Python-based analysis toolkit for finding periodic patterns, normalizing the lengths, and averaging across line profiles of viscoelastic (here, sarcomeres) resulting from AFM mapping.  
Briefly, the following analyses are implemented in the jupyter notebooks:

The example input dataset for this averaging is resulting line profiles drawn over Scale factor for stiffness (E0) and fluidity index (beta) images, using an [ImageJ/Fiji macro](AFM_extract_line_profiles.ijm), from the AFM-based viscoelastic analysis from the [PyFMLab](https://github.com/DyNaMo-INSERM/PyFMlab_DyNaMo). This example excel file consists of E0 and beta of a single sarcomere extracted from a mechananical map from an AFM of the indirect flight muscle of a fly of sls\[23-YPet\] genotype, used in the manuscript [`example_data`](data/example_data.csv). In this `.csv` file, there is the distance, E0 and beta and file_sac number and date to identify uniquely the experimental file id.

- *Normalization and Averaging*: This involves aligning individual sarcomere patterns in E0, normalizing their lengths, and averaging the viscoelastic properties (E0 & beta) over all probed sarcomeres. An example implementation for the test dataset is provided in the [norm_averaging_VE.ipynb](norm_averaging_VE.ipynb) Jupyter notebook.

## Step-by-Step Guide

- Clone the repository

```
git clone 
cd SarcoMechProfiler
```

- Create an environment with python 3.10

```
conda create -n yourenvname python=3.10 
conda activate yourenvname
```

- Install the dependencies from requirements.txt

```
pip install -r requirements.txt
```

- Launch Jupyter Notebook

```
jupyter lab
```

- Open and run the notebooks. The repository contains the example dataset and the following analysis notebooks: [norm_averaging_VE.ipynb](norm_averaging_VE.ipynb) – Normalizing and averaging analysis

## Citation

If you use **SarcoMechProfiler** in your research, please cite the associated preprint:

Purohit D, Saravanan Y, Loreau V, Rico F, Schnorrer F **Titin PEVK length determines the viscoelasticity of Drosophila sarcomeres**. (2026).

## Contact

If you have any suggestions, comments or encounter any issues, please write to anyone of us: yogesh.saravanan@inserm.fr, devam.PUROHIT@univ-amu.fr, felix.rico@inserm.fr or frank.SCHNORRER@univ-amu.fr