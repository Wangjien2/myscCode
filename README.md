# Dara Type Convert

# Generate Seurat object

# Raw Filter

## Doublet Cells Predict

## Scrublet

## DoubletFinder

# Remove batch and run cluster

## BBKNN in R

To run the bbknn algorithm in the R environment, ensure you have installed bbknn in your Python environment. In this step, I recommend doing this in a conda virtual environment.

```bash
conda create -n bbknn_env python==3.8 -y -c conda-forge
conda activate bbknn_env
pip install bbknn
```

## BBKNN in python

To run BBKNN batch correction in Python, first of all, we should create a virtual environment. Here are the steps:

```powershell
# use conda create a virtual environment
conda create -n scanpy_env python==3.8  -c conda-forge -y
conda activate scanpy_env
# install python moudle
pip install bbknn==1.6.0
pip install scanpy==1.9.5
pip install scikit-learn==1.3.2
pip install leidenalg==0.10.1 matplotlib==3.6.3 numba==0.55.2 numpy==1.21.6 palantir==1.0.1 pandas==1.3.5


```


# Cell Type Annotation

# Cell Fraction In Different Group
