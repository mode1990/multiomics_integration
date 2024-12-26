# Multiomics Integration

## Overview
This repository contains a pipeline for integrating multi-modal single-cell data, specifically RNA-seq and ATAC-seq, using the MultiVI framework from scvi-tools.

## Pipeline Features
- Data loading and preprocessing for RNA-seq and ATAC-seq data
- Integration of multiple modalities using MultiVI
- Quality control and filtering
- Dimensionality reduction and visualization
- Comprehensive analysis of integrated data

## Prerequisites
```python
import scanpy as sc
import pandas as pd
import numpy as np
import scvi
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
