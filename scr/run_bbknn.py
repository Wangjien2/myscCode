'''
Author: wangje
Date: 2025-08-11 10:12:52
'''
import scanpy as sc
import bbknn
import pandas as pd
import numpy as np
import os
import argparse
import sys

## 读入数据
adata = sc.read_h5ad(args.input)

