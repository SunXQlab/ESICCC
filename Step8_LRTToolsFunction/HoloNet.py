#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import HoloNet as hn

import os
import pandas as pd
import random
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import torch
import sys

import warnings
warnings.filterwarnings('ignore')
hn.set_figure_params(tex_fonts=False)
sc.settings.figdir = './figures/'


# In[ ]:


random.seed(123) 
np.random.seed(123)

counts_path = sys.argv[1]
meta_path = sys.argv[2]
img_path = sys.argv[3]
goi_fpath = sys.argv[4]
output_path = sys.argv[5]


# In[ ]:


# load ST data (gene expression matrix + celltype + spatial)
adata = sc.read_csv(counts_path).T
meta = pd.read_csv(meta_path,  index_col="Barcodes")
adata.obs = meta
adata = sc.read_visium(adata, path = img_path)


# In[ ]:


# run normalization
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)


# In[ ]:


# load LR database and filter genes epressed by less than 5% of cells
LR_df = hn.pp.load_lr_df()
expressed_LR_df = hn.pp.get_expressed_lr_df(LR_df, adata, expressed_proportion=0.05)


# In[ ]:


# calculate w_best
w_best = hn.tl.default_w_visium(adata)
w_best


# In[ ]:


# construct CE
CE_tensor = hn.tl.compute_ce_tensor(adata, lr_df=expressed_LR_df, w_best=w_best)
CE_tensor_filtered = hn.tl.filter_ce_tensor(CE_tensor, adata,
                                            lr_df=expressed_LR_df, w_best=w_best)


# In[ ]:


# Selecting the target gene (icgs) to be predicted
target_all_gene_expr, used_gene_list = hn.pr.get_gene_expr(adata, expressed_LR_df,
                                                           max_sparse = 0.05)


# In[ ]:


goi = pd.read_csv(goi_fpath)
goi = list(set(goi["x"]).intersection(set(used_gene_list)))
len(goi)


# In[ ]:


X, cell_type_names = hn.pr.get_one_hot_cell_type_tensor(adata, categorical_cell_type_col = 'celltype')
adj = hn.pr.adj_normalize(adj=CE_tensor_filtered, cell_type_tensor=X, only_between_cell_type=True)

for gene in goi:
  target = hn.pr.get_one_case_expr(target_all_gene_expr, cases_list=used_gene_list,
                                   used_case_name= gene)
  trained_MGC_model_list = hn.pr.mgc_repeat_training(X, adj, target, device='gpu')
  predict_result = hn.pl.plot_mgc_result(trained_MGC_model_list, adata, X, adj)
  for lr in list(expressed_LR_df['LR_Pair']):
    _ = hn.pl.fce_cell_type_network_plot(trained_MGC_model_list, expressed_LR_df, X, adj,
                                       cell_type_names, plot_lr=lr, edge_thres=0.2,
                                       )
    res = pd.DataFrame(_[0])
    file_path = output_path + lr + "_" + gene + ".csv"
    res.to_csv(file_path)

