#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import scanpy as sc
import scConnect as cn
import matplotlib
import matplotlib.pyplot as plt
import sys
import random
random.seed(123)

counts_path = sys.argv[1]
meta_path = sys.argv[2]
species = sys.argv[3]
output_path = sys.argv[4]
adata = sc.read_csv(counts_path).T
meta = pd.read_csv(meta_path,  index_col="Cell")
adata.obs = meta

adata_tissue = cn.genecall.meanExpression(adata, groupby="Annotation", normalization=False, use_raw=False, transformation="log1p")
adata_tissue = cn.connect.ligands(adata_tissue, organism=species)
adata_tissue = cn.connect.receptors(adata_tissue, organism=species)
adata_tissue = cn.connect.specificity(adata_tissue, n=100, groupby="Annotation", organism=species)

edges = cn.connect.interactions(emitter=adata_tissue, target=adata_tissue, self_reference=True, organism=species)
nodes = cn.connect.nodes(adata_tissue)

edge = pd.DataFrame(edges)
col_new = list(edge.iloc[0,2].keys())
edge.rename(columns={0:'sender', 1:'reciever', 2:'info'}, inplace = True)
for each in col_new:
    edge[each] = edge["info"].map(lambda x:x[each])
    
edge.to_csv(output_path)

