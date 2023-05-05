#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import cell2cell as c2c
import scanpy as sc
import pandas as pd
import sys
import random
random.seed(123)

counts_path = sys.argv[1]
meta_path = sys.argv[2]
lr_path = sys.argv[3]
output_path = sys.argv[4]

adata = sc.read_csv(counts_path).T
meta = pd.read_csv(meta_path, index_col = "Cell")
adata.obs = meta

lr_pairs = pd.read_csv(lr_path)
lr_pairs = lr_pairs.astype(str)
meta = adata.obs.copy()

interactions = c2c.analysis.SingleCellInteractions(rnaseq_data=adata.to_df().T,
                                                   ppi_data=lr_pairs,
                                                   metadata=meta,
                                                   interaction_columns=('ligand_symbol', 'receptor_symbol'),
                                                   communication_score='expression_thresholding',
                                                   expression_threshold=0.05, # values after aggregation
                                                   cci_score='bray_curtis',
                                                   cci_type='undirected',
                                                   aggregation_method='nn_cell_fraction',
                                                   barcode_col='Cell',
                                                   celltype_col='Annotation',
                                                   complex_sep='&',
                                                   verbose=False)
interactions.compute_pairwise_communication_scores()
# interactions.compute_pairwise_communication_scores() return result: 
ccc_matrix_path = output_path + '/communication_matrix.csv'
interactions.interaction_space.interaction_elements['communication_matrix'].to_csv(ccc_matrix_path)
ccc_pvals_path = output_path + '/ccc_pval.csv'
ccc_pvals = interactions.permute_cell_labels(evaluation='communication',verbose=True, random_state = 123)
ccc_pvals.to_csv(ccc_pvals_path)

