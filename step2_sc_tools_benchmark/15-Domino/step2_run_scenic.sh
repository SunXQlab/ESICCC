pyscenic grn --num_workers 30 -o ./input_for_Domino/adj.tsv ./input_for_SCENIC/count.csv ./input_for_SCENIC/hs_hgnc_curated_tfs.txt --seed 123

pyscenic ctx ./input_for_Domino/adj.tsv ./input_for_SCENIC/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather ./input_for_SCENIC/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather --annotations_fname ./input_for_SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname ./input_for_SCENIC/count.csv --mode "dask_multiprocessing" --output ./input_for_Domino/reg.csv --num_workers 30

pyscenic aucell ./input_for_SCENIC/count.csv ./input_for_Domino/reg.csv -o ./input_for_Domino/auc_mtx.csv --num_workers 30 --seed 123
