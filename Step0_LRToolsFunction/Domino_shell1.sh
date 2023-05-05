source /home/ljx/miniconda3/bin/activate pyscenic

pyscenic grn --num_workers 30 -o $1 $2 $5 --seed 123

pyscenic ctx $1 $6 $7 $8 $9 ${10} ${11} --annotations_fname ${12} --expression_mtx_fname $2 --mode "dask_multiprocessing" --output $3 --num_workers 30

pyscenic aucell $2 $3 -o $4 --num_workers 30 --seed 123