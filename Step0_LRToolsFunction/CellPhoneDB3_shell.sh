source /home/ljx/miniconda3/bin/activate cpdb
cellphonedb method degs_analysis $1 $2 $3 --threshold 0.05 --output-path $4 --threads 10 --counts-data gene_name --debug-seed 123