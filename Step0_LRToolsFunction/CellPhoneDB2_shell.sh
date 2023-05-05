source /home/ljx/miniconda3/bin/activate cpdb2
cellphonedb method statistical_analysis $1 $2 --threshold 0.05 --output-path $3 --threads 10 --counts-data gene_name --debug-seed 123
