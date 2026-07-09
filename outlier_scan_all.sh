#!/bin/bash
#SBATCH -p glacier
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

source /etc/profile.d/huji-lmod.sh
module load miniconda3/24.3.0-gcc-iqeknet
eval "$(conda shell.bash hook)"
conda activate adna

WORK_DIR=/sci/labs/gilig/lab_share/adna_db/ohad_mhc

BINS=(1000_2000 2000_3000 3000_4000 4000_5000 5000_6000 6000_7000 7000_8000 8000_9000 9000_10000)


for bin_label in "${BINS[@]}"; do
	python $WORK_DIR/outlier_scan.py \
		--fst-dir $WORK_DIR/fst/100k_window \
		--bin $bin_label \
		--chromosomes $(seq 1 22) \
		--window-size 100000 \
		--pair-index 0 \
		--stat-index 2 \
		--min-snps 20 \
		--q-threshold 0.05 \
		--output-path $WORK_DIR/fst/100k_window/signific/outliers_${bin_label}.csv
done
