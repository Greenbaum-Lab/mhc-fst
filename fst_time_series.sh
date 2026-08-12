#!/bin/bash
#SBATCH -p glacier
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

source /etc/profile.d/huji-lmod.sh
module load miniconda3/24.3.0-gcc-iqeknet
eval "$(conda shell.bash hook)"
conda activate adna

WORK_DIR=/sci/labs/gilig/lab_share/adna_db/ohad_mhc

python $WORK_DIR/run_fst_time_series.py \
	--config $WORK_DIR/config_fst_time.json

python $WORK_DIR/plot_fst_time_series.py \
	--results $WORK_DIR/results/fst_time_series.csv \
	--output-dir $WORK_DIR/results
