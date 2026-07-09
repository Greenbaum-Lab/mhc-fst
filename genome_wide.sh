#!/bin/bash
#SBATCH -p glacier
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-8

source /etc/profile.d/huji-lmod.sh
module load miniconda3/24.3.0-gcc-iqeknet
eval "$(conda shell.bash hook)"
conda activate adna

WORK_DIR=/sci/labs/gilig/lab_share/adna_db/ohad_mhc
COMB=$WORK_DIR/fast_fst-master/comb.py

BINS=(1000_2000 2000_3000 3000_4000 4000_5000 5000_6000 6000_7000 7000_8000 8000_9000 9000_10000)
bin_label=${BINS[$SLURM_ARRAY_TASK_ID]}


python $COMB \
	--prefix $WORK_DIR/fst/10k_window/${bin_label}/fst_${bin_label} \
	--window-size 1000000 \
	--chromosomes $(seq 1 22) \
	--output $WORK_DIR/fst/10k_window/gw/gw_${bin_label}.npy
