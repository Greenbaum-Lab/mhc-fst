#!/bin/bash
#SBATCH -p salmon
#SBATCH --gres=gpu:l40s:1
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-197

source /etc/profile.d/huji-lmod.sh
module load miniconda3/24.3.0-gcc-iqeknet
module load cuda/12.8.1
eval "$(conda shell.bash hook)"
conda activate adna

WORK_DIR=/sci/labs/gilig/lab_share/adna_db/ohad_mhc
FAST_FST=$WORK_DIR/fast_fst-master/fast_fst.py

BINS=(1000_2000 2000_3000 3000_4000 4000_5000 5000_6000 6000_7000 7000_8000 8000_9000 9000_10000)

bin_index=$((SLURM_ARRAY_TASK_ID / 22))
chromosome=$((SLURM_ARRAY_TASK_ID % 22 + 1))
bin_label=${BINS[$bin_index]}


python $FAST_FST \
	--bed-prefix $WORK_DIR/subsets/subset_${bin_label} \
	--output-prefix $WORK_DIR/fst/100k_window/${bin_label}/fst_${bin_label} \
	--chromosome $chromosome \
	--window-size 100000 \
	--threads 8
