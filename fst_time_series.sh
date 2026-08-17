#!/bin/bash
#SBATCH -p glacier
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

source /etc/profile.d/huji-lmod.sh
module load miniconda3/24.3.0-gcc-iqeknet
eval "$(conda shell.bash hook)"
conda activate adna

REPO_DIR=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
cd $REPO_DIR

python run_fst_time_series.py \
	--config config_fst_time.json

python plot_fst_time_series.py \
	--results results/fst_time_series.csv \
	--gene-background results/gene_background.csv \
	--output-dir results
