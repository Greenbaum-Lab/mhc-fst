#!/bin/bash
#SBATCH -p glacier
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

source /etc/profile.d/huji-lmod.sh
module load miniconda3/24.3.0-gcc-iqeknet
module load cuda/12.8.1
eval "$(conda shell.bash hook)"
conda activate adna

DATA_DIR=/sci/labs/gilig/lab_share/adna_db/poseidon_aadr
WORK_DIR=/sci/labs/gilig/lab_share/adna_db/ohad_mhc
REPO_DIR=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
cd $REPO_DIR

python -m window_scan.build_subsets \
	--janno $DATA_DIR/Poseidon_AADR_v62.janno \
	--source-bed-prefix $DATA_DIR/Poseidon_AADR_v62 \
	--config window_scan/config_full.json \
	--subset-prefix-template $WORK_DIR/subsets/subset_{low}_{high}
