#!/bin/bash -l
#SBATCH -p ag2tb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=36:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256g
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=he000176@umn.edu
#SBATCH --job-name=0.1_LS

cd /home/panwei/he000176/jinwen/gln
conda activate eir-gln

eirtrain --global_configs configs/global_config.yaml --input_configs configs/input_config.yaml --output_configs configs/output_config.yaml 
