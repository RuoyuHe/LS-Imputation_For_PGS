#!/bin/bash -l
#SBATCH -p ag2tb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=36:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256g
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=he000176@umn.edu

cd /home/panwei/he000176/jinwen/gln
conda activate eir-gln

BASE_PATH="/home/panwei/he000176/jinwen/gln/"

eirpredict --global_configs ${BASE_PATH}configs/global_config.yaml --input_configs ${BASE_PATH}configs/test/input_config.yaml --output_configs ${BASE_PATH}configs/test/output_config.yaml --model_path ${BASE_PATH}model_output/0.1/saved_models/0.1_model_450_perf-average=0.0380.pt --evaluate --output_folder /home/panwei/he000176/jinwen/results/gln
