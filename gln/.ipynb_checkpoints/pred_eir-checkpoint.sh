#!/bin/bash -l
#SBATCH -p ag2tb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=80g
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=he000176@umn.edu

cd 
module load conda
conda activate eir-gln

python3 -m eirtrain --global_configs configs/global_config.yaml --input_configs configs/input_config.yaml --output_configs configs/output_config.yaml
