#!/bin/bash -l

# the following are command-line arguments specifying the data path (path to phenotypes); genotype path; output path; pheno name isn't necessary
DATA_PATH=$1
BED_PATH=$2
OUT_PATH=$3
PHENO_NAME=$4

#DATA_PATH="/home/panwei/he000176/jinwen/data/pheno/sample_split1"
#BED_PATH="/home/panwei/he000176/jinwen/data/bed/lipA/pruned"
#OUT_PATH="/home/panwei/he000176/jinwen/results/GWAS/GWAS_data/lipA"

echo "Starting GWAS SLURM job chr${SLURM_ARRAY_TASK_ID} with phenotype: $PHENO_NAME" | tee -a $LOG_FILE

plink2 --bfile ${BED_PATH}/chr${SLURM_ARRAY_TASK_ID} --keep ${DATA_PATH}/gwas_ids_forGWAS.txt --pheno ${DATA_PATH}/truepheno_gwas_adjusted.txt --pheno-name pheno --glm hide-covar --out ${OUT_PATH}/chr${SLURM_ARRAY_TASK_ID}

echo "Completed GWAS SLURM job chr${SLURM_ARRAY_TASK_ID} with phenotype: $PHENO_NAME" | tee -a $LOG_FILE

#module purge
