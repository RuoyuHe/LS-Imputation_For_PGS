#!/bin/bash -l

# the following are command line arguments

PHENOTYPE_PATH=$1
BED_PATH=$2
LOG_FILE=$3
PHENO_NAME=$4

echo "Starting pruning SLURM job ${SLURM_ARRAY_TASK_ID} with phenotype path: $PHENOTYPE_PATH, bed file path: $BED_PATH" | tee -a $LOG_FILE

module load plink/2.00-alpha-091019; plink2 --pfile /home/panwei/shared/UKBiobankIndiv/imputed/pgen/ukbb_chr${SLURM_ARRAY_TASK_ID}_1 --chr ${SLURM_ARRAY_TASK_ID} --maf 0.05 --geno 0.1 --hwe 0.001 --indep-pairwise 50 5 0.8 --keep ${PHENOTYPE_PATH}/all_ids.txt --make-bed --out ${BED_PATH}/chr${SLURM_ARRAY_TASK_ID}

plink2 --bfile ${BED_PATH}/chr${SLURM_ARRAY_TASK_ID} --chr ${SLURM_ARRAY_TASK_ID} --extract ${BED_PATH}/chr${SLURM_ARRAY_TASK_ID}.prune.in --make-bed --out ${BED_PATH}/pruned/chr${SLURM_ARRAY_TASK_ID}

echo "Completed pruning SLURM job ${SLURM_ARRAY_TASK_ID} creating bed files and SNP pruning" | tee -a $LOG_FILE
