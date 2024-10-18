# Trait imputation enhances nonlinear genetic prediction for some traits
This repository contains the code used for the analysis in [He, R., Fu, J., Ren, J., & Pan, W. (2024). Trait imputation enhances nonlinear genetic prediction for some traits. Genetics, iyae148.](https://doi.org/10.1093/genetics/iyae148)


## Description
The code is organized into three main parts:

1. **Preprocessing and sample splitting**: This part splits the data into GWAS data and analysis data, and further divides the analysis data into training, ensemble, and test sets (see the paper for details).
2. **Trait imputation**: We provide the code for the relatively new method [LS-imputation](https://doi.org/10.1016/j.xhgg.2023.100197).
3. **PGS model fitting**: We include the code for fitting Polygenic Score (PGS) models, specifically the tuning parameters for XGBoost and GLN.

For PRS-CS, we used the original authors' [code](https://github.com/getian107/PRScs).

## Data preprocessing
Detailed descriptions of the preprocessing steps can be found in [He, R., Fu, J., Ren, J., & Pan, W. (2024). Trait imputation enhances nonlinear genetic prediction for some traits. Genetics, iyae148.](https://doi.org/10.1093/genetics/iyae148)

The data preprocessing code is in the `./data_preprocessing/` folder:

- `get_phenotype.R`: Preprocesses phenotypes (including covariates) and performs sample splitting.
- `get_genotype.sh`: Contains PLINK commands for genotype preprocessing.


### GWAS
Both LS-imputation and PRS-CS use GWAS summary statistics for imputation. The code for running GWAS using PLINK is included.


### LS-imputation
LS-imputation is a new nonparametric method that uses GWAS summary statistics and genotype data to impute traits. Conceptually, it solves a "reverse regression" problem where, given genotype data (`X`) and GWAS coefficients (marginal projections), we ask the question: ask the question: what `Y` would give these projections? For more details and why LS-imputation might preserve some nonlinear genetic associations, see [Ren, J., Lin, Z., He, R., Shen, X., & Pan, W. (2023). Using GWAS summary data to impute traits for genotyped individuals. Human Genetics and Genomics Advances, 4(3).](https://doi.org/10.1016/j.xhgg.2023.100197)

The implementation can be found in `./LS-impute/LSimp.py`. Key functions include:

- `imputeY(snp, beta, batch_size=20000)`: Takes as input the entire genotype matrix (`snp`), GWAS coefficients (`beta`), and `batch_size`. The genotype matrix should be centered (columns with mean 0).
- `imputeY_batch(snp_batch, beta)`: Batch version of `imputeY`, performing one iteration over the batches.
- `imputeY_batch_cholesky(snp_batch, beta)`: Batch version using Cholesky decomposition for faster matrix inversion.

Additionally, the following files support LS-imputation:

- `read_data.R`: Contains R code to read the SNP matrix from `.bed` files.
- `impute_pheno.py`: Provides an example of using the `imputeY_batch_cholesky()` function from `LSimp.py`.


### XGBoost
The XGBoost implementation is located in `./xgboost/xgb.py`. The code includes all relevant tuning parameters required for training.


### GLN
A DL based model proposed in [Sigurdsson, A. I., Louloudis, I., Banasik, K., Westergaard, D., Winther, O., Lund, O., ... & Rasmussen, S. (2023). Deep integrative models for large-scale human genomics. Nucleic Acids Research, 51(12), e67-e67.](https://doi.org/10.1093/nar/gkad373)
The guide to using their code can be found [here](https://eir.readthedocs.io/en/latest/index.html).
In this repository, we provide configuration files for model training (in `./gln/configs`) and example scripts for training and predicting with the GLN model in the `./slurm-files/` directory.