import sys
import math
import scipy
import scipy.linalg 
import numpy as np
import json
import os
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
from pandas_plink import get_data_folder
from LSimp import imputeY, imputeY_batch_cholesky
# from importlib import reload

impute_data = sys.argv[1]
batch_size_idx = int(sys.argv[2])
batch_idx = int(sys.argv[3])

batch_sizes = [7500,15000,25000,35000]
batch_size = batch_sizes[batch_size_idx-1]

save_path = '/home/panwei/he000176/jinwen/data/pheno/'

r = robjects.r
r['source']('/home/panwei/he000176/jinwen/LS-impute/read_data.R')

# Loading the function we have defined in R.
read_data = robjects.globalenv['read_data']

with localconverter(robjects.default_converter + pandas2ri.converter):
    data = read_data(impute_data,batch_size_idx,batch_idx)
    sample_id = robjects.conversion.rpy2py(data[0])
    snp = robjects.conversion.rpy2py(data[1])
    rs_id = robjects.conversion.rpy2py(data[2])
    gwas = robjects.conversion.rpy2py(data[3])

print('sample id length: {}'.format(sample_id.shape[0]))
print('snp dim: {}'.format(snp.shape))
print('snp column matches with gwas? {}'.format(all(rs_id==gwas.SNP.values)))

# replace negative entried in snp matrix with np.nan
snp = snp.astype(np.float32)
snp[snp<0] = np.nan
print('all rs_id equal to gwas SNP ids? {}'.format(all(rs_id==gwas.SNP.values)))
if all(rs_id==gwas.SNP.values):
    print('Start ls imputation')
    # impute missing values for snps
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    snp_imputed = imp_mean.fit_transform(snp)
    
    yhat, _ = imputeY_batch_cholesky(snp_imputed, gwas.BETA.values)
    print('ls imputation done, start saving results')
    
    results = pd.DataFrame({'ID': sample_id, 'pheno': yhat.squeeze()})
    results.to_csv(save_path + 'LSimp_' + impute_data + '_{}_{}_cholesky.csv'.format(batch_size,batch_idx), index=False)
    print('Finished, results saved to ' + save_path + 'LSimp_' + impute_data + '_{}_{}_cholesky.csv'.format(batch_size,batch_idx))
else:
    print('asdfasdfasdf!!!!!!')
    raise ValueError("gwas snps and input bim not aligned!")
   