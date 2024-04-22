import sys
import scipy
import scipy.linalg 
import numpy as np
import os
import pandas as pd
import time
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from keras.utils import to_categorical
from sklearn.impute import SimpleImputer
import csv

dataset = sys.argv[1]

print(dataset)

save_path = '/home/panwei/he000176/jinwen/xgboost/'

r = robjects.r
r['source']('/home/panwei/he000176/jinwen/xgboost/read_snp.R')
read_snp = robjects.globalenv['read_snp']

with localconverter(robjects.default_converter + pandas2ri.converter):
    sample_id, snp, rs_id = read_snp(dataset)

print('snp dim: {}'.format(snp.shape))

snp = snp.astype(np.float32)
snp[snp<0] = np.nan

imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
snp_imputed = imp_mean.fit_transform(snp)

np.save(save_path + dataset + '.npy',snp_imputed)
np.save(save_path + dataset + '_id.npy',sample_id)
np.save(save_path + dataset + '_rs.npy',rs_id)