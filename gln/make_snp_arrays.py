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
import csv

dataset = sys.argv[1]

print(dataset)

save_path = '/home/panwei/he000176/jinwen/data/arrays/'

r = robjects.r
r['source']('/home/panwei/he000176/jinwen/gln/read_snp.R')
read_snp = robjects.globalenv['read_snp']

with localconverter(robjects.default_converter + pandas2ri.converter):
    sample_id, snp, rs_id = read_snp(dataset)

print('snp dim: {}'.format(snp.shape))
    
n=snp.shape[0]
for i in range(n):
    print(i)
    snp1=snp[i]
    snp1[np.where(snp1<0)[0]]=3
    snp1=to_categorical(snp1).T
    np.save(save_path + dataset + '/{}'.format(sample_id[i]) +'.npy',snp1)