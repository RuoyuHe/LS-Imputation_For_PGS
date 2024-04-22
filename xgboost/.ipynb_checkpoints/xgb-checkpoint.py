import os
import gc
import sys 
import xgboost as xgb
import matplotlib.pyplot as plt
import pandas as pd, numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from sklearn.impute import SimpleImputer
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split


p = int(sys.argv[1])
datatype = sys.argv[2]

wd = '/home/panwei/he000176/jinwen/xgboost/'
data_path = '/home/panwei/he000176/jinwen/data/'
model_save_path = '/home/panwei/he000176/jinwen/results/xgboost/models/'
results_save_path = '/home/panwei/he000176/jinwen/results/xgboost/'

def get_data(wd,data_path,dataset):
    # ID_imp and ID_real need to be pd.DataFrame 
    snp = np.load(wd + dataset + '.npy')
    sample_id = np.load(wd + dataset + '_id.npy')
    cov = pd.read_csv(data_path + 'pheno/' + dataset + '_cov.csv')
    snp = pd.DataFrame(snp)
    snp['ID'] = sample_id
    X = pd.merge(snp, cov, on = 'ID')
    y = pd.read_csv(data_path + 'pheno/' + dataset + '_pheno.csv')
    all_data = pd.merge(X, y, on='ID')
    if(dataset=='base'):
        lsimp = pd.read_csv(data_path + 'pheno/LSimp_base_7500.csv')
        lsimp.columns = ['ID','LSimp']
        prs = pd.read_csv(data_path + 'pheno/prs_pheno_on_analysis_gln.txt',sep=' ')
        prs.columns = ['ID','prs']
        all_data = pd.merge(all_data,lsimp,on='ID')
        all_data = pd.merge(all_data,prs,on='ID')
    #
#     y = all_data['pheno'].values
#     y_ls = all_data['LSimp'].values
#     y_prs = all_data['pheno'].values
#     X = all_data.drop(columns=['ID','pheno','LSimp','prs'])
    #
    return all_data


######## base ########
ID_imp = pd.read_csv('/home/panwei/he000176/jinwen/data/pheno/training_data_for_base_models/ID_0.{}.txt'.format(p),sep=' ')
ID_imp.columns = ['ID','IID']
ID_imp = ID_imp.drop(columns=['IID'])
if p!=5:
    ID_real = pd.read_csv('/home/panwei/he000176/jinwen/data/pheno/training_data_for_base_models/ID_0.{}.txt'.format(p),sep=' ')
else:
    print('real 0.5B')
    ID_real = pd.read_csv('/home/panwei/he000176/jinwen/data/pheno/training_data_for_base_models/ID_0.{}B.txt'.format(p),sep=' ')

ID_real.columns = ['ID','IID']
ID_real = ID_real.drop(columns=['IID'])
    
base_data = get_data(wd,data_path,'base')
base_imp = pd.merge(base_data, ID_imp, on='ID')
X_base_imp = base_imp.drop(columns=['ID','pheno','LSimp','prs']).values
y_base_ls = base_imp['LSimp'].values
y_base_prs = base_imp['prs'].values

base_real = pd.merge(base_data, ID_real, on='ID')
X_base_real = base_real.drop(columns=['ID','pheno','LSimp','prs']).values
y_base_real = base_real['pheno'].values

X_base_imp_train, X_base_imp_val, \
y_base_ls_train, y_base_ls_val, \
y_base_prs_train, y_base_prs_val = train_test_split(X_base_imp, y_base_ls, y_base_prs, test_size=0.1, random_state=0)

X_base_real_train, X_base_real_val, \
y_base_real_train, y_base_real_val = train_test_split(X_base_real, y_base_real, test_size=0.1, random_state=1)


dtrain_ls = xgb.DMatrix(X_base_imp_train, label=y_base_ls_train)
dval_ls = xgb.DMatrix(X_base_imp_val, label=y_base_ls_val)

dtrain_prs = xgb.DMatrix(X_base_imp_train, label=y_base_prs_train)
dval_prs = xgb.DMatrix(X_base_imp_val, label=y_base_prs_val)

dtrain_real = xgb.DMatrix(X_base_real_train, label=y_base_real_train)
dval_real = xgb.DMatrix(X_base_real_val, label=y_base_real_val)

######## ensemble ########
ensemble_data = get_data(wd,data_path,'ensemble')
ensemble_ID = ensemble_data['ID'].values
X_ensemble = ensemble_data.drop(columns=['ID','pheno']).values
y_ensemble = ensemble_data['pheno'].values

densemble = xgb.DMatrix(X_ensemble, label=y_ensemble)

######## test ########
test_data = get_data(wd,data_path,'test')
test_ID = test_data['ID'].values
X_test = test_data.drop(columns=['ID','pheno']).values
y_test = test_data['pheno'].values

dtest = xgb.DMatrix(X_test, label=y_test)


######## start training ########
params = {
    'eta': 0.1,
    'max_depth': 6,
    'subsample': 0.9,
    'colsample_bytree': 0.5
}

if datatype=='LSimp':
    bst = xgb.train(params, dtrain_ls, num_boost_round=1000, early_stopping_rounds=20, 
                    evals=[(dval_ls, 'eval')])
elif datatype=='prs':
    bst = xgb.train(params, dtrain_prs, num_boost_round=1000, early_stopping_rounds=20, 
                    evals=[(dval_prs, 'eval')])
elif datatype=='real':
    bst = xgb.train(params, dtrain_real, num_boost_round=1000, early_stopping_rounds=20, 
                    evals=[(dval_real, 'eval')])

bst.save_model(model_save_path + datatype + '_0.{}.model'.format(p))

pred_ensemble = bst.predict(densemble)
pred_test = bst.predict(dtest)

pred_ensemble_df = pd.DataFrame({
    'ID':ensemble_ID,
    'pheno':pred_ensemble
})

pred_test_df = pd.DataFrame({
    'ID':test_ID,
    'pheno':pred_test
})

pred_ensemble_df.to_csv(results_save_path + 'ensemble/' + datatype + '_0.{}.csv'.format(p),index=False)
pred_test_df.to_csv(results_save_path + 'test/' + datatype + '_0.{}.csv'.format(p),index=False)
