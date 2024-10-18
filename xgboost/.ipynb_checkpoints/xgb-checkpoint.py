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
import argparse
import logging

# Setup argument parser
parser = argparse.ArgumentParser(description='Process some impute data.')
parser.add_argument('--datatype', type=str, required=True, help='LS-imputed trait, prs-imputed trait, or observed trait?')
parser.add_argument('--p', type=str, required=True, help='number of folds used as the imputed data')
parser.add_argument('--xgb_data_path', type=str, required=True, help='Path to the xgb .npy genotype data')
parser.add_argument('--pheno_path', type=str, required=True, help='path to phenotype')
parser.add_argument('--model_path', type=str, required=True, help='path to save model')
parser.add_argument('--out_path', type=str, required=True, help='path to save results')
parser.add_argument('--log_file', type=str, required=True, help='path to log file for debugging')

# Parse arguments
args = parser.parse_args()
datatype = args.datatype
p = args.p
wd = args.xgb_data_path
data_path = args.pheno_path
model_save_path = args.model_path
results_save_path = args.out_path
log_file = args.log_file


logging.basicConfig(
    filename = log_file,  # File to write logs to
    level = logging.DEBUG,  # Level of logging
    format = '%(asctime)s - %(levelname)s - %(message)s',  # Format of logs
    datefmt = '%Y-%m-%d %H:%M:%S'  # Date format in logs
)



def get_data(wd,data_path,dataset,p):
    # ID_imp and ID_real need to be pd.DataFrame 
    snp = np.load(wd + '/' + dataset + '.npy')
    sample_id = np.load(wd + '/' + dataset + '_id.npy')
    # cov = pd.read_csv(data_path + 'pheno/' + dataset + '_cov.csv')
    snp = pd.DataFrame(snp)
    snp['ID'] = sample_id
    # X = pd.merge(snp, cov, on = 'ID')
    X = snp.copy()
    if(dataset=='base'):
        if p==5:
            y = pd.read_csv(data_path + '/training_data_for_base_models/pheno_real_0.{}B_adjusted.csv'.format(p))
        else:
            y = pd.read_csv(data_path + '/training_data_for_base_models/pheno_real_0.{}_adjusted.csv'.format(p))
        data_real = pd.merge(X, y, on='ID')
        lsimp = pd.read_csv(data_path + '/training_data_for_base_models/pheno_LSimp_0.{}.csv'.format(p))
        lsimp.columns = ['ID','LSimp']
        prs = pd.read_csv(data_path + '/training_data_for_base_models/pheno_prs_0.{}.csv'.format(p))
        prs.columns = ['ID','prs']
        data_imp = pd.merge(X,lsimp,on='ID')
        data_imp = pd.merge(data_imp,prs,on='ID')
    else:
        y = pd.read_csv(data_path + dataset + '_pheno_adjusted.csv')
        data_real = pd.merge(X, y, on='ID')
        data_imp = None
    #
    return data_real, data_imp


######## training data ########
logging.info("Reading base data")
base_real, base_imp = get_data(wd,data_path,'base',p)

X_base_imp = base_imp.drop(columns=['ID','LSimp','prs']).values
y_base_ls = base_imp['LSimp'].values
y_base_prs = base_imp['prs'].values

X_base_real = base_real.drop(columns=['ID','pheno']).values
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
logging.info("Reading ensemble data")
ensemble_data, _ = get_data(wd,data_path,'ensemble',p)
ensemble_ID = ensemble_data['ID'].values
X_ensemble = ensemble_data.drop(columns=['ID','pheno']).values
y_ensemble = ensemble_data['pheno'].values

densemble = xgb.DMatrix(X_ensemble, label=y_ensemble)

######## test ########
logging.info("Reading test data")
test_data, _ = get_data(wd,data_path,'test',p)
test_ID = test_data['ID'].values
X_test = test_data.drop(columns=['ID','pheno']).values
y_test = test_data['pheno'].values

dtest = xgb.DMatrix(X_test, label=y_test)


######## start training ########
params = {
    'eta': 0.1,
    'max_depth': 6,
    'subsample': 0.9,
    'colsample_bytree': 0.5,
    'colsample_bylevel':0.7, 
    'colsample_bynode':0.7,
    'alpha' : 20,
    'lambda' : 50
}

logging.info("Train base model")
if datatype=='LSimp':
    bst = xgb.train(params, dtrain_ls, num_boost_round=1000, early_stopping_rounds=20, 
                    evals=[(dval_ls, 'eval')])
elif datatype=='prs':
    bst = xgb.train(params, dtrain_prs, num_boost_round=1000, early_stopping_rounds=20, 
                    evals=[(dval_prs, 'eval')])
elif datatype=='real':
    bst = xgb.train(params, dtrain_real, num_boost_round=1000, early_stopping_rounds=20, 
                    evals=[(dval_real, 'eval')])



bst.save_model(model_save_path + '/' + datatype + '_0.{}.model'.format(p))
logging.info("Base model saved to " + model_save_path + '/' + datatype + '_0.{}.model'.format(p)) 


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

pred_ensemble_df.to_csv(results_save_path + '/ensemble/' + datatype + '_0.{}.csv'.format(p),index=False)
logging.info("Predicted values saved to " + results_save_path + '/ensemble/' + datatype + '_0.{}.csv'.format(p)) 
pred_test_df.to_csv(results_save_path + '/test/' + datatype + '_0.{}.csv'.format(p),index=False)
logging.info("Predicted values saved to " + results_save_path + '/test/' + datatype + '_0.{}.csv'.format(p)) 

logging.info('test R^2 {}'.format(r2_score(test_data.pheno.values, pred_test)))