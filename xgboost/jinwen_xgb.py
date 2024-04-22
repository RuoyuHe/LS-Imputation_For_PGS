import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV
import random
import gc
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt


from sklearn.linear_model import LinearRegression
from scipy import stats
from sklearn.ensemble import RandomForestRegressor

model_save_path = '/home/panwei/he000176/jinwen/results/xgboost/jinwen_rep/models/'
results_save_path = '/home/panwei/he000176/jinwen/results/xgboost/jinwen_rep/'

wba_unrelated = pd.read_csv('/home/panwei/he000176/jinwen/data/pheno/WBA_unrelated_ID.txt',sep=' ')

########## prepare ensemble and test data ##########
# snp=np.load('/home/panwei/fu000217/dl/LipA/7_apply/realdata/all_snp_matrix/test/test_snp_final.npy')
# covariates=pd.read_csv('/home/panwei/fu000217/dl/LipA/7_apply/realdata/extract_covariates/test/covariates.csv')
# snp = pd.DataFrame(snp)
# y_val=pd.read_csv('/home/panwei/fu000217/dl/LipA/7_apply/realdata/phenos/test/pheno_real_test.csv')
# y_val=y_val.drop(columns='ID')

snp = pd.read_csv('/home/panwei/fu000217/dl/LipA/7_apply/realdata/all_snp_matrix/unrelated/test/snp.csv')
covariates = pd.read_csv('/home/panwei/fu000217/dl/LipA/7_apply/realdata/all_snp_matrix/unrelated/test/cov.csv')
y_val = pd.read_csv('/home/panwei/fu000217/dl/LipA/7_apply/realdata/all_snp_matrix/unrelated/test/pheno_real.csv')
print('IDs aligned? {}'.format(all(snp.ID.values==covariates.ID.values) and all(snp.ID.values==y_val.ID.values)))
y_val = y_val.drop(columns='ID')
snp = snp.drop(columns='ID')
covariates = covariates.drop(columns='ID')


##### KEEP ONLY WBA UNRELATED INDIVIDUALS
# mask = covariates['ID'].isin(wba_unrelated['ID'])
# keep_idx = mask.to_numpy().nonzero()[0]
# covariates=covariates.drop(columns='ID')
# covariates = covariates.iloc[keep_idx,:]
# y_val = y_val.iloc[keep_idx,:]
# snp = snp.iloc[keep_idx,:]

X = pd.concat([snp,covariates],axis=1)
X_val_train, X_val_test, y_train_true, y_test_true = train_test_split(X, y_val, test_size=0.2, random_state=42)
sample_weights_val_train = np.ones(X_val_train.shape[0])
sample_weights_val_test = np.ones(X_val_test.shape[0])
densemble = xgb.DMatrix(X_val_train, label=y_train_true, weight=sample_weights_val_train)
dr2test = xgb.DMatrix(X_val_test, label=y_test_true, weight=sample_weights_val_test)

# collect garbage
del snp, X, y_val
del covariates
gc.collect()

########## prepare training data ##########
# snp = np.load('/home/panwei/fu000217/dl/LipA/7_apply/realdata/all_snp_matrix/train/train_snp_final.npy')
# covariates = pd.read_csv('/home/panwei/fu000217/dl/LipA/7_apply/realdata/extract_covariates/train/covariates.csv')
# snp = pd.DataFrame(snp)
# Y = pd.read_csv('/home/panwei/fu000217/dl/LipA/7_apply/realdata/phenos/training/pheno_real_train.csv')
# Y_impute=pd.read_csv('/home/panwei/fu000217/dl/LipA/7_apply/realdata/phenos/training/pheno_pred_train.csv')

snp = pd.read_csv('/home/panwei/fu000217/dl/LipA/7_apply/realdata/all_snp_matrix/unrelated/train/snp.csv')
covariates = pd.read_csv('/home/panwei/fu000217/dl/LipA/7_apply/realdata/all_snp_matrix/unrelated/train/cov.csv')
Y = pd.read_csv('/home/panwei/fu000217/dl/LipA/7_apply/realdata/all_snp_matrix/unrelated/train/pheno_real.csv')
Y_impute = pd.read_csv('/home/panwei/fu000217/dl/LipA/7_apply/realdata/all_snp_matrix/unrelated/train/pheno_pred.csv')
print('IDs aligned? {}'.format(all(snp.ID.values==covariates.ID.values) and all(snp.ID.values==Y.ID.values) and all(snp.ID.values==Y_impute.ID.values)))
Y = Y.drop(columns='ID')
Y_impute = Y_impute.drop(columns='ID')
snp = snp.drop(columns='ID')
covariates = covariates.drop(columns='ID')

##### KEEP ONLY WBA UNRELATED INDIVIDUALS
# mask = covariates['ID'].isin(wba_unrelated['ID'])
# keep_idx = mask.to_numpy().nonzero()[0]
# covariates = covariates.drop(columns='ID')
# covariates = covariates.iloc[keep_idx,:]
# Y = Y.iloc[keep_idx,:]
# Y_impute = Y_impute.iloc[keep_idx,:]
# snp = snp.iloc[keep_idx,:]

random.seed(42)
# use_index = random.sample(range(covariates.shape[0]),(covariates.shape[0])//5)
use_index = random.sample(range(covariates.shape[0]),int(covariates.shape[0]//2.15))
snp = snp.iloc[use_index]
covariates = covariates.iloc[use_index]
X = pd.concat([snp,covariates],axis=1)

# all_ids=Y["ID"]
Y = Y.iloc[use_index]
Y_impute = Y_impute.iloc[use_index]


del snp
del covariates
gc.collect()



########## start training ##########
p=0.4
random.seed(42)
real_ind=random.sample(range(len(use_index)),int(len(use_index)*p))

######
#train with real traits
X_train, X_test, y_train, y_test = train_test_split(X.iloc[real_ind], Y.iloc[real_ind], test_size=0.2, random_state=42)
sample_weights_train=np.ones(X_train.shape[0])
sample_weights_test=np.ones(X_test.shape[0])
dtrain = xgb.DMatrix(X_train, label=y_train, weight=sample_weights_train)
dtest = xgb.DMatrix(X_test, label=y_test, weight=sample_weights_test)
print('here')
def custom_loss(y_pred, dtrain):
    y_true = dtrain.get_label()
    sample_weight = dtrain.get_weight()
    # Gradient
    grad = -2 * sample_weight * (y_true - y_pred)
    # Hessian (constant in this case)
    hess = 2 * sample_weight * np.ones_like(y_true)
    return grad, hess

def custom_eval_metric(y_pred, dtrain):
    y_true = dtrain.get_label()
    sample_weight = dtrain.get_weight()
    # Gradient
    custom_metric_value = np.sum(sample_weight * (y_true - y_pred)**2)
    return 'CustomMetricName', custom_metric_value
params = {
    'eta': 0.1,
    'max_depth': 6,
    'subsample': 0.9,
    'colsample_bytree': 0.9,
}

# Train with custom objective
bst_real = xgb.train(params, dtrain, num_boost_round=1000, early_stopping_rounds=20, 
                evals=[(dtest, 'eval')], obj=custom_loss, feval=custom_eval_metric)
bst_real.save_model(model_save_path + f'real_{p}_unrelated.model')
y_pred = bst_real.predict(dr2test)
np.savetxt(results_save_path + f'test/real_{p}_unrelated.txt', y_pred)
np.savetxt(results_save_path + f'ensemble/real_{p}_unrelated.txt', bst_real.predict(densemble))
print('R2 of predicted real vs observed: {}'.format(r2_score(dr2test.get_label(),y_pred)))
del dtrain
del dtest, X_train, X_test, y_train, y_test
gc.collect()

######
# train with imputed traits
impute_ind=list(set(range(len(use_index)))-set(real_ind))
X_train, X_test, y_train, y_test = train_test_split(X.iloc[impute_ind], Y_impute.iloc[impute_ind], test_size=0.2, random_state=42)
sample_weights_train=np.ones(X_train.shape[0])
sample_weights_test=np.ones(X_test.shape[0])
dtrain = xgb.DMatrix(X_train, label=y_train, weight=sample_weights_train)
dtest = xgb.DMatrix(X_test, label=y_test, weight=sample_weights_test)
print('here')
bst_impute = xgb.train(params, dtrain, num_boost_round=1000, early_stopping_rounds=20, 
                evals=[(dtest, 'eval')], obj=custom_loss, feval=custom_eval_metric)
bst_impute.save_model(model_save_path + f'LSimp_{1-p}_unrelated.model')
y_pred = bst_impute.predict(dr2test)
np.savetxt(results_save_path + f'test/LSimp_{1-p}_unrelated.txt', y_pred)
np.savetxt(results_save_path + f'ensemble/LSimp_{1-p}_unrelated.txt', bst_impute.predict(densemble))
print('R2 of predicted imputed vs observed: {}'.format(r2_score(dr2test.get_label(),y_pred)))
del dtrain
del dtest, X_train, X_test, y_train, y_test
gc.collect()

#calculate combined R^2
y_train_real=bst_real.predict(densemble)
y_test_real=bst_real.predict(dr2test)
y_train_impute=bst_impute.predict(densemble)
y_test_impute=bst_impute.predict(dr2test)

regressor_train=np.array([y_train_real,y_train_impute]).T
regressor_test=np.array([y_test_real,y_test_impute]).T
regressor_train1, regressor_val,y_train_true1, y_val_true=train_test_split(regressor_train,y_train_true,test_size=0.2,random_state=42)

#add start here
model = LinearRegression()

model.fit(regressor_train1, y_train_true1)

y_pred_comb=model.predict(regressor_test).flatten()
print('R2 of ensemble model:',r2_score(y_test_true,y_pred_comb))
print('R2 of single model:',r2_score(y_test_true,y_test_real))
print('improvement:',(r2_score(y_test_true,y_pred_comb)-r2_score(y_test_true,y_test_real))/r2_score(y_test_true,y_test_real))

residuals = np.array(y_test_true).flatten() - y_pred_comb
stderr_residuals = np.std(residuals)
X_with_const = np.column_stack([np.ones(regressor_train.shape[0]), regressor_train])
stderr_coeff = stderr_residuals * np.sqrt(np.linalg.inv(np.dot(X_with_const.T, X_with_const)).diagonal())

# Compute the t-values
t_values = model.coef_ / stderr_coeff[1:]

# Compute the p-values
p_values = (2 * (1 - stats.t.cdf(np.abs(t_values), df=len(X) - X.shape[1] - 1)))
print('pvalue:',p_values)
print('coefficients:',model.coef_)
print('intercept:',model.intercept_)
#add end here
param_grid = {
    'learning_rate': [0.01, 0.1, 0.3],
    'max_depth': [3, 5, 7],
    'subsample': [0.5, 0.7, 1],
    'colsample_bytree': [0.5, 0.7, 1],
    'gamma': [0, 0.1, 0.2],
    'objective': ['reg:squarederror']  # Objective based on the problem type
}

model = xgb.XGBRegressor()
grid_search = GridSearchCV(model, param_grid, cv=3, scoring='neg_mean_squared_error', verbose=1)

# Train the model
grid_search.fit(regressor_val, y_val_true)
best_params = grid_search.best_params_
print(best_params)
#r2_score(y_test_true,model.predict(regressor_test))

model = xgb.XGBRegressor(**best_params)
model.fit(regressor_train1,y_train_true1)
print('combined R2 with XGB:',r2_score(y_test_true,model.predict(regressor_test)))
print('improvement:',(r2_score(y_test_true,model.predict(regressor_test))-r2_score(y_test_true,y_test_real))/r2_score(y_test_true,y_test_real))
# Plot feature importance based on F-score (weight)
xgb.plot_importance(model, importance_type='weight', title='Feature Importance based on Weight')
plt.show()
