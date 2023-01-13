import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from sklearn.utils import resample
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import pearsonr
import xgboost as xgb
from sklearn.metrics import mean_absolute_error

df = pd.read_csv(r'data\imputed_ML_df.csv', sep='\t')

df = df.drop(['gene_biotype'], axis=1)
df.dropna(inplace=True)
df.reset_index(drop=True, inplace=True)


df['ensembl_gene_id'] = df['ensembl_gene_id'].astype('string')
df = df.set_index('ensembl_gene_id')

X = df.drop('log_ae', axis = 1)
features = list(X.columns)
y = df['log_ae']

train_X, test_X, train_y, test_y = train_test_split(X, y, test_size = 0.25, random_state = 42)



def XG(X_train, y_train, grid_search = False, specified_params = False, params = dict(), n_cores=4, n_iter = 2000):

   
    steps = [('scaler', StandardScaler()),('regr',xgb.XGBRegressor(seed=20, verbosity = 2))]
    pipe = Pipeline(steps)

    xgb_params = {
        'regr__max_depth': [3, 5, 6, 10, 15, 20,30],
        'regr__learning_rate': [0.001, 0.01, 0.1, 0.2], # aka eta or epsilon
        'regr__subsample': np.arange(0.5, 1.0, 0.1),
        'regr__colsample_bytree': np.arange(0.4, 1.0, 0.1),
        'regr__colsample_bylevel': np.arange(0.4, 1.0, 0.1),
        'regr__n_estimators': [100, 500, 1000],
        'regr__min_split_loss' : [0, 1, 2, 5],
        'regr__reg_lambda' : [1, 2, 5], # gamma
        'regr__reg_alpha' : [0, 1 , 2, 5]
        }
        # lambda ==> L2 regularization parameter ==> prevent overfitting (default = 1)
        # alpha ==> L1 regularization parameter (default = 0)
        # gamma ==> aka min_split_loss ==> important for pruning  (default = 0)
        # for all of the three parameters: ==> the larger the more conservative and range[0, inf]           
    if grid_search:
        search = GridSearchCV(
            pipe,
            param_grid = xgb_params,
            scoring = 'neg_root_mean_squared_error',
            n_jobs=n_cores,
            verbose=1
            )
    elif specified_params:
        
        search = GridSearchCV(
            pipe,
            param_grid = params,
            n_jobs=n_cores,
            scoring = 'neg_root_mean_squared_error',
            verbose=1
            )
    else:
        search =  RandomizedSearchCV(
            pipe,
            n_iter = n_iter,
            param_distributions= xgb_params,
            scoring='neg_root_mean_squared_error',
            random_state=42,
            n_jobs=n_cores,
            verbose = 1
            )
    
    search.fit(X_train, y_train)
    return search

search = XG(train_X, train_y,n_iter=15000)
search_df = pd.DataFrame(search.cv_results_)
search_df.sort_values(by='rank_test_score')
# search_df.to_csv(r'data\xgb_randomsearch_15k_iter.csv', sep = '\t',index = False)

importances = list(search.best_estimator_[1].feature_importances_)
importances
imp_df = pd.DataFrame({'features': features,'importance' : importances})
imp_df.sort_values(by='importance',axis=0, ascending=False, inplace= True)
plt.bar(imp_df['features'], imp_df['importance'])
plt.xticks(rotation=90)
plt.xlabel('features')
plt.ylabel('importance')
plt.show()
# plt.savefig(r'data\xgb_importance_plot.png', dpi=300)

train_pred = search.predict(train_X)
test_pred = search.predict(test_X)
fig, axs = plt.subplots(1, 2)

p_rho_train, p_p_train = pearsonr(train_pred, train_y)
pearson_train = 'pearson corr: '+str(round(p_rho_train,3)) #+'\n'+'p value: ' +str(p_p_train)
rmse_train = 'RMSE: '+str(round(np.sqrt(mean_absolute_error(train_y, train_pred)), 4))

p_rho_test, p_p_test = pearsonr(test_pred, test_y)
p_p_test = str(p_p_test)[0:3]+str(p_p_test)[-5:]
pearson_test = 'pearson corr: '+str(round(p_rho_test, 3)) #+'\n' +'p value: '+str(p_p_test)
rmse_test = 'RMSE: '+str(round(np.sqrt(mean_absolute_error(test_y, test_pred)),4))

axs[0].scatter(train_pred, train_y, s = 2, alpha = 0.5)
axs[0].add_artist(AnchoredText(rmse_train, loc=2))
axs[0].set_xlabel('log(VG_AE) training')

axs[1].scatter(test_pred, test_y, s = 2, alpha = 0.5)
axs[1].add_artist(AnchoredText(rmse_test, loc=2))
axs[1].set_xlabel('log(VG_AE) test')

fig.supylabel('predicted log(VG_AE)')
plt.show()
# fig.savefig(r'data\xgb_performance.png', dpi = 300)

