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

df = df.drop(['Avg_VG_AE','Avg_VG_eQTL'], axis=1)
df = df.drop(['gene_biotype','num_enhancers'], axis=1)
df.dropna(inplace=True)
df.reset_index(drop=True, inplace=True)


df['ensembl_gene_id'] = df['ensembl_gene_id'].astype('string')
df = df.set_index('ensembl_gene_id')

X = df.drop('log_ae', axis = 1)
features = list(X.columns)
y = df['log_ae']

train_X, test_X, train_y, test_y = train_test_split(X, y, test_size = 0.25, random_state = 42)
def GB(X_train, y_train, grid_search = False, specified_params = False, params = dict(), n_cores=4):

   
    steps = [('scaler', StandardScaler()),('regr',GradientBoostingRegressor())]
    pipe = Pipeline(steps)

    pipeline_params =  {
        'regr__loss' : ['squared_error'],
        'regr__max_depth': [3, 5, 6, 10, 15, 20],
        'regr__learning_rate': [0.01, 0.1, 0.2, 0.3],
        'regr__subsample': np.arange(0.5, 1.0, 0.1),
        'regr__min_samples_split': np.arange(0.3, 1.0, 0.1),
        'regr__min_samples_leaf':np.arange(0.3, 1.0, 0.1),
        'regr__max_features':np.arange(0.3, 1.0, 0.2),
        'regr__n_estimators': [100, 500, 1000]    
        }
    
    if grid_search:
        search = GridSearchCV(
            pipe,
            param_grid = pipeline_params,
            scoring = 'neg_root_mean_squared_error',
            n_jobs=n_cores,
            verbose=1
            )
    elif specified_params:
        pipeline_params = params
        search = GridSearchCV(
            pipe,
            param_grid = pipeline_params,
            n_jobs=n_cores,
            scoring = 'neg_root_mean_squared_error',
            verbose=1
            )
    else:
        search =  RandomizedSearchCV(
            pipe,
            n_iter = 25,
            param_distributions= pipeline_params,
            scoring='neg_root_mean_squared_error',
            random_state=42,
            n_jobs=n_cores,
            verbose = 1
            )
    
    search.fit(X_train, y_train)
    return search



search = GB(train_X, train_y,grid_search=True)
search_df = pd.DataFrame(search.cv_results_)
search_df.sort_values(by='rank_test_score')
# search_df.to_csv(r'data/gb_randsearch_grid.csv', sep = '\t',index = False)

importances = list(search.best_estimator_[1].feature_importances_)
importances
imp_df = pd.DataFrame({'features': features,'importance' : importances})
imp_df.sort_values(by='importance',axis=0, ascending=False, inplace= True)
plt.bar(imp_df['features'], imp_df['importance'])
plt.xticks(rotation=90)
plt.xlabel('features')
plt.ylabel('importance')
plt.show()
# plt.savefig(r'data\gb_importance_plot.png', dpi=300)

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
# fig.savefig(r'data\gb_performance.png', dpi = 300)

