import pandas as pd
import numpy as np
# from matplotlib.offsetbox import AnchoredText
from sklearn.utils import resample
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV


df = pd.read_csv(r'data\imputed_ML_df.csv', sep='\t')

df = df.drop(['Avg_VG_AE','Avg_VG_eQTL'], axis=1)
df = df.drop(['gene_biotype','num_enhancers'], axis=1)
df.dropna(inplace=True)
df.reset_index(drop=True, inplace=True)


df['ensembl_gene_id'] = df['ensembl_gene_id'].astype('string')
df = df.set_index('ensembl_gene_id')

X = df.drop('log_ae', axis = 1)
features = list(X.columns)
# X = X.drop('log_eQTL', axis=1)
y = df['log_ae']
# print(df.head())

train_X, test_X, train_y, test_y = train_test_split(X, y, test_size = 0.25, random_state = 42)
train_X, valid_X, train_y, valid_y = train_test_split(train_X, train_y, test_size = 0.25, random_state = 42)
# print('\n','Training Features Shape:', train_X.shape)
# print('Training Labels Shape:', train_y.shape)
# print('Validation Features Shape:', valid_X.shape)
# print('Validation Labels Shape:', valid_y.shape)
# print('Testing Features Shape:', test_X.shape)
# print('Testing Labels Shape:', test_y.shape)
# print(sum([train_X.shape[0],valid_X.shape[0], test_X.shape[0]])==len(df))


def SV(X_train, y_train, grid_search = False, specified_params = False, params = dict(), n_cores=4):

    steps = [('scaler', StandardScaler()),('regr',SVR())]
    pipe = Pipeline(steps)

    pipeline_params = {
        'svm__C': 1*10**np.arange(0, 5.0),
        'svm__gamma':np.arange(0.0005,0.5,0.05),
        'svm__epsilon': 1*10**np.arange(-4, 4.0) ,
        'svm__kernel': ['rbf'] #linear and poly are more options
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

search = SV(train_X, train_y, grid_search=True)

search_df = pd.DataFrame(search.cv_results_)
search_df.sort_values(by='rank_test_score')
search_df.to_csv(r'data\svm_gridsearch_uppm.csv', sep = '\t',index = False)