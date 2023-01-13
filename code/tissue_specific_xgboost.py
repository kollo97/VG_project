import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.colors as colors
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
from sklearn.metrics import mean_squared_error
from matplotlib_venn import venn2

## load all necessary files
df = pd.read_csv(r'data\imputed_ML_df.csv', sep='\t')
df = df.drop(['log_ae','log_eqtl','median_tpm','gene_biotype'], axis = 1)

ae = pd.read_excel(r'data\VG_AE.xlsx', sheet_name=0, header=1, skiprows=0)
eqtl= pd.read_excel(r'data\VG_eqtl.xlsx', sheet_name=0, header=1, skiprows=0)
expression_data = pd.read_csv(r'data\GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct', skiprows=2, sep='\t')

tissue_list = set(ae.columns[1:len(ae.columns)]).intersection(set(eqtl.columns[1:len(eqtl.columns)]))
tissue_list = list(tissue_list) # important for iteration later!
print(len(tissue_list))



ae['ensembl_gene_id'] = ae['IDs']
ae.set_index('ensembl_gene_id', inplace=True)
ae.drop('IDs', axis=1, inplace=True)

cols = list(ae.columns)
for ind, col in enumerate(cols):
    new_col = 'log_ae_'+str(col)
    cols[ind] = new_col
    # print(new_col)
ae = np.log(ae)
ae.columns = cols

ae.reset_index(inplace=True)




eqtl['ensembl_gene_id'] = eqtl['IDs']
eqtl.set_index('ensembl_gene_id', inplace=True)
eqtl.drop('IDs', axis=1, inplace=True)

cols = list(eqtl.columns)
for ind, col in enumerate(cols):
    new_col = 'log_eqtl_'+str(col)
    cols[ind] = new_col
eqtl = np.log(eqtl)
eqtl.columns = cols
eqtl.reset_index(inplace=True)


expression_data.head()
expression_data['ensembl_gene_id'] = expression_data['gene_id']
expression_data.drop(['gene_id','Description'], axis=1, inplace=True)
expression_data['ensembl_gene_id'] = expression_data['ensembl_gene_id'].apply(lambda x: str(x).split('.')[0])
expression_data.set_index('ensembl_gene_id', inplace=True)
tissues = pd.read_excel(r'data\VG_AE.xlsx', sheet_name=4, header=1, skiprows=0)

tissues_dict = {}
for i, r in tissues.iterrows():
    key = r.TISSUE_NAME_Original
    value = r.TISSUE_ABBRV
    tissues_dict[key] = value
cols = list(expression_data.columns)
for ind, col in enumerate(cols):
    if col in tissues_dict.keys():
        cols[ind] = tissues_dict[col]
expression_data.columns = cols

cols = list(expression_data.columns)
for ind, col in enumerate(cols):
    new_col = 'median_tpm_'+str(col)
    cols[ind] = new_col
expression_data.columns = cols
expression_data.reset_index(inplace=True)

## get best params from RandomizedSearchCV dataframe
def extract_params(search_df_path, sep = '\t'):
    search_df = pd.read_csv(search_df_path, sep = sep)
    params =  search_df.loc[search_df['rank_test_score']==1,'params'].values[0] #take the row of the best performing model and extract the 'params' column (it's a dict) 
    params = eval(str(params))
    for key, value in params.items(): # put the single values into a list, because the GridSearchCV needs an iterable. I use gridsearch because it performs 5-fold cross-validation on a selected validation set by default 
        params[key] = [value]
    return params

params = extract_params(r'data\xgb_randomsearch_20k_iter.csv', sep = '\t')

def XG(X_train, y_train, grid_search = False, specified_params = False, params=dict(), n_cores=4, n_iter = 2000):

   
    steps = [('scaler', StandardScaler()),('regr',xgb.XGBRegressor(seed=20, verbosity = 1))]
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
        xgb_params = params
        search = GridSearchCV(
            pipe,
            param_grid = xgb_params,
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

# summarized lists of following for loop
#maybe just use the dict and the importnaces_df
bar_plot_metric = []
bar_plot_eqtl = []
performances_xgboost = {'tissue':[],'rmse_train':[],'rmse_test':[],'rmse_eqtl':[],'pearson_train': [],'pearson_test':[],'pearson_eqtl':[],'num_ae_values':[],'num_predicted_values':[],'num_eqtl_predicted_values':[]}
print(performances_xgboost.keys())
performances_eqtl = []
importances_df = pd.DataFrame(columns=['features','importance','tissue'])


for tissue in tissue_list:
    print(tissue)
    ae_tissue = ae.columns[ae.columns.str.endswith(tissue)].values[0]
    eqtl_tissue = eqtl.columns[eqtl.columns.str.endswith(tissue)].values[0]
    median_tpm_tissue = expression_data.columns[expression_data.columns.str.endswith(tissue)].values[0]
    df1 = pd.merge(df,eqtl[['ensembl_gene_id',eqtl_tissue]],on='ensembl_gene_id', how='left')
    df1 = pd.merge(df1,expression_data[['ensembl_gene_id',median_tpm_tissue]], on='ensembl_gene_id', how='left')
    df1 = df1.dropna()
    df1 = pd.merge(df1, ae[['ensembl_gene_id', ae_tissue]], on = 'ensembl_gene_id', how='left')

    ml_df = df1.drop('ensembl_gene_id', axis=1)
    ml_df = ml_df.dropna()
    ml_df = ml_df.reset_index(drop=True)

    X = ml_df.drop(ae_tissue, axis = 1)
    # X = X.drop('log_eqtl', axis=1)
    features = list(X.columns)
    y = ml_df[ae_tissue]
    train_X, test_X, train_y, test_y = train_test_split(X, y, test_size = 0.25, random_state = 42)
    print('Training Features Shape:', train_X.shape)
    print('Training Labels Shape:', train_y.shape)
    print('Testing Features Shape:', test_X.shape)
    print('Testing Labels Shape:', test_y.shape)

## Venn diagram
    metrics = df1.drop(ae_tissue,axis=1).dropna()
    metrics = metrics.loc[:,'ensembl_gene_id']
    fig, ax = plt.subplots(1,1)
    # aes = df1.loc[:,['ensembl_gene_id',ae_tissue]].dropna()
    # aes = aes.loc[:,'ensembl_gene_id']
    aes=    set(ae[ae[ae_tissue].notna()]['ensembl_gene_id'])
    venn = venn2([set(metrics), set(aes)],['metrics','AE data'], ax =ax)
    plt.title(tissue)
    
    # fig.savefig(r'tissue_specific\venn_metric\{}_venn_metrics.png'.format(tissue), dpi=300)
    
    bar_plot_metric.append(metrics)
    bar_plot_metric.append(aes)


## perform XGBoost
    search = XG(train_X, train_y, specified_params =True, params=params)


## feature importance plot
    fig, ax = plt.subplots(1,1)
    importances = list(search.best_estimator_[1].feature_importances_)
    
    imp_df = pd.DataFrame({'features': features,'importance' : importances})
    imp_df.sort_values(by='importance',axis=0, ascending=False, inplace= True)
    ax.bar(imp_df['features'], imp_df['importance'])
    plt.xticks(rotation=90)
    ax.set_xlabel('features')
    ax.set_ylabel('importance')

    # plt.show()
    # fig.savefig(r'tissue_specific\feature_importance\{}_feature_importance.png'.format(tissue), dpi=300)
    tissue_series = (tissue +' ')*len(X.columns)
    tissue_series = tissue_series.split(' ')
    tissue_series = pd.Series(tissue_series) 
    imp_df['tissue'] = tissue_series
    importances_df = pd.concat([importances_df,imp_df])

## model performance
    train_pred = search.predict(train_X)
    test_pred = search.predict(test_X)
    fig, axs = plt.subplots(1, 2)
    p_rho_train, p_p_train = pearsonr(train_pred, train_y)
    pearson_train = 'pearson corr: '+str(round(p_rho_train,3)) #+'\n'+'p value: ' +str(p_p_train)
    rmse_train = 'RMSE: '+str(round(np.sqrt(mean_squared_error(train_y, train_pred)), 4))

    p_rho_test, p_p_test = pearsonr(test_pred, test_y)
    p_p_test = str(p_p_test)[0:3]+str(p_p_test)[-5:]
    pearson_test = 'pearson corr: '+str(round(p_rho_test, 3)) #+'\n' +'p value: '+str(p_p_test)
    rmse_test = 'RMSE: '+str(round(np.sqrt(mean_squared_error(test_y, test_pred)),4))

    axs[0].scatter(train_pred, train_y, s = 2, alpha = 0.5)
    axs[0].add_artist(AnchoredText(rmse_train, loc=2))
    axs[0].set_xlabel('log(VG_AE) training')


    axs[1].scatter(test_pred, test_y, s = 2, alpha = 0.5)
    axs[1].add_artist(AnchoredText(rmse_test, loc=2))
    axs[1].set_xlabel('log(VG_AE) test')

    fig.supylabel('predicted log(VG_AE)')
    # plt.show()
    # fig.savefig(r'tissue_specific\model_performance\xgb\{}_xgb_performance.png'.format(tissue), dpi = 300)

    

## comparison to only eQTL
    fig, ax = plt.subplots(1,1)
    rho, p = pearsonr(ml_df[eqtl_tissue], ml_df[ ae_tissue])
    rho_text = 'pearson corr: ' + str(round(rho,4))
    rmse = round(np.sqrt(mean_squared_error(ml_df[eqtl_tissue], ml_df[ae_tissue])),4)
    rmse = 'RMSE: '+ str(rmse)
    ax.scatter(ml_df[eqtl_tissue], ml_df[ ae_tissue])
    ax.add_artist(AnchoredText(rmse, loc=2))
    ax.set_xlabel('log(VG_AE)')
    ax.set_ylabel('log(VG_eQTL)')
    # fig.savefig(r'tissue_specific\model_performance\eqtl\{}_eqtl_correlation.png'.format(tissue), dpi = 300)

## eQTL venn diagram
    fig, ax = plt.subplots(1,1)
    eqtls = set(eqtl[eqtl[eqtl_tissue].notna()]['ensembl_gene_id'])
    # set2 = set(ae[ae[ae_tissue].notna()]['ensembl_gene_id']) ==
    venn2([eqtls, aes], ['eqtl_data','ae_data'], ax=ax)
    plt.title(tissue)
    # fig.savefig(r'tissue_specific\venn_eqtl\{}_venn_eQTL.png'.format(tissue), dpi=300)

    performances_xgboost['rmse_train'].append(np.sqrt(mean_squared_error(train_y, train_pred)))
    performances_xgboost['rmse_test'].append(np.sqrt(mean_squared_error(test_y, test_pred)))
    performances_xgboost['pearson_train'].append(p_rho_train)
    performances_xgboost['pearson_test'].append(p_rho_test)
    performances_xgboost['tissue'].append(str(tissue))
    performances_xgboost['num_ae_values'].append(len(aes))
    performances_xgboost['num_predicted_values'].append(len(aes.union(metrics)))
    performances_xgboost['num_eqtl_predicted_values'].append(len(aes.union(eqtls)))
    performances_xgboost['rmse_eqtl'].append(np.sqrt(mean_squared_error(ml_df[eqtl_tissue], ml_df[ae_tissue])))
    performances_xgboost['pearson_eqtl'].append(rho)

## wirte prediction to csv
    y_predicted_complete = df1['ensembl_gene_id'] #for df1 the na values were already dropped 
    # make a series of all the ensembl gene ids that we have metrics for, including those that we have no AE_VG for 
    X_complete = df1.drop('ensembl_gene_id', axis=1)
    X_complete = X_complete.drop(ae_tissue, axis=1)
    vg_values = pd.Series(search.predict(X_complete), name='log_ae_vg')
    y_predicted_complete = pd.concat([y_predicted_complete, vg_values], axis=1) # concatenate the ensembl_gene_ids with the predicted log VG_AE values

    # y_predicted_complete.to_csv(r'tissue_specific\predicted_vg\only_predicted\{}_pred_AE_VG.tsv'.format(tissue), sep = '\t',index = False)

    plt.close('all')

performances_df = pd.DataFrame.from_dict(performances_xgboost)
performances_df.to_csv(r'data\summary_df.tsv', sep='\t', index=False)
importances_df.to_csv(r'data\importances_df.tsv',sep='\t', index=False)