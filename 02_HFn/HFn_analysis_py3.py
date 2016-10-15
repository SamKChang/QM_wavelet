
# coding: utf-8

# In[ ]:

#import qctoolkit as qtk
#import qctoolkit.projects.Basel.p14_stml.scatteringTransform as qst
#import qctoolkit.projects.Basel.p14_stml.estimators as qes

import sys
import os
script_path = os.path.join(os.path.split(os.getcwd())[0], 'script')
sys.path.insert(0, script_path)
import scatteringTransform as qst
import estimators as qes
import tools as qtl
import pickle

import numpy as np
from sklearn.cross_validation import ShuffleSplit
from sklearn.cross_validation import cross_val_score
from sklearn.linear_model import Ridge
from sklearn.kernel_ridge import KernelRidge
import time
import glob

def pload(fname):
    with open(fname, 'rb') as f:
        u = pickle._Unpickler(f)
        u.encoding = 'latin1'
        data = u.load()
        return data
#        return pickle.load(f)


# # Construct ST regression matrix and OLS path

# In[ ]:

signal_setting = {
    'padding': 10,
    'grid_step': 0.01,
    'components': ['all', 'valence', 'core']
}

filter_setting = {
    'derivatives': np.array([0, 1, 2]),
    'scales': np.arange(3, 10),
    #'derivatives': np.array([0,1]),
    #'scales': np.arange(6,8),
    'wavelet': qst.gabor_k,
}
fname_list = sorted(glob.glob('../data/data_m?.pkl'))


# In[ ]:

st_matrix_list = []
E_list = []

for fname in fname_list:
    print("processing %s" % fname)

    ti = time.time()
    st_matrix, E = qst.stModel_1d(
        fname, 
        batch = 5,
        signal_setting = signal_setting, 
        filter_setting = filter_setting
    )
    tf = time.time()
    print("time: %f" % (tf - ti))
    st_matrix_list.append(st_matrix)
    E_list.append(E)

np.savez('st_models_HFn.npz', st=st_matrix_list, E=E_list)


# In[ ]:

#st_model = np.load('../data/st_models_HFn.npz')
#st_matrix_list = st_model['st']
#E_list = st_model['E']


# In[ ]:

mae_list = []
mse_list = []
components_list = []
for i in range(len(fname_list)):
    E = np.array(E_list[i])
    stm = st_matrix_list[i]
    cv = ShuffleSplit(len(E), n_iter=10, test_size=.1, random_state=42)
    MSE, MAE, components = qes.ols_path_cv(stm, E.ravel(), cv=cv, ridge_penalty=1e-6, max_num_atoms=500, return_components=True)
    mse_list.append(MSE)
    mae_list.append(MAE)
    components_list.append(components)
    
#np.savez('st_models_HFn_ols.npz', maes=mae_list, mses=mse_list, components_list = components_list)


# # ST learning curves

# In[ ]:

#st_data = np.load('../data/st_models_HFn.npz')
#ols_data = np.load('../data/st_models_HFn_ols.npz')
data_list = [pload(f) for f in sorted(glob.glob('../data/data_m?.pkl'))]
st_matrix_list = st_data['st']
E_list = st_data['E']
components_list = ols_data['components_list']

n_samples_list = range(10,100,10) + range(100, 1000, 100)
alphas = [1e-6]
n_components_list = [1, 2, 5, 10, 20, 50, 100, 200, 300, 500]

def getData_st(i):
    return data_list[i], st_matrix_list[i], components_list[i]
def getData_kr(i):
    return data_list[i]


# In[ ]:

st_scores = []
for i in range(len(data_list)):
    data, st_matrix, components = getData_st(i)
    E = np.array(data['E'])
    cv = ShuffleSplit(len(E), n_iter=10, test_size=.1, random_state=42)
    stsetting = {
        'regression_matrix': st_matrix,
        'n_samples_list': n_samples_list,
        'alphas': alphas,
        'n_components_list': n_components_list,
        'ols_components': components,
        'cv': cv,
        'threads': 10,
    }
    st_scores.append(qst.stScore(data, **stsetting))
    
st_scores = np.stack(st_scores)
np.savez('st_scores_HFn.npz',
         score = st_scores,
         n_samples = n_samples_list,
         components = n_components_list)


# In[ ]:

descriptors = {
    'distance': {'n': 1, 'sort':False},
    'distance_2': {'n': 2, 'sort':False},
    'coulomb_nosort_nocharge': {'n': -1, 'sort':False},
    'coulomb_nocharge': {'n': -1},
    'coulomb_nosort': {'n': -1, 'sort':False, 'nuclear_charges':True},
    'coulomb': {'n': -1, 'nuclear_charges':True},
    'coulomb_inv2': {'n': -2, 'nuclear_charges':True},
}


# In[ ]:

gammas = [.01, .02, .05, .1, .2, .5, 1]
alphas = [1e-6, 1e-9, 1e-12]
n_components_list = [1, 2, 5, 10, 20, 50, 100, 200, 300, 500]
n_samples_list = list(range(10, 100, 10)) +     list(range(100, 1000, 100)) +     list(range(1000, 3000, 200)) +     list(range(3000, 5000, 500)) +     list(range(5000, 7500, 1000))


# In[ ]:

kr_all_scores = []
for i in range(len(data_list)):
    kr_scores = []
    kr_all_scores.append(kr_scores)
    data = getData_kr(i)
    E = np.array(data['E'])
    cv = ShuffleSplit(len(E), n_iter=10, test_size=.1, random_state=42)
    krrsetting = {
        'kernel': 'rbf',
        'gammas': gammas,
        'alphas': alphas,
        'n_samples_list': n_samples_list,
        'cv': cv,
    }
    for name, descriptor in descriptors.iteritems():
        if 'nuclear_charges' in descriptor and descriptor['nuclear_charges'] is not None:
            descriptor['nuclear_charges'] = data['Z']
        krrsetting['descriptor_setting'] = descriptor
        print(krrsetting)
        scores = qtl.krrScore(data, **krrsetting)
        kr_scores.append(scores)
kr_all_scores = np.stack(kr_all_scores)
np.savez('st_scores_krr_rbf_HFn.npz',
         score = st_scores,
         n_samples = n_samples_list,
         gammas = gammas)

