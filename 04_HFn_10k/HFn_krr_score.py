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
from cheml import datasets

data_list = [
   datasets.load_HF3(large=True),
   datasets.load_HF4(large=True),
   datasets.load_HF5(large=True),
   datasets.load_HF6(large=True),
]

descriptors = [
    {'n': 1, 'sort':False},
    {'n': -1, 'sort':False, 'nuclear_charges':True},
    {'n': -6, 'sort':False, 'nuclear_charges':True},
]

#gammas = [.01, .02, .05, .1, .2, .5, 1]
#gammas = [1E-12, 1E-9, 1E-6, 1E-3, 1E-1]
gammas = [1E-8, 1E-4, 1E-1]
alphas = [1e-11]
n_samples_list = list(range(10, 100, 10)) + list(range(100, 1000, 100)) + list(range(1000, 3000, 200)) + list(range(3000, 5000, 500)) + list(range(5000, 7500, 1000))

kr_all_scores = []
for i in range(len(data_list)):
    data = data_list[i]
    E = np.array(data['E'])
    cv = ShuffleSplit(len(E), n_iter=10, test_size=.1, random_state=42)
    scores = qtl.krrScore(data, 
               kernels = 'laplacian',
               gammas = gammas,
               alphas = alphas,
               n_samples = n_samples_list,
               cv = cv,
               descriptors = descriptors,
               report=True,
             )
    kr_all_scores.append(scores)

    # step-wise backup
    try:
        tmp = np.stack(kr_all_scores)
        np.savez('HFn_krr_score.npz',
                 score = tmp,
                 n_samples = n_samples_list,
                 gammas = gammas)
    except Exception as e:
        print("save attempt failed with err: %s" % str(e))

kr_all_scores = np.stack(kr_all_scores)
np.savez('st_scores_krr_HFn_lps.npz',
         score = kr_all_scores,
         n_samples = n_samples_list,
         gammas = gammas)
