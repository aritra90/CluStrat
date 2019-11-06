from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import pandas as pd
import scipy, time, math
from sklearn import metrics

import _utils

def se(X1, y, ridgeinv, mse, Sig, sketch_flag): 

    n, m = X1.shape 
    
    if int(sketch_flag) == 1:
    	########## SKETCH X
        time0 = time.time()
    	# O(n) sketch dimension
        s = 2*X1.shape[0]
        #S = np.random.normal(0.0, 1.0, (X1.shape[1],s))/math.sqrt(s)
        #C = np.matmul(X1,S)
        C = countSketch(X1,s)
        print('sketch size: ', C.shape)
        print('Sketching time (mins): ', (time.time()-time0)/60.0)
    else:
        C = X1
  
    ###########################################################
    print(C.shape)
    coeff = np.zeros(m)
    
    for i in range(m):
        b1 = ridgeinv.dot(X1[:,i])
        coeff[i] = np.linalg.norm(b1)**2
   
    #print("\n -----  Sig  -----\n")
    #print(Sig)
    compdof = lambda x: x/(x + np.median(Sig))
    #print("\n ----- effdof ----- \n") 
    
    effdof = list(map(compdof, Sig))
    #print(effdof)
    
    est_mse = mse/(n - ((1.25)*sum(effdof)+0.5))
    #print(est_mse) 
    #print("\n ----- ########## -----\n")
    se_coeff = np.sqrt(est_mse*coeff)
    
    return se_coeff

def coef_tval(X, y, ridgeinv, beta, mse, Sig, sketch_flag):
    """Calculate t-statistic for beta coefficients.

    Parameters
    ----------
   
    X : numpy.ndarray
        Training data used to fit the classifier.
    y : numpy.ndarray
        Target training values, of shape = [n_samples].

    Returns
    -------
    numpy.ndarray
        An array of t-statistic values.
    """
    #a = np.array(clf.intercept_ / coef_se(clf, X, y, alpha)[0])
    #b = np.array(clf.coef_ / coef_se(clf, X, y, alpha)[1:])
    denom = se(X, y, ridgeinv, mse, Sig, sketch_flag) 
    #print(denom) 
    print(denom[:10])
    print(beta[:10])
    a = np.array(beta[0]/denom[0])
    b = np.array(beta[1:] / denom[1:])
    #print(a) 
    #print(b) 	
 	
    return np.append(a,b)


def coef_pval(X, y, ridgeinv, beta, mse, Sig, sketch_flag):
    """Calculate p-values for beta coefficients.

    Parameters
    ----------
    clf : sklearn.linear_model
        A scikit-learn linear model classifier with a `predict()` method.
    X : numpy.ndarray
        Training data used to fit the classifier.
    y : numpy.ndarray
        Target training values, of shape = [n_samples].

    Returns
    -------
    numpy.ndarray
        An array of p-values.
    """
    n = X.shape[0]
    t = coef_tval(X, y, ridgeinv, beta, mse, Sig, sketch_flag)
    
    p  = 2 * (1 - scipy.stats.t.cdf(abs(t), n-1))
    #print(*p, sep='\t')
    #print(np.nansum(p))
    return p

