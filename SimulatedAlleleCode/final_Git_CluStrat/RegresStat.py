# -*- coding: utf-8 -*-

"""This module contains functions for calculating various statistics and
coefficients."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import pandas as pd
import scipy, time, math
from sklearn import metrics
from sklearn.decomposition import PCA

import _utils

def countSketch(matrixA, s):
    m, n = matrixA.shape
    matrixC = np.zeros([m, s])
    hashedIndices = np.random.choice(s, n, replace=True)
    randSigns = np.random.choice(2, n, replace=True) * 2 - 1 # a n-by-1{+1, -1} vector
    matrixA = matrixA * randSigns.reshape(n,1) # flip the signs of 50% columns of A
    for i in range(s):
        idx = (hashedIndices == i)
        matrixC[i,:] = np.sum(matrixA[:,idx], 1)
    return matrixC


def residuals(clf, X, y, r_type='standardized'):
    """Calculate residuals or standardized residuals.

    Parameters
    ----------
    clf : sklearn.linear_model
        A scikit-learn linear model classifier with a `predict()` method.
    X : numpy.ndarray
        Training data used to fit the classifier.
    y : numpy.ndarray
        Target training values, of shape = [n_samples].
    r_type : str
        Type of residuals to return: 'raw', 'standardized', 'studentized'.
        Defaults to 'standardized'.

        * 'raw' will return the raw residuals.
        * 'standardized' will return the standardized residuals, also known as
          internally studentized residuals, which is calculated as the residuals
          divided by the square root of MSE (or the STD of the residuals).
        * 'studentized' will return the externally studentized residuals, which
          is calculated as the raw residuals divided by sqrt(LOO-MSE * (1 -
          leverage_score)).

    Returns
    -------
    numpy.ndarray
        An array of residuals.
    """
    # Make sure value of parameter 'r_type' is one we recognize
    assert r_type in ('raw', 'standardized', 'studentized'), (
        "Invalid option for 'r_type': {0}".format(r_type))
    y_true = y.view(dtype='float')
    # Use classifier to make predictions
    y_pred = clf.predict(X)
    # Make sure dimensions agree (Numpy still allows subtraction if they don't)
    assert y_true.shape == y_pred.shape, (
        "Dimensions of y_true {0} do not match y_pred {1}".format(y_true.shape,
                                                                  y_pred.shape))
    # Get raw residuals, or standardized or standardized residuals
    resids = y_pred - y_true
    if r_type == 'standardized':
        resids = resids / np.std(resids)
    elif r_type == 'studentized':
        # Prepare a blank array to hold studentized residuals
        studentized_resids = np.zeros(y_true.shape[0], dtype='float')
        # Calcluate hat matrix of X values so you can get leverage scores
        hat_matrix = np.dot(
            np.dot(X, np.linalg.inv(np.dot(np.transpose(X), X))),
            np.transpose(X))
        # For each point, calculate studentized residuals w/ leave-one-out MSE
        for i in range(y_true.shape[0]):
            # Make a mask so you can calculate leave-one-out MSE
            mask = np.ones(y_true.shape[0], dtype='bool')
            mask[i] = 0
            loo_mse = np.average(resids[mask] ** 2, axis=0)  # Leave-one-out MSE
            # Calculate studentized residuals
            studentized_resids[i] = resids[i] / np.sqrt(
                loo_mse * (1 - hat_matrix[i, i]))
        resids = studentized_resids
    return resids


def sse(clf, X, y):
    """Calculate the standard squared error of the model.

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
    float
        The standard squared error of the model.
    """
    y_hat = clf.predict(X)
    sse = np.sum((y_hat - y) ** 2)
    return sse


def adj_r2_score(clf, X, y):
    """Calculate the adjusted :math:`R^2` of the model.

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
    float
        The adjusted :math:`R^2` of the model.
    """
    n = X.shape[0]  # Number of observations
    p = X.shape[1]  # Number of features
    r_squared = metrics.r2_score(y, clf.predict(X))
    return 1 - (1 - r_squared) * ((n - 1) / (n - p - 1))


def coef_se(clf, X, y, alpha, Sig, sketch_flag):
    """Calculate standard error for beta coefficients.

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
        An array of standard errors for the beta coefficients.
    """
    n = X.shape[0]

    #print("MSE is ")
    #print(metrics.mean_squared_error(y, clf.predict(X)))
    #print("Solution Beta")

    X1 = np.hstack((np.ones((n, 1)), np.matrix(X)))
    m = X1.shape[1]
    #print(X1.shape)

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
    # print(C.shape)
    coeff = np.zeros(m)
    ridgeinv = np.linalg.inv(C.dot(C.T) + alpha*np.eye(n))
    for i in range(m):
        b1 = ridgeinv.dot(X1[:,i])
        coeff[i] = np.linalg.norm(b1)**2

    #effdof = sum([i/(i + alpha) for i in Sig])
    effdof = n - (1.25)*np.trace((X1.T).dot(ridgeinv.dot(X1))) + 0.5
    #print(effdof)
    est_mse = (metrics.mean_squared_error(y,clf.predict(X)))/effdof
    # print(est_mse)
    se_coeff = np.sqrt(est_mse*coeff)

    return se_coeff


def coef_tval(clf, X, y, alpha, Sig, sketch_flag):
    """Calculate t-statistic for beta coefficients.

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
        An array of t-statistic values.
    """
    #a = np.array(clf.intercept_ / coef_se(clf, X, y, alpha)[0])
    #b = np.array(clf.coef_ / coef_se(clf, X, y, alpha)[1:])
    denom = coef_se(clf, X, y, alpha, Sig, sketch_flag)
    #print(denom)
    a = np.array(clf.intercept_ / denom[0])
    b = np.array(clf.coef_ / denom[1:])
    #print(a)
    #print(b)

    return np.append(a, b), denom


def coef_pval(clf, X, y, alpha, Sig, sketch_flag):
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
    t = coef_tval(clf, X, y, alpha, Sig, sketch_flag)
    p, SE = 2 * (1 - scipy.stats.t.cdf(abs(t), n-1))
    print('p values: ', p)
    #print(*p, sep='\t')
    #print(np.nansum(p))
    return p, SE


def f_stat(clf, X, y):
    """Calculate summary F-statistic for beta coefficients.

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
    float
        The F-statistic value.
    """
    n = X.shape[0]
    p = X.shape[1]
    r_squared = metrics.r2_score(y, clf.predict(X))
    return (r_squared / p) / ((1 - r_squared) / (n - p - 1))


def summary(clf, X, y, xlabels=None):
    """
    Output summary statistics for a fitted regression model.

    Parameters
    ----------
    clf : sklearn.linear_model
        A scikit-learn linear model classifier with a `predict()` method.
    X : numpy.ndarray
        Training data used to fit the classifier.
    y : numpy.ndarray
        Target training values, of shape = [n_samples].
    xlabels : list, tuple
        The labels for the predictors.
    """
    # Check and/or make xlabels
    ncols = X.shape[1]
    if xlabels is None:
        xlabels = np.array(
            ['x{0}'.format(i) for i in range(1, ncols + 1)], dtype='str')
    elif isinstance(xlabels, (tuple, list)):
        xlabels = np.array(xlabels, dtype='str')
    # Make sure dims of xlabels matches dims of X
    if xlabels.shape[0] != ncols:
        raise AssertionError(
            "Dimension of xlabels {0} does not match "
            "X {1}.".format(xlabels.shape, X.shape))
    # Create data frame of coefficient estimates and associated stats
    coef_df = pd.DataFrame(
        index=['_intercept'] + list(xlabels),
        columns=['Estimate', 'Std. Error', 't value', 'p value']
    )
    coef_df['Estimate'] = np.concatenate(
        (np.round(np.array([clf.intercept_]), 6), np.round((clf.coef_), 6)))
    coef_df['Std. Error'] = np.round(coef_se(clf, X, y), 6)
    coef_df['t value'] = np.round(coef_tval(clf, X, y), 4)
    coef_df['p value'] = np.round(coef_pval(clf, X, y), 6)
    # Create data frame to summarize residuals
    resids = residuals(clf, X, y, r_type='raw')
    resids_df = pd.DataFrame({
        'Min': pd.Series(np.round(resids.min(), 4)),
        '1Q': pd.Series(np.round(np.percentile(resids, q=25), 4)),
        'Median': pd.Series(np.round(np.median(resids), 4)),
        '3Q': pd.Series(np.round(np.percentile(resids, q=75), 4)),
        'Max': pd.Series(np.round(resids.max(), 4)),
    }, columns=['Min', '1Q', 'Median', '3Q', 'Max'])
    # Output results
    print("Residuals:")
    print(resids_df.to_string(index=False))
    print('\n')
    print('Coefficients:')
    print(coef_df.to_string(index=True))
    print('---')
    print('R-squared:  {0:.5f},    Adjusted R-squared:  {1:.5f}'.format(
        metrics.r2_score(y, clf.predict(X)), adj_r2_score(clf, X, y)))
    print('F-statistic: {0:.2f} on {1} features'.format(
        f_stat(clf, X, y), ncols))
