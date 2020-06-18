#######################################################################
## Wrapper for running CluStrat. Please refer to the Readme to see how
## to run CluStrat/
##
## The code performs a structure informed clustering based stratification
## technique to detect causal SNPs

### Run:
### For simulated data
### python or python3 CluStrat_wrapper.py --sim/-s 1
###				-model BN/PSD/TGP -m 100 -n 5000 -tf 1 or 0 -pf 1 or 0 -prop
###
### 1 or 0 for trait flags whether you want binary or continuous traits
########################################################################
## Authors:
##            Aritra Bose, IBM Research, Yorktown Heights, NY
##            Myson Burch, Purdue University, West Lafayette, IN
## Contact:
##            a.bose@ibm.com; mcburch@purdue.edu
#######################################################################
#########################       IMPORT      ###########################
#######################################################################
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy import *
import scipy as sp
import sys, csv, array, os, math, subprocess, itertools, time, random, warnings
warnings.filterwarnings(action="ignore")
import fastcluster, normalize, getSE
import scipy.stats
from sklearn.model_selection import KFold
from scipy.sparse.linalg import svds
import scipy.cluster.hierarchy as sch
from sklearn import linear_model, svm
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.linear_model import RidgeCV, Ridge, LassoCV
from sklearn.model_selection import train_test_split, GridSearchCV
import sklearn.metrics as metrics
from sklearn.preprocessing import QuantileTransformer, quantile_transform
from scipy.stats import t
from sklearn.kernel_ridge import KernelRidge
# import regressors as rgrs
from regressors import stats
from sklearn.feature_selection import chi2
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression, LogisticRegression


def intersect(*d):
    sets = iter(map(set, d))
    result = sets.next()
    for s in sets:
        result = result.intersection(s)
    return result

def scipy_cluster(data, method):
    clusters = sch.linkage(data, method=method, metric='euclidean')
    return clusters

def fast_cluster(data, method):
    clusters = fastcluster.linkage(data, method=method, metric='euclidean', preserve_input=False)
    return clusters

def ridge_pvals(X, Y, sketch_flag):
####################################################################
    #lambdas = np.logspace(-2.5,2.5,25)
    Y = Y - np.mean(Y)
    m, n = X.shape
    print(X.shape)
    k = 10
    if k >= min(m,n):
       _,Sig,_ = np.linalg.svd(X.dot(X.T))
    else:
       if m < n:
          _,Sig,_ = svds(X.dot(X.T), k)
       else:
          _,Sig,_ = svds((X.T).dot(X), k)
    #Sig[::-1].sort()
    medSig = np.median(Sig)
    maxSig = max(Sig)
    #alphas = np.append(np.append(Sig[:3], medSig), min(Sig))
    #regr = RidgeCV(alphas = Sig, cv = 5)
    #regr.fit(X,Y)
    #y_pred = regr.predict(X)
    X1 = np.hstack((np.ones((m, 1)), np.matrix(X)))
    ridgeinv = np.linalg.inv(X1.dot(X1.T) + medSig*np.eye(m)) # medSig*np.eye(m)
    # print("Ridge inverse diagonal: ",ridgeinv.diagonal()[:10])
    betahat = ((X1.T).dot(ridgeinv)).dot(Y)
    y_pred = X1.dot(betahat.T)
    # print("true y (mu centered): ", Y[Y.shape[0]-25:])
    # print("y pred: ", y_pred[y_pred.shape[0]-25:])
    betahat = betahat.flatten().tolist()[0]
    mse = metrics.mean_squared_error(y_pred, Y)

    print("\n=============== \n")
    #print("\n Lambda is : "+str(regr.alpha_))
    #print("\n R2 score : "+str(regr.score(X,Y)))
    print("MSE : " + str(mse))
    print("Median Eigenvalue : " + str(medSig))
    print("Max Eigenvalue : " + str(maxSig))
    print("\n =============== \n")

    pvals = getSE.coef_pval(X1, Y, ridgeinv, betahat, mse, Sig, sketch_flag)
    print('CluStrat pvals: ',pvals)
    print('mean: ',np.nanmean(pvals))
    print('min: ',min(pvals))
    print('max: ',max(pvals))

    return pvals


def chi2_pvals(X,Y):
    print('\n\n\n\n\n')
    print("ASSOCIATION TEST OUTPUT:")
    print('Dimension of X: ',X.shape)
    
    scores, pvals = chi2(X,Y)
    print(pvals)
    
    return pvals
    
def linear_pvals(X, Y, sketch_flag, lmbda=0):
####################################################################
    #lambdas = np.logspace(-2.5,2.5,25)
    print('\n\n\n\n\n')
    print("LINEAR REGRESSION OUTPUT:")
    Y = Y - np.mean(Y)
    m, n = X.shape
    print('Dimension of X: ',X.shape)

    lm = LinearRegression()
    #lm = LogisticRegression(random_state=0, penalty='l2')
    lm.fit(X,Y)
    params = np.asarray(lm.coef_)
    print('Linear regression estimated coeffs: ', params)
    predictions = lm.predict(X)
    print('Original signal: ',Y)
    print('Linear regression predictions: ',predictions)
    # predictions = predictions + np.random.uniform(min(predictions),max(predictions),predictions.shape[0])/2
    # print(predictions)
    MSE = metrics.mean_squared_error(predictions, Y) 
    #print('DEBUGGING')
    #MSE = 10**16*MSE
    print('Mean-squared error: ',MSE)
    #print('DEBUGGING')
    #sys.exit(0)
    temp = np.linalg.pinv(np.dot(X.T,X)).diagonal()
    print('Diagonal elements of pseudo inverse: ',temp)
    var_b = MSE*(temp)
    print('var_b: ', var_b)
    sd_b = np.sqrt(var_b)
    print('sd_b: ', var_b)
    ts_b = params/ sd_b
    print('ts_b: ', ts_b)
    p_values =np.asarray([2*(1-t.cdf(np.abs(i),(m-1))) for i in ts_b])
    print('Linear regression pvalues: ', p_values)
    dbg_nnzidx = np.nonzero(p_values)[0]
    print('pvalue mean (0s/nans excluded): ',np.nanmean(p_values[dbg_nnzidx]))
    print('pvalue min (0s/nans excluded): ',np.nanmin(p_values[dbg_nnzidx]))
    print('pvalue max (0s/nans excluded): ',np.nanmax(p_values[dbg_nnzidx]))
    print('Length of pvalues: ', p_values.shape[0])
    print('Number of zeros in pvalues: ', np.sum(p_values == 0))
    print('Number of nans in pvalues: ', np.isnan(p_values).sum())


    # print('No. of significant SNPs by LinearRegression: ',len([i for i in dbg_nnzidx if p_values[i] < 0.0001]))
    # print('mean: ',np.nanmean(p_values[dbg_nnzidx]))
    # print('min: ',np.nanmin(p_values[dbg_nnzidx]))
    # print('max: ',np.nanmax(p_values[dbg_nnzidx]))

    # print('\nend of DEBUGGING')

    return p_values #pvals

def cluster(R, D, status, depth, pvalue, numclust, sketch_flag, verbose, ids_and_chroms=None):

    clustcount = np.zeros((len(numclust),1))
    CS = np.zeros((len(numclust),1))
    SP = np.zeros((len(numclust),1))
    allidx = []
    Z = scipy_cluster(D,'ward')
    for k in range(0,len(numclust)):

        den2 = sch.dendrogram(Z, leaf_rotation=90.,orientation='left', leaf_font_size=2)
        plt.xlabel('Distance')
        plt.savefig('Z_dendogram_'+str(k)+'.png')
        #sys.exit(0)
	
        plt.clf()
        plt.cla()
        plt.close()
        ClustMem = sch.fcluster(Z, numclust[k], depth=depth)

        oldidx = []
        allpvals = []
        newpval = []
        intidx = []
        sti = []
        combset = []
        m, n = R.shape
        if verbose == "1":
            print('The number of clusters is : ', np.unique(ClustMem).shape[0]);
        clustcount[k] = np.unique(ClustMem).shape[0]
        ct = 0

        # for each cluster...
        for i in range(0, np.unique(ClustMem).shape[0]):
            print('cluster ID is: ',i)
            Ri = R[ClustMem == i+1,:]
            normRi,fr = normalize.norm(Ri,0)
            #normRi = Ri
            sti = status[ClustMem == i+1]
            pvalc = np.zeros(n)
            ####################################################################
            # get weights of all SNPs using Ridge regression (output are SNP ids?)
            pvals = ridge_pvals(normRi,sti, sketch_flag)
            ridge_pval = np.array(pvals[1:])
            ridge_pval[isnan(ridge_pval)] = 0

            if sum(ridge_pval) != 0:
                nnzidx = np.nonzero(ridge_pval)[0]
                pvidx = [i for i in nnzidx if ridge_pval[i] < (pvalue)]
                # if verbose == "1":
                #     print("number of significant SNPs : " + str(len(pvidx)))
            else:
                if verbose == "1":
                    print("\n All Ridge pvals are ZERO")
                pvidx = np.zeros(normRi.shape[1])
            ####################################################################

            if ct == 0:
                finalpvals_pt1 = ridge_pval
            else:
                if verbose == "1":
                    print("Computing minimum pvalues from previous iteration...")
                for i in range(finalpvals_pt1.shape[0]):
                    if finalpvals_pt1[i] != 0 and ridge_pval[i] != 0:
                        finalpvals_pt1[i] = min(finalpvals_pt1[i], ridge_pval[i])
                    else:
                        finalpvals_pt1[i] = max(finalpvals_pt1[i], ridge_pval[i])

            if sum(ridge_pval) != 0:
                if len(oldidx) == 0:
                  oldidx = pvidx
                  ct = ct + 1
                elif ct == 1:
                  intidx = np.intersect1d(oldidx,pvidx)
                  combset = np.unique(np.append(oldidx,pvidx))
                  ct = ct + 1
                else:
                  comb_intidx = np.intersect1d(combset,pvidx)
                  #Update set of significant intersections
                  intidx = np.unique(np.append(intidx,comb_intidx))
                  combset = list(np.unique(np.append(combset,pvidx)))

            # print('(running) Number of significant SNPs: ',len(combset))

        # print('final pvalues shape (pt1)...')
        # print(finalpvals_pt1.shape)
        # np.savetxt("CluStrat_pvals_pt1_"+str(k)+".txt", finalpvals_pt1)
        # np.savetxt("CluStrat_idx_pt1_"+str(k)+".txt", combset)

        #Select the significant SNPs from the original data
        Rpv = R[:,list(np.array(intidx).astype(int))]
        #Normalize the data
        normRpv,_ = normalize.norm(Rpv,0)
        #normRpv = Rpv 
# ################################################################################
# ################################################################################
#         print("DEBUGGING")
#         regressor_OLS = sm.OLS(status,normRpv).fit()
#         debug_pvals = np.asarray(regressor_OLS.pvalues.astype(float))
#         print(debug_pvals)
#         debug_pvals[isnan(debug_pvals)] = 0
#         dbg_nnzidx = np.nonzero(debug_pvals)[0]
#         dbg_pvidx = [i for i in dbg_nnzidx if debug_pvals[i] < pvalue]
#         print(len(dbg_pvidx))
#         np.savetxt("debugging_pvals.out", debug_pvals[dbg_pvidx])
# ################################################################################
# ################################################################################

        #Rpvals = linear_pvals(normRpv, status, sketch_flag, 0)
        Rpvals = ridge_pvals(normRpv, status, sketch_flag)
        Rpvals = np.array(Rpvals[1:])
        Rpvals[isnan(Rpvals)] = 0
        print('P-value threshold: ', pvalue)
        if sum(ridge_pval) != 0:
            nnzidx = np.nonzero(Rpvals)[0]
            Rpvidx = [i for i in nnzidx if Rpvals[i] < pvalue]
        else:
            if verbose == "1":
                print("All pvals are ZERO\n")
            #Store Junk
            Rpvidx = [-9]
        Rpvidx = np.array(Rpvidx).astype(int)
        if Rpvidx.all() != -9:
            sig3 = len(np.where(intidx[Rpvidx].astype(int) < 10)[0])
            SP[k] = abs(len(Rpvidx) - sig3)
            CS[k] = sig3
        else:
            SP[k] = 0
            CS[k] = 0

        # print("Number of causal SNPs in CluStrat = " + str(CS[k]))
        print("Number of final associations " + str(SP[k]))

        # Grab and convert final significant SNP indices
        Rpvidx = np.asarray(Rpvidx.astype(int))
        combset = np.asarray(combset).astype(int)
        final_idx = np.asarray(combset[Rpvidx].astype(int))
        # final_idx = np.asarray(final_idx.astype(int))
        # pvals = finalpvals_pt1[final_idx]
        # np.savetxt("CluStrat_pvals_final_"+str(k)+".txt", pvals)
        # np.savetxt("CluStrat_idx_final_"+str(k)+".txt", final_idx)

        ############################ Write Output ############################
        # Writing final output (SNP IDs, chromosome label and p-value)
        if ids_and_chroms is None:
            pass
        else:
            chromids = np.asarray(ids_and_chroms[1])
            SNPids = np.asarray(ids_and_chroms[0])

            data_frame = pd.DataFrame(chromids[final_idx],columns=['chrom'])
            data_frame['SNPs'] = pd.Series(SNPids[final_idx], index=data_frame.index)
            # VERSION 3 of GITHUB CODE FIXING PVAL BUG
            # data_frame['p-values'] = pd.Series(finalpvals_pt1[final_idx], index=data_frame.index)
            data_frame['p-values'] = pd.Series(Rpvals[Rpvidx], index=data_frame.index)
            data_frame = data_frame.sort_values(by='p-values')
            data_frame.to_csv('CluStrat_assoc_numclust'+str(k)+'.txt', index = False, sep = ' ', header=True)

    return  CS, clustcount, SP#, finalpvals_pt1, final_idx
