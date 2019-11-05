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
import RegresStat as RS
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

import ArmitChisq

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
    ridgeinv = np.linalg.inv(X1.dot(X1.T) + medSig*np.eye(m))
    betahat = ((X1.T).dot(ridgeinv)).dot(Y)
    y_pred = X1.dot(betahat.T)
    betahat = betahat.flatten().tolist()[0]
    mse = metrics.mean_squared_error(y_pred, Y)
    print("\n=============== \n")
    #print("\n Lambda is : "+str(regr.alpha_))
    #print("\n R2 score : "+str(regr.score(X,Y)))
    print("\n MSE : " + str(mse))
    print("\n Median Eigenvalue : " + str(medSig))
    print("\n Max Eigenvalue : " + str(maxSig))
    print("\n =============== \n")


    pvals = getSE.coef_pval(X1, Y, ridgeinv, betahat, mse, Sig, sketch_flag)
    #SE = getSE.se(X1, Y, ridgeinv, mse, Sig, sketch_flag)
    print(pvals)


    return pvals




def cluster(R, D, pops, status, pvalue, dele, sketch_flag):

    clustcount = np.zeros((len(dele),1))
    CS = np.zeros((len(dele),1))
    SP = np.zeros((len(dele),1))
    allidx = []
    Z = scipy_cluster(D,'ward')
    for k in range(0,len(dele)):
        numclust = pops+dele[k]


        den2 = sch.dendrogram(Z, leaf_rotation=90.,orientation='left', leaf_font_size=2)
        #plt.title('Dendrogram for PSD model \{$\alpha=0.1$\}')
        plt.xlabel('Distance')
        # plt.ylabel('Distance')
        plt.savefig('Z_dendogram_'+str(k)+'.png')#,format='png', dpi=1000)

        # sys.exit(0)

        plt.clf()
        plt.cla()
        plt.close()
        # dendrogram(Z,'Orientation','left','ColorThreshold','default');
        ClustMem = sch.fcluster(Z, numclust, depth=85000)

        oldidx = []
        allpvals = []
        newpval = []
        intidx = []
        sti = []
        combset = []
        m, n = R.shape
        print('The number of clusters is : ', np.unique(ClustMem).shape[0]);
        clustcount[k] = np.unique(ClustMem).shape[0]
        # print(clustcount)
        # print(' ')
        ct = 0
        #se_clust = np.zeros(np.unique(ClustMem).shape[0])
        # for each cluster...
        for i in range(0, np.unique(ClustMem).shape[0]):
            print('cluster ID is: ',i)
            Ri = R[ClustMem == i+1,:]
            normRi,fr = normalize.norm(Ri,0)
            #normRi = Ri
            sti = status[ClustMem == i+1]
            pvalc = np.zeros(n)
            ####################################################################
            # get wieghts of all SNPs using Ridge regression (output are SNP ids?)
            pvals = ridge_pvals(normRi,sti, sketch_flag)
            ridge_pval = np.array(pvals[1:])
            ridge_pval[isnan(ridge_pval)] = 0


            if sum(ridge_pval) != 0:
                nnzidx = np.nonzero(ridge_pval)[0]

                pvidx = [i for i in nnzidx if ridge_pval[i] < (pvalue)]

                print("number of significant SNPs : " + str(len(pvidx)))
            else:
                print("\n All Ridge pvals are ZERO")
                pvidx = np.zeros(normRi.shape[1],1)
            ####################################################################

            if ct == 0:
                # old_nnzidx = nnzidx
                # finalpvals_pt1 = ridge_pvals
                finalpvals_pt1 = np.zeros(ridge_pval.shape)
            # else:
            #     print("Computing minimum pvalues...")
            #     # problem here is the zeros...
            #     finalpvals_pt1 = np.minimum(finalpvals_pt1, ridge_pvals)

            if sum(ridge_pval) != 0:
                if len(oldidx) == 0:
                  oldidx = pvidx
                  ct = ct + 1
                elif ct == 1:
                  intidx = np.intersect1d(oldidx,pvidx)
                  combset = np.unique(np.append(oldidx,pvidx))

                  ############ Grabbing final pvals (this code just updates the new ones coming in)
                  if len(intidx) > 0:
                      repeatidx = np.intersect1d(intidx,combset)
                      if len(repeatidx) > 0:
                          for elem in repeatidx:
                              finalpvals_pt1[elem] = np.min(finalpvals_pt1[elem],ridge_pval[elem])
                  #########################

                  ct = ct + 1
                else:
                  comb_intidx = np.intersect1d(combset,pvidx)
                  #Update set of significant intersections
                  intidx = np.unique(np.append(intidx,comb_intidx))
                  combset = list(np.unique(np.append(combset,pvidx)))

                  ############ Grabbing final pvals
                  if len(intidx) > 0:
                      repeatidx = np.intersect1d(comb_intidx,combset)
                      if len(repeatidx) > 0:
                          for elem in repeatidx:
                              finalpvals_pt1[elem] = np.min(finalpvals_pt1[elem],ridge_pval[elem])
                  #########################

        print('final pvalues shape (pt1)...')
        print(finalpvals_pt1.shape)
        np.savetxt("CluStrat_pvals_pt1_"+str(k)+".txt", finalpvals_pt1)
        np.savetxt("CluStrat_idx_pt1_"+str(k)+".txt", combset)

        #Select the significant SNPs from the original data
        Rpv = R[:,intidx.astype(int)]
        #Normalize the data
        normRpv,_ = normalize.norm(Rpv,0)
        Rpvals = ridge_pvals(normRpv, status, sketch_flag)
        Rpvals = np.array(Rpvals[1:])
        Rpvals[isnan(Rpvals)] = 0
        if sum(ridge_pval) != 0:
            nnzidx = np.nonzero(Rpvals)[0]
            Rpvidx = [i for i in nnzidx if Rpvals[i] < pvalue]
        else:
            print("All pvals are ZERO\n")
            #Store Junk
            Rpvidx = [-9]
        Rpvidx = np.array(Rpvidx)
        if Rpvidx.all() != -9:
            sig3 = len(np.where(intidx[Rpvidx].astype(int) < 10)[0])
            SP[k] = abs(len(Rpvidx) - sig3)
            CS[k] = sig3
        else:
            SP[k] = 0
            CS[k] = 0

        # finalpvals_pt2 = Rpvals
        print('shape of ridge pvalues (pt2)...')
        print(Rpvals.shape)
        print('shape of significant SNPs...')
        print(Rpvidx.shape)
        np.savetxt("CluStrat_pvals_pt2_"+str(k)+".txt", Rpvals)
        np.savetxt("CluStrat_idx_pt2_"+str(k)+".txt", Rpvidx)

        print("number of causal SNPs in CluStrat = " + str(CS[k]))
        print("Number of final associations " + str(SP[k]))

    return  CS, clustcount, SP, Rpvals, Rpvidx

#



######################################################################
################# JUNK
######################################################################


     # if (len(sti) > 1):
                # for j in range(0,normR.shape[1]):
                    # pvalc[j] = ArmitChisq.chisq(Ri[:,j], sti.T, len(sti))

            # print(pvalc)
            #sys.exit(0)

            #pvidx = pvalc[np.where(pvalc < pvalue)[0]]
            #chisqst = t.ppf(1-pvalue, df =1 )
            #pvidx = np.where(ridge_pvals < pvalue)[0]
            # 'nan's coming from 'corrcoef' since the stddev is 0
            # print(pvalc[np.where(np.isnan(pvalc) == True)[0]])
            # these also occur in the Matlab data but is ignored (no warnings thrown)
            #print(pvidx)
            #nnzidx = nnzidx[0].astype(int)
            #print(pvidx)

            #print(nnzidx)
            #print(ridge_pval[nnzidx[pvidx]])


            #print(combset)


        # print(intidx)


# def ridge_reg_pvals(X, status):
# ####################################################################
    # #X1 = np.hstack((np.ones((n, 1)), np.matrix(X)))
    # #X_tr, X_te, y_tr, y_te = train_test_split(X, status, test_size=0.3,
    # #                                        random_state=0)
    # alphas = [0.0001,0.001,0.01,0.05,0.1,1.0,5.0,10.0,50.0,100.0]
    # #X = normalize.norm(X,0)

    # # print(X_tr.shape)
    # # print(X_te.shape)
    # # X_tr = normalize.norm(X_tr,0)
    # # X_te = normalize.norm(X_te,0)
    # # y_tr = y_tr - np.mean(y_tr)
    # # y_te = y_te - np.mean(y_te)


    # regr = RidgeCV(alphas=alphas, cv = 5)
    # regr.fit(X_tr, y_tr)
    # y_pred = regr.predict(X_te)

    # print("\n Train R2 score: " + str(regr.score(X_tr,y_tr)))
    # print("\n Test R2 score: " + str(regr.score(X_te,y_te)))
    # print("\n MAE score: " + str(metrics.median_absolute_error(y_te, y_pred)))

    # # regr_trans = TransformedTargetRegressor(regressor=RidgeCV(alphas=alphas, cv = 5),
                                            # # transformer=QuantileTransformer(output_distribution='normal'))

    # # regr_trans.fit(X_tr, y_tr)
    # # y_pred = regr_trans.predict(X_te)
    # # ax1.scatter(y_te, y_pred)
    # # ax1.plot([0, 2000], [0, 2000], '--k')
    # # ax1.set_ylabel('Target predicted')
    # # ax1.set_xlabel('True Target')
    # # ax1.set_title('Ridge regression \n with target transformation')
    # # ax1.text(100, 1750, r'$R^2$=%.2f, MAE=%.2f' % (
    # # metrics.r2_score(y_te, y_pred), metrics.median_absolute_error(y_te, y_pred)))
    # # ax1.set_xlim([0, 2000])
    # # ax1.set_ylim([0, 2000])

    # # f.suptitle("Synthetic data", y=0.035)
    # # f.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])

    # # f.savefig("scatter.png", dpi = 600, bbox_inches = "tight")

    # # print("\n R2 score: " + str(metrics.r2_score(y_te, y_pred)))
    # # print("\n ABS score: " + str(metrics.median_absolute_error(y_te, y_pred)))
    # #U,S = np.linalg.eig(np.matmul(X_tr,X_tr.T))
    # #print(np.diagonal(S))
    # # kf = KFold(n_splits = 2)
    # # model = RidgeCV(alphas=alphas,cv=5)
    # # RLM =  model.fit(X_tr, y_tr)
    # # print("\n ========== \n")
    # alpha = regr.alpha_
    # # y_pred = RLM.predict(X_te)
    # # print("R-squared training " + str(RLM.score(X_tr, y_tr)))
    # # print("alpha is " + str(alpha))
    # # MSE = np.sqrt(metrics.mean_squared_error(y_te,y_pred))
    # # print("MSE is " + str(MSE))
    # # scores = cross_val_score(RLM, X, status, cv = 5)
    # # print ("Cross-validated scores : ")
    # # print(scores)

    # # fig = plt.scatter(y_te, y_pred)
    # # plt.xlabel("True Values")
    # # plt.ylabel("Predictions")
    # #plt.savefig("pred.png", bbox_inches = "tight", dpi = 600)
    # # n_folds = 5
    # # Ridge_model = Ridge(random_state = 0, max_iter = 10000)
    # # tuned_params = [{'alpha': alphas}]
    # # clf = GridSearchCV(Ridge_model, tuned_params, cv=n_folds, refit=False)
    # # clf.fit(X, status)
    # # scores = clf.cv_results_['mean_test_score']
    # # print(scores)
    # # scores_std = clf.cv_results_['std_test_score']
    # # fig = plt.figure().set_size_inches(8, 6)
    # # plt.semilogx(alphas, scores)

    # # # plot error lines showing +/- std. errors of the scores
    # # std_error = scores_std / np.sqrt(n_folds)

    # # plt.semilogx(alphas, scores + std_error, 'b--')
    # # plt.semilogx(alphas, scores - std_error, 'b--')
    # # plt.ylabel('CV score +/- std error')
    # # plt.xlabel('alpha')
    # # plt.axhline(np.max(scores), linestyle='--', color='.5')
    # # plt.xlim([alphas[0], alphas[-1]])
    # # plt.savefig("alpha_scores.png",dpi = 600, bbox_inches = "tight")

    # #print(X.shape)
    # #print("Determinant of cluster data: ", np.linalg.det(X))
    # #sys.exit(0)
    # # print("R-squared for test : " + str(RLM.score(X_te, y_te)))
    # # print ("R-squared for test metrics : " + str(metrics.r2_score(y_te, y_pred)))
    # # print("\n ========== \n")
    # pvals = RS.coef_pval(regr, X_te, y_te, alpha)
    # return pvals
####################################################################
######################## META ANALYSIS TO REVIVE?
####################################################################
 #   ct += 1
        # print(pd_beta)
        # print("\n =============== \n")
        # print(pd_se)
        # pvals = metal(pd_beta, pd_se)

        #assoc = np.where(pvals < pvalue)[0]
        #print(assoc[:15])

        #CS[k] = len(np.where(assoc < 10))
        #SP[k] = len(assoc) - CS[k]
        #print(pd_pval)
        #finpv = np.zeros(n)
        #for i in range(n):
           # print(sp.stats.combine_pvalues(pd_pval.iloc[i,:],method='fisher')[1])
            #finpv[i] = sp.stats.combine_pvalues(pd_pval.iloc[i,:],method='mudholkar_george')[1]
        #finpv[isnan(finpv)] = 0

        #print(finpv[:11])
        #if sum(finpv) != 0:
        #   nnzidx = np.nonzero(finpv)[0]
        #   finidx = [i for i in nnzidx if finpv[i] < pvalue]
        #else:
        #   print("All pvals are ZERO\n")
        #   finidx = [-9]
        #finidx = np.array(finidx)
        #if finidx.all() != -9:
         #  sig3 = len(np.where(finidx.astype(int) < 10)[0])
         #  SP[k] = abs(len(finidx) - sig3)
        #   CS[k] = sig3
        #else:
         #  SP[k] = 0
         #  CS[k] = 0
        #intidx.sort()
        #print(len(intidx))
        #print(intidx[:10])

        #CS[k] = len(np.where(intidx < 10)[0])
        #SP[k] = len(intidx) - CS[k]
 ### for Meta analysis purposes
# def metal(beta, se):
    # m,n = beta.shape
    # pdwt = se.apply(np.reciprocal)
    # pdwt = pdwt.pow(2)
    # print(pdwt)
    # sumwt = list(pdwt.sum(axis=1))
    # print(sumwt)
    # newSE = list(np.sqrt(np.reciprocal(sumwt)))
    # print(newSE)

    # betawt = beta.multiply(pdwt)
    # print(betawt)
    # sum_betawt = list(betawt.sum(axis=1))
    # beta_q = [b/m for b,m in zip(sum_betawt, sumwt)]
    # print(beta_q[:10])
    # print(newSE[:10])
    # Zsc = [b/m for b,m in zip(beta_q,newSE)]
    # Zsc = np.array(Zsc)
    # pvals = scipy.stats.norm.sf(abs(Zsc))*2

    # return pvals
