import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random, warnings
warnings.filterwarnings(action="ignore")
import fastcluster
import normalize
from sklearn.model_selection import KFold 
import RegresStat as RS
import scipy.cluster.hierarchy as sch
from sklearn import linear_model, svm
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.linear_model import RidgeCV, Ridge, LassoCV
from sklearn.model_selection import train_test_split, GridSearchCV
import sklearn.metrics as metrics
from sklearn.compose import TransformedTargetRegressor
from sklearn.preprocessing import QuantileTransformer, quantile_transform
from scipy.stats import t
from sklearn.kernel_ridge import KernelRidge
import array 

# import regressors as rgrs
from regressors import stats

import ArmitChisq

def scipy_cluster(data, method):
    clusters = sch.linkage(data, method=method, metric='euclidean')
    return clusters

def fast_cluster(data, method):
    clusters = fastcluster.linkage(data, method=method, metric='euclidean', preserve_input=False)
    return clusters

def ridge_pvals(X, Y): 
####################################################################
    lambdas = np.logspace(-2.5,2.5,25)
    regr = RidgeCV(alphas=lambdas, cv = 5)
    regr.fit(X,Y)
    y_pred = regr.predict(X) 
    print("\n =============== \n")
 
    #print("\n Lambda is : "+str(regr.alpha_))
    print("\n R2 score : "+str(regr.score(X,Y)))
    print("\n MSE : " + str(metrics.mean_squared_error(y_pred, Y)))
    print("\n =============== \n")
    #alpha = regr.alpha_ 
    
    pvals = RS.coef_pval(regr, X, Y, regr.alpha_)	
    return pvals
	


def cluster(normR, U, D, pops, status, pvalue, dele):

    clustcount = np.zeros((len(dele),1))
    CS = np.zeros((len(dele),1))
    SP = np.zeros((len(dele),1))
    allidx = []
    Z = scipy_cluster(D,'ward')

    for k in range(0,len(dele)):
        numclust = pops+dele[k]

        
        den2 = sch.dendrogram(Z, leaf_rotation=90.,orientation='left', leaf_font_size=2)
        #plt.title('Dendrogram for PSD model \{$\alpha=0.1$\}')
        #plt.xlabel('Distance')
        # plt.ylabel('Distance')
        plt.savefig('Z_dendogram_'+str(k)+'.png')#,format='png', dpi=1000)

        plt.clf()
        plt.cla()
        plt.close()
        # dendrogram(Z,'Orientation','left','ColorThreshold','default');
        ClustMem = sch.fcluster(Z, numclust, depth=150)
        # ClustMem = cluster(Z,'Cutoff',numclust,'Depth',2e5);
        #print(ClustMem.shape)
        # print(np.unique(ClustMem))
        oldidx = []
        intidx = []
        sti = []
        combset = []
        print('The number of clusters is : ', np.unique(ClustMem).shape[0]);
        clustcount[k] = np.unique(ClustMem).shape[0]
        # print(clustcount)
        # print(' ')
        ct = 0
        # for each cluster...
        for i in range(0, np.unique(ClustMem).shape[0]):
            print('cluster ID is: ',i)
            Ri = normR[ClustMem == i+1,:]
            sti = status[ClustMem == i+1]
            pvalc = np.zeros((normR.shape[1],1))

            ####################################################################
            # get wieghts of all SNPs using Ridge regression (output are SNP ids?)
            ridge_pval = ridge_pvals(Ri,sti)
            print("\n Min Ridge p-value \n")
            print(min(ridge_pval))
            if sum(ridge_pval) != 0:
                nnzidx = np.nonzero(ridge_pval)
                #pvidx = [i for i in ridge_pvals if i < pvalue]
                pvidx = np.where( ridge_pval[nnzidx] < pvalue)[0]
                #print(pvidx)
            else:
                print("\n All Ridge pvals are ZERO")
            ####################################################################

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
            print(pvidx)
            nnzidx = nnzidx[0].astype(int)
            print("Cluster " + str(i) + " -- number of significant SNPs : " 
			                 + str(len(pvidx)))
            #print(nnzidx) 
            #print(ridge_pval[nnzidx[pvidx]])
            
            if len(oldidx) == 0:
                oldidx = pvidx
                ct = ct + 1
            elif ct == 1:
                #intidx = np.intersect1d(oldidx,pvidx);
                #oldidx = intidx
                intidx = np.intersect1d(oldidx,pvidx)
                #print("intidx : " + str(len(intidx)))
                combset = np.unique(np.append(oldidx, pvidx))
                #print("combset : " + str(len(combset)))
                ct = ct + 1
				#if len(np.intersect1d(oldidx, pvidx)) > len(intidx):
            else:   
                comb_intidx = np.intersect1d(combset,pvidx)
                #print("comb_intidx : " + str(len(comb_intidx)))
                if len(intidx) < len(comb_intidx):
                   intidx = comb_intidx
                combset = np.unique(np.append(combset,pvidx))
                print("combset : " + str(len(combset)))
                ct = ct + 1 				
				#oldidx = np.append(oldidx,pvidx)
        print(ct)				
        #allidx = [allidx, -9, oldidx.T]
        #oldidx = np.unique(oldidx)
        print(len(intidx))
        print(' ')
        sig3 = [i for i in intidx if i < 10]
        SP[k] = len(intidx) - len(sig3)
        CS[k] = len(sig3);
        #print(sig3)
        print("number of causal SNPs in CluStrat = " + str(len(sig3)))
        print("Number of spurious associations " + str(len(intidx) - len(sig3)))

    return CS, clustcount, SP


# def ridge_reg_pvals(X, status):
# ####################################################################
	# #X1 = np.hstack((np.ones((n, 1)), np.matrix(X)))
    # #X_tr, X_te, y_tr, y_te = train_test_split(X, status, test_size=0.3,
    # #	                                    random_state=0)
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
if __name__ == '__main__':
    print("HEY HEY HEY")
