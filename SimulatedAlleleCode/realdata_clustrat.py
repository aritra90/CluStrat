from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from plinkio import plinkfile
import parseplink as pp
from scipy.sparse.linalg import svds
import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random
from datetime import datetime
import statsmodels.api as sm
from scipy.stats.distributions import chi2
from scipy.spatial.distance import pdist, squareform
from sklearn.model_selection import train_test_split
import normalize, ArmitChisq, EigStrat, CluStrat, getMH

def read_handlers(file_handle):
    print(file_handle)
    plink_file = plinkfile.open(file_handle)
    X, pheno = pp.read_n_parse(plink_file)
    return X, pheno

if __name__ == '__main__':

    file_handle = str(sys.argv[1])
    sketch_flag = sys.argv[2]
    plot_flag = sys.argv[3]
    begin_time = time.time()
    print('##################### Loading in data...')
    load_time = time.time()
    X, pheno = read_handlers(file_handle)
    print("Loaded genotype matrix of dimension ", X.shape)
    print('Loading time (secs): ', (time.time()-load_time))
    print(' ')
    
    #print(pheno.shape)
    # we need SNPs x inidvs here
    print('##################### Transposing data...')
    #transp_time = time.time()
    #X = X.T
    #pheno = np.reshape(pheno, (pheno.shape[0], 1))
    #print('Transposing time (secs): ', (time.time()-transp_time))
    #print(' ')
    m = X.shape[0]
    n = X.shape[1]
    pvalue = 1e-3

    print('##################### Normalizing data...')
    norm_time = time.time()
    normX,_ = normalize.norm(X,0)
    print(normX.shape)
    print('Normalizing time (secs): ', (time.time()-norm_time))
    print(' ')
    # traits, status = traitSim.simulate(normX,S,v,m,n,d)
    Y = pheno
    
    #print(np.unique(Y, return_counts=True))
    
    if (plot_flag == 1):
    	svd_time = time.time()
    	temp = np.matmul(normX,normX.T)
    	U, Sig, _ = svds(temp, k =10)
    	print('SVD time (mins): ', (time.time()-svd_time)/60.0)
    	print(' ')
    	print(np.sqrt(Sig))
    	# PCA PLOTS
    	# grab indices for cases and controls
    	idxcase = np.where(Y == 1)[0]
    	idxcontr = np.where(Y == 0)[0]

    	plt.style.use('seaborn-darkgrid')
    	fig = plt.figure()
    	ax = fig.add_subplot(111)
    	scale = 200.0 * np.random.rand(750)
    	print('# Cases: '+str(len(idxcase)))
    	print('# Controls: '+str(len(idxcontr)))
    	ax.scatter(U[idxcase, 0], U[idxcase, 1], marker='o', color='red', s=10, label='case')
    	ax.scatter(U[idxcontr, 0], U[idxcontr, 1], marker='*', color='blue', s=10, label='control')
    	ax.legend()
    	ax.set_xlabel('PC 1')
    	ax.set_ylabel('PC 2')
    	plt.savefig('clustrat_top2PCs.png', bbox_inches='tight', dpi=600)


    print('##################### Getting distance matrix...')
    dist_time = time.time()
    # 1000 choose 2 pairs (calculating distance between each pair)
    #D = squareform(pdist(normX))
    D = getMH.MH(normX)
    print('Calculating distance matrix time (mins): ', (time.time()-dist_time)/60.0)
    print(' ')
    print(D)
    dele = [8, 12]
    d = 2
    # normR = normX.T

    print('##################### Running CluStrat...')
    clu_time = time.time()
    SP, CS, clustcount = CluStrat.cluster(X, D, d, Y, pvalue, dele, sketch_flag)
    print('CluStrat time (mins): ', (time.time()-clu_time)/60.0)
    #print(' ')
    #print(sigsnp_indices)
    # np.savetxt('CluStrat_pvals.txt', pvals, delimiter =', ')
    CS = list(CS)
    SP = list(SP)

    # basically take the max "top PCs" from the clusterings and that is how
    # many more candidate SNPs we have (also record for which cluster gave the max)

    #maxidx = CS.index(max(CS))
    SPmax = max(SP)

    SP3 = SPmax
    #CS3 = CS[maxidx]
    print('Total time (mins): ', (time.time() - begin_time)/60.0)
