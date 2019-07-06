import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random

import fastcluster
import scipy.cluster.hierarchy as sch

import ArmitChisq

def scipy_cluster(data, method):
    clusters = sch.linkage(data, method=method, metric='euclidean')
    return clusters

def fast_cluster(data, method):
    clusters = fastcluster.linkage(data, method=method, metric='euclidean', preserve_input=False)
    return clusters

def cluster(normR, U, D, pops, status, chisqst, dele):

    clustcount = np.zeros((len(dele),1))
    CS = np.zeros((len(dele),1))
    allidx = []
    Z = scipy_cluster(D,'ward')

    for k in range(0,len(dele)):
        numclust = pops+dele[k]

        den2 = sch.dendrogram(Z, p=30, truncate_mode='lastp', show_leaf_counts=False,
                leaf_rotation=90.,orientation='left', leaf_font_size=10)
        plt.title('Dendrogram for the clustering')
        plt.xlabel('Distance')
        # plt.ylabel('Distance')
        plt.savefig('Z_dendogram_'+str(k)+'.png')#,format='png', dpi=1000)

        plt.clf()
        plt.cla()
        plt.close()
        # dendrogram(Z,'Orientation','left','ColorThreshold','default');
        ClustMem = sch.fcluster(Z, numclust, depth=2e5)
        # ClustMem = cluster(Z,'Cutoff',numclust,'Depth',2e5);
        # print(ClustMem.shape)
        # print(np.unique(ClustMem))
        oldidx = []
        sti = []

        print('The number of clusters is : ', np.unique(ClustMem).shape[0]);
        clustcount[k] = np.unique(ClustMem).shape[0]
        # print(clustcount)
        # print(' ')

        for i in range(0, np.unique(ClustMem).shape[0]):
            print('cluster ID is: ',i);
            Ri = normR[ClustMem == i+1,:];
            sti = status[ClustMem == i+1];
            pvalc = np.zeros((normR.shape[1],1));

            if (len(sti) > 1):
                for j in range(0,normR.shape[1]):
                    pvalc[j] = ArmitChisq.chisq(Ri[:,j], sti.T, len(sti))

            pvidx = pvalc[np.where(pvalc > chisqst)[0]]

            # 'nan's coming from 'corrcoef' since the stddev is 0
            # print(pvalc[np.where(np.isnan(pvalc) == True)[0]])
            # these also occur in the Matlab data but is ignored (no warnings thrown)

            if len(oldidx) == 0:
                oldidx = pvidx;
            else:
                intidx = np.intersect1d(oldidx,pvidx);
                oldidx = intidx

        allidx = [allidx, -9, oldidx.T]
        CS[k] = len(oldidx);
        print(' ')

    return CS, clustcount, allidx

if __name__ == '__main__':
    print("HEY HEY HEY")
