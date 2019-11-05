from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random
from datetime import datetime

def manhattan_plot(data, prefix, group=None):
    plt.style.use('seaborn-darkgrid')
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plotted_data = (-1)*np.log10(data.iloc[:,-1])

    # Mark sigSNPs based on some criteria
    #color_arr = list(np.where(plotted_data>5,'b','r'))
    print(data.shape)
    # pvalue = 1e-3
    # thresh = (-1)*math.log10(25.0/data.shape[0])
    thresh = (-1)*math.log10(1e-3)
    print(thresh)

    sigSNPs_idx = np.where( plotted_data>thresh )
    nonsigSNPs_idx = np.where(plotted_data <= thresh)

    print(sigSNPs_idx)
    print(nonsigSNPs_idx)
    print(' ')

    #print(sigSNPs_idx)
    #print(group)
    #print(group[sigSNPs_idx])

    ax.scatter(sigSNPs_idx[0] , plotted_data[sigSNPs_idx[0]], marker='d', color='b', s=8)
    # cmap = 'magma' is cool, 'plasma' aight too
    if group is None:
        ax.scatter(nonsigSNPs_idx[0] , plotted_data[nonsigSNPs_idx[0]], marker='o', s=8)
    else:
        ax.scatter(nonsigSNPs_idx[0] , plotted_data[nonsigSNPs_idx[0]], marker='o', c=group[nonsigSNPs_idx[0]], cmap='inferno', s=8)

    #ax.scatter(np.arange(0,data.shape[0]) , plotted_data, marker='o', c=group, s=8)
    ax.plot(np.arange(0,data.shape[0]) , thresh*np.ones((data.shape[0],1)), '--', color='k')
    # ax.legend()
    ax.set_xlabel('SNPs')
    ax.set_ylabel('-log(p-value)')
    plt.savefig('mnhtn_'+str(prefix)+'.png', bbox_inches='tight', dpi=600)

if __name__ == '__main__':
    # Grab the eignenvector file (with first column and row removed)
    # COMMAND = "sed '1d' output_singularVectors.txt | awk '{$1=\"\"}1' > singVectors.txt"
    # subprocess.call(COMMAND, shell=True)
    # X = np.loadtxt("singVectors.txt")

    COMMAND = "awk '{print $1}' PKHCGSRRS_pruned.bim > chromosomeLabel.txt"
    subprocess.call(COMMAND, shell=True)
    group = np.loadtxt("chromosomeLabel.txt")

    #data = pd.read_csv('PKHCGSRRS_pruned_assoc_outfile.ps',header=None,delimiter='\t')
    #print(data.head())
    #print(data.shape)
    #print(' ')
    #manhattan_plot(data, 'emmax', group)

    #data = pd.read_csv('output/PKHCGSRRS_pruned_LMM.assoc.txt',delimiter='\t')
    #print(data.head())
    #print(data.shape)
    #print(' ')
    #manhattan_plot(data, 'gemma', group)

    #data = pd.read_csv('ArmitChisq_pvals.txt',header=None,delimiter='\t')
    #print(data.head())
    #print(data.shape)
    #print(' ')
    #manhattan_plot(data, 'armitchisq', group)

    #data = pd.read_csv('EigStrat_pvals.txt',header=None,delimiter='\t')
    #print(data.head())
    #print(data.shape)
    #print(' ')
    #manhattan_plot(data, 'eigstrat', group)

    # file_list = ["CluStrat_pvals_pt1_0.txt", "CluStrat_pvals_final.txt"]

    file_list = ["CluStrat_pvals_pt1_0.txt", "CluStrat_pvals_final_0.txt"]

    other_list = ["CluStrat_pvals_final_0.txt"]

    for elem in file_list:
        # e.g. CluStrat_24_pvals.txt
        data = pd.read_csv(elem,header=None,delimiter=',')
        data = data.astype(float)
        # data = data.T
        data = data.replace(np.nan, 1.0)
        data = data.replace(0.0, 1.0)
        # print(data.info())
        print(data.head())
        # print(' ')
        if elem in other_list:
            manhattan_plot(data, elem.split('.')[0])
        else:
            manhattan_plot(data, elem.split('.')[0], group)
