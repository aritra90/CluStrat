import matplotlib
matplotlib.use('Agg')
import matplotlib as pyplot
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random, warnings, array, shutil, glob, argparse
warnings.filterwarnings(action="ignore")
import scipy.cluster.hierarchy as sch
import sklearn.metrics as metrics

from pysnptools.snpreader import Bed, Pheno
from pysnptools.kernelreader import SnpKernel
from pysnptools.standardizer import Unit
from pysnptools.util import intersect_apply

pyplot.rcParams['figure.dpi']=600
pyplot.rcParams.update({'font.size': 10})

warnings.filterwarnings("ignore")
sns.set(style='darkgrid')
# get_ipython().run_line_magic('matplotlib', 'inline')
# get_ipython().run_line_magic('load_ext', 'memory_profiler')
# get_ipython().run_line_magic('load_ext', 'line_profiler')

import getMH

def msg(name=None):
    return '''
         >> python3 CluStrat_wrapper.py --sim 1 --prop 10,20,70 --trait
        '''

def parse_arguments():
    parser = argparse.ArgumentParser(usage=msg())

    parser.add_argument("-f", "--fn", dest='data_filename', action='store', help="Put the full path where the dataset is stored.",
                        metavar="FN")

    parser.add_argument("-ct", "--maxclust", dest='maxclust', action='store', help="Enter clustering threshold (max inter-cluster distance allowed for scipy.cluster.hierarchy.fcluster).",
                        metavar="NC")

    parser.add_argument("-k", "--mhdistk", dest='k', action='store', help="Enter value of k for Mahalanobis distance approximation.",
                        metavar="MH")

    args = parser.parse_args()

    return args

def get_data(bed_fn):

    ignore_in = None
    #store snp data from Bedfile
    snp_on_disk = Bed(bed_fn, count_A1=False)
    kernel_in = SnpKernel(Bed(bed_fn, count_A1=False), Unit())

    #intersect and sort IDs from bed, pheno returning new data with same IDs
    ignore_out, bed_out = intersect_apply([ignore_in, snp_on_disk])

    sample_id = bed_out.iid
    num_samples = bed_out.iid_count
    # print("Number of samples: ", bed_out.iid_count)
    snp_id = bed_out.sid
    snp_chrom = bed_out.pos
    snp_chrom[0][1] = 0; snp_chrom[0][2] = 0; snp_ids_and_chrom = pd.DataFrame(np.hstack([snp_id.reshape(snp_id.shape[0], 1),snp_chrom]))
    snp_ids_and_chrom.columns = ['SNP','Chr','ChrPos', 'ChrPos2']
    num_snps = bed_out.sid_count
    normZ = bed_out.read().standardize().val

    return normZ, snp_ids_and_chrom, num_samples, num_snps

def scipy_cluster(data, method):
    clusters = sch.linkage(data, method=method, metric='euclidean');
    return clusters

def cluster(D, t, depth = 150):
    Z = scipy_cluster(D,'ward');
    den = sch.dendrogram(Z, leaf_rotation=90.,orientation='left', leaf_font_size=2)
    plt.savefig('Z_dendogram.png',format='png', dpi=600)

    # maxclust :
    #     Finds a minimum threshold r so that the cophenetic distance between any two original observations in the same flat cluster is no more than r and no more than t flat clusters are formed.

    ClustMem = sch.fcluster(Z, t, criterion = "maxclust")
    return ClustMem

if __name__ == '__main__':
    args = parse_arguments()
    data_fn = args.data_filename
    bed_fn = data_fn + ".bed"
    fam_fn = data_fn + ".fam"
    pheno_fn = data_fn + ".phen"
    k = int(args.k)
    maxclust = int(args.maxclust)
    print("mahalanobis distance k:", k)
    print("clutering threshold:", maxclust)


    # read in fam ids and assign to each corresponding cluster
    FIDs = []
    IIDs = []
    status = []
    with open(fam_fn, 'r') as f:
        for line in f:
            elems = line.split()
            # print(elems)
            FIDs.append(elems[0].strip())
            IIDs.append(elems[1].strip())
            # status.append(elems[5].strip())
    # status = np.array(status).astype(int)
    FIDs = np.asarray(FIDs)
    IIDs = np.asarray(IIDs)

    with open(pheno_fn, 'r') as f:
        for line in f:
            elems = line.split()
            # print(elems[2].strip() == '-9')
            if elems[2].strip() == '-9':
                status.append(1)
            else:
                status.append(elems[2].strip())
    status = np.array(status).astype(int)

    # read in the data
    print("reading in the data...")
    normZ, snp_ids_and_chrom, num_samples, num_snps = get_data(bed_fn)

    # # row indices of what samples are cases
    # case_idx = np.where(status != 1)[0]

    # mahalanobis distance
    print("computing mahalanobis distance...")
    D = getMH.MH(normZ, k = k)
    row_sums = D.sum(axis=1)
    D = D / row_sums[:, np.newaxis]

    # clustering
    print("clustering...")
    clustmem = cluster(D, t = maxclust)
    # print(clustmem)
    num_clusters = np.unique(clustmem).shape[0]
    print('the number of clusters is:', num_clusters);
    np.savetxt('original_clustmem', clustmem, delimiter=',', fmt='%f')

    # print("DEBUGGING RANDOM CLUSTERING")
    # clustmem = np.random.randint(maxclust, size = normZ.shape[0])
    # # print(clustmem)
    # num_clusters = np.unique(clustmem).shape[0]
    # print('the number of clusters is:', num_clusters);
    # np.savetxt('random_clustmem', clustmem, delimiter=',', fmt='%f')

    print("generating output...")
    for i in range(num_clusters):
        print("cluster",i+1)

        # df = pd.DataFrame()
        # df['FID'] = FIDs[np.where(clustmem == i+1)[0]]
        # df['IID'] = IIDs[np.where(clustmem == i+1)[0]]
        # df['status'] = status[np.where(clustmem == i+1)[0]]
        #
        # iids_df = df.drop(columns=['status'])
        # pheno_df = df.sort_values('IID')
        #
        # iids_df.to_csv("cluster_iids_" + str(i+1), sep='\t', index = False, header = None)
        # pheno_df.to_csv("pheno_cluster_iids_" + str(i+1), sep='\t', index = False)



        # for binary trait and balancing
        if np.unique(status[np.where(clustmem == i+1)[0]]).shape[0] == 2 and np.unique(status[np.where(clustmem == i+1)[0]], return_counts = True)[1][1] > 5:
            df = pd.DataFrame()
            df['FID'] = FIDs[np.where(clustmem == i+1)[0]]
            df['IID'] = IIDs[np.where(clustmem == i+1)[0]]
            df['status'] = status[np.where(clustmem == i+1)[0]]

            # remove controls from dataframe proportional to the number of cases here
            num_cases = np.unique(df['status'], return_counts = True)[1][1]
            num_ctrls = np.unique(df['status'], return_counts = True)[1][0]
            print("number of cases:", num_cases)
            print("number of controls:", num_ctrls)


            print("no balancing...")
            cases = df.loc[df['status'] == 2]
            controls = df.loc[df['status'] == 1]

            # print("balancing...")
            # if num_ctrls > num_cases:
            #     cases = df.loc[df['status'] == 2]
            #     controls = df.loc[df['status'] == 1].head(num_cases)
            # else:
            #     # in case controls out-represent cases (which is never but I need this here)
            #     cases = df.loc[df['status'] == 2].head(num_ctrls)
            #     controls = df.loc[df['status'] == 1]
            # # print("number of cases & controls after balancing:", cases.shape[0], controls.shape[0])

            # print("number of controls:", controls.shape[0])
            print(" ")
            frames = [cases, controls]
            df = pd.concat(frames)
            iids_df = df.drop(columns=['status'])
            # pheno_df = df.drop(columns=['FID'])
            pheno_df = df.sort_values('IID')

            iids_df.to_csv("cluster_iids_" + str(i+1), sep='\t', index = False, header = None)
            pheno_df.to_csv("pheno_cluster_iids_" + str(i+1), sep='\t', index = False)
            # print(df)
        else:
            print("NOT ENOUGH CASES IN CLUSTER", i+1)
            print(np.unique(status[np.where(clustmem == i+1)[0]], return_counts = True))
            print(" ")








# end of code
