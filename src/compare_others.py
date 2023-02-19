import random, time, itertools, subprocess, math, os, normalize, csv, sys, pdb, argparse
import EigStrat, ArmitChisq

from scipy.spatial.distance import pdist
from scipy.stats.distributions import chi2
from sklearn.cluster import KMeans
from subprocess import Popen, PIPE
from datetime import datetime
from operator import add

import numpy as np
import pandas as pd

from pysnptools.snpreader import Bed, Pheno
from pysnptools.kernelreader import SnpKernel
from pysnptools.standardizer import Unit
from pysnptools.util import intersect_apply

####################################################################
###################   DEFINE VARIABLES   ###########################
####################################################################
PATH_TO_GEMMA = '/depot/pdrineas/data/gwas_software/gemma/gemma-0.98.1'
PATH_TO_EMMAX_KIN = '/depot/pdrineas/data/gwas_software/emmax/emmax-kin'
PATH_TO_EMMAX = '/depot/pdrineas/data/gwas_software/emmax/emmax'
PATH_TO_PLINK = '/depot/pdrineas/data/gwas_software/plink2/plink'
# PRE_PATH = subprocess.call(['pwd'])

def msg(name=None):
    return '''
         >> python3 compare_others.py --fn /depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/MODEL_simdata/simdata
        '''

def parse_arguments():
    parser = argparse.ArgumentParser(usage=msg())

    parser.add_argument("-f", "--fn", dest='data_filename', action='store', help="Put the full path where the dataset is stored.",
                        metavar="FN")

    parser.add_argument("-p", "--pval", dest='pval', action='store', help="Enter p-value threshold.",
                        metavar="PVAL")

    parser.add_argument("-m", "--model", dest='model', action='store', help="Enter model flag.",
                        metavar="MODEL")

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

if __name__ == '__main__':

    ##############################################################
    ####################       READ DATA       ###################
    ##############################################################
    args = parse_arguments()
    data_fn = args.data_filename
    bed_fn = data_fn + ".bed"
    fam_fn = data_fn + ".fam"
    pheno_fn = data_fn + ".phen"
    pval_thresh = float(args.pval)
    model = str(args.model)
    causal_df = pd.read_csv("/depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/"+model+"_simdata/causal_variants.txt", header = None)
    thesnps = np.reshape(causal_df.values, causal_df.values.shape[0])

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
    # print("reading in the data...")
    normZ, snp_ids_and_chrom, num_samples, num_snps = get_data(bed_fn)
    n, m = normZ.shape


    ###########################################################################
    ####################      RUN ARMITAGE TREND CHISQ      ###################
    ###########################################################################
    print("running armitage-trend...")
    # st2 = time.time()
    # % Get p-values based on correlation between each individual and the status
    # % (for every SNP)
    # % essentially getting a measurement of how well each SNP is correlated to
    # % the status
    pvstat = np.zeros((m, 1))
    for i in range(0, m):
        pvstat[i] = ArmitChisq.chisq(normZ[:, i], status, n)
    # Find a value that exceeds (1-0.0025)% of the samples from chi-sq with 1
    # degree of freedom
    # Only going to observe values greater than 'chisqst' 0.0025% of the time
    chisqst = chi2.ppf(1-pval_thresh, df=1)
    # Grab pvals that are larger than 'chisqst'
    sigset = np.where(pvstat > chisqst)[0]
    sigsetinfo = snp_ids_and_chrom.iloc[sigset]
    causal_hits = np.isin(sigsetinfo.SNP,thesnps)

    # print(thesnps)
    #
    # print(sigsetinfo.SNP)
    #
    # print(list(causal_hits))

    causal = sum(bool(x) for x in causal_hits)
    spur = len(causal_hits) - causal

    # causal
    print("causal:", causal)
    # spurious
    print("spurious:", spur)
    print(" ")

    with open('causal_armitage', 'a+') as f:
        f.write(str(causal)+'\n')

    with open('spurious_armitage', 'a+') as f:
        f.write(str(spur)+'\n')



    ###########################################################################
    ####################           RUN EIGENSTRAT           ###################
    ###########################################################################
    print("running eigenstrat...")
    ##PCA step on indivs by indivs matrix
    # st3 = time.time()
    temp = np.matmul(normZ, normZ.T)
    U, Sig, _ = np.linalg.svd(temp)
    # project out the top 'K' PCs (in this case 10)
    K = 10
    adjR, adjstat = EigStrat.project(normZ, U[:, 0:K], status)
    # Get p-values based on correlation between each individual and the status
    # (for every SNP) *now on the adjusted data
    adjpvstat = np.zeros((m, 1))
    for i in range(0, m):
        adjpvstat[i] = ArmitChisq.chisq(adjR[:, i], adjstat.T, n-K-1)
    # Grab pvals that are larger than 'chisqst'
    adjsigset = np.where(adjpvstat > chisqst)[0]
    sigsetinfo = snp_ids_and_chrom.iloc[adjsigset]
    causal_hits = np.isin(sigsetinfo.SNP,thesnps)

    causal = sum(bool(x) for x in causal_hits)
    spur = len(causal_hits) - causal

    # causal
    print("causal:", causal)
    # spurious
    print("spurious:", spur)
    print(" ")

    with open('causal_eigenstrat', 'a+') as f:
        f.write(str(causal)+'\n')

    with open('spurious_eigenstrat', 'a+') as f:
        f.write(str(spur)+'\n')




    # # NOT GIVING ANY CAUSAL OR SPURIOUS SO TOSSING OUT FOR NOW CAUSE TAKES TOO LONG
    #
    # ######################################################################
    # ####################           RUN GEMMA           ###################
    # ######################################################################
    # print("running gemma...")
    # # st4 = time.time()
    # GRM_PATH = str(os.path.basename(data_fn))+'_GRM'
    # ASSOC_PATH = str(os.path.basename(data_fn))+'_ASSOC'
    # subprocess.call([PATH_TO_GEMMA, '-bfile', data_fn,
    #                  '-p', pheno_fn, '-gk', '1', '-o', GRM_PATH],stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # GRM_PATH = 'output/'+str(GRM_PATH)+'.cXX.txt'
    # subprocess.call([PATH_TO_GEMMA, '-bfile', data_fn, '-p', pheno_fn, '-k',
    #                  GRM_PATH, '-lmm', '2', '-o', ASSOC_PATH],stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # ASSOC_PATH = 'output/'+str(ASSOC_PATH)+'.assoc.txt'
    # gemmadf = pd.read_csv(ASSOC_PATH, sep = '\t')
    # gemmadf = gemmadf.sort_values('p_lrt')
    # # print(gemmadf)
    # hits = gemmadf.loc[gemmadf['p_lrt'] < pval_thresh]
    # causal_hits = hits.loc[hits['rs'].isin(thesnps)]
    # # causal
    # print("causal:", causal_hits.shape[0])
    # # spurious
    # print("spurious:", hits.shape[0] - causal_hits.shape[0])
    # print(" ")
    #
    #
    #
    # ######################################################################
    # ####################           RUN EMMAX           ###################
    # ######################################################################
    # print("running emmax...")
    # # st5 = time.time()
    # subprocess.call([PATH_TO_PLINK, '--bfile', data_fn, '--recode', '12',
    #                  '--output-missing-genotype', '0', '--transpose', '--out', data_fn],stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # subprocess.call([PATH_TO_EMMAX_KIN, '-v', '-d', '10', data_fn],stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # KIN_PATH = str(data_fn)+'.BN.kinf'
    # EMMAXOUT_PATH = 'output/' + \
    #     str(os.path.basename(data_fn))+'_EMMAX'
    # subprocess.call([PATH_TO_EMMAX, '-v', '-d', '10', '-t', data_fn,
    #                  '-p', pheno_fn, '-k', KIN_PATH, '-o', EMMAXOUT_PATH],stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # # output file format: SNPID SE P
    # EMMAXOUT_PATH = EMMAXOUT_PATH + '.ps'
    # emmaxdf = pd.read_csv(EMMAXOUT_PATH, header = None, names = ['SNP', 'SE', 'P'], sep = '\t')
    # emmaxdf = emmaxdf.sort_values('P')
    # # print(emmaxdf)
    # hits = emmaxdf.loc[emmaxdf['P'] < pval_thresh]
    # causal_hits = hits.loc[hits['SNP'].isin(thesnps)]
    # # causal
    # print("causal:", causal_hits.shape[0])
    # # spurious
    # print("spurious:", hits.shape[0] - causal_hits.shape[0])
    # print(" ")


















    # end of file
