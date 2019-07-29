from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random
from operator import add
from datetime import datetime

from sklearn.cluster import KMeans
from scipy.stats.distributions import chi2
from scipy.spatial.distance import pdist
import normalize, traitSim, ArmitChisq, EigStrat, getMH, CluStrat
from plinkio.plinkfile import WritablePlinkFile, Sample, Locus

def matrix_sketch(X, k):
    # # Calculate sketch parameters based on design matrix
    #
    # # Gaussian normal test matrices
    # mu, sigma = 0, 1
    # # This is the sketch we are concerned about because we are reducing the larger dimension
    # Omega = np.random.normal(mu, sigma, (X.shape[0],k))
    # sketch_X = np.matmul(X,Omega)
    #
    # # Psi = np.random.normal(mu, sigma, (l, X.shape[1]))
    # # sketch_X = np.matmul(Psi, X)

    # ############################################################################
    # # Set sketch parameters
    # # sufficiently large constant (c2 > 3330*15e2)
    # C = 4000*15e2
    # # 0 < epsilon < 1/3
    # epsilon = np.random.uniform(low = 2.22045e-16, high = .3333333)
    # # number of artificial features (do we guarantee that r < m for X in nxm?)
    # # m > k so maybe this is the only constraint needed
    # r = int(C*(k/epsilon**2))
    # # Fill the random matrix
    # Omega = np.random.choice([1/math.sqrt(r),-1/math.sqrt(r)],p=[0.5,0.5], size=(X.shape[0],r))
    # # Compute the sketch
    # sketch_X = np.matmul(X, Omega)

    ############################################################################
    # Set sketch parameters
    r = k
    # Fill the random matrix
    mu = 0
    sigma = 1.0/r
    Omega = np.random.normal(mu, sigma, (X.shape[1],r))
    # Compute the sketch
    sketch_X = np.matmul(X, Omega)

    return sketch_X.T

def convert2plink(X, status, model_flag):

    #create fake sample IDs
    sampleIDs = [str(model_flag)+"000"+str(s) for s in range(1,X.shape[1]+1)]
    #create fake rsIDs
    rsIDs = ["rs000"+str(s) for s in range(1,X.shape[0]+1)]
    #create fake positions increasing randomly between 200 and 20k centimorgans
    snpPositions = [float(s) for s in range(70000,70000+X.shape[0])]
    for k in range(1,X.shape[0]):
        snpPositions[k] = snpPositions[k-1] + random.randint(200,1000)

    # print(snpPositions)
    # print(' ')
    # print(sampleIDs)
    # print(" ")
    # print(rsIDs)

    list_of_Samples = []
    # Create samples
    for i in range(0,X.shape[1]):
        list_of_Samples.append(Sample(sampleIDs[i], sampleIDs[i], 0, 0, -9, status[i]))

    # print(list_of_Samples)

    list_of_Loci = []
    # Create loci
    for j in range(0,X.shape[0]):
        rng = random.randint(0,1)
        if rng == 0:
            alleles = ['A','T']
        else:
            alleles = ['C','G']

        # Choosing to encode 0 as 'A'/'C' and 2 as 'T'/'G'
        # Get the allele frequencies
        allele_counts = np.unique(X[j,:], return_counts = True)
        # print(allele_counts)
        # if 0,1,2 all occur in the SNP...
        if allele_counts[0].shape[0] == 3:
            if allele_counts[1][0] < allele_counts[1][2]:
                list_of_Loci.append(Locus(0, rsIDs[j], 0, snpPositions[j], alleles[0], alleles[1]))
            else:
                list_of_Loci.append(Locus(0, rsIDs[j], 0, snpPositions[j], alleles[1], alleles[0]))
        else:
            if 0 in allele_counts[0] and 1 in allele_counts[0]:
                # print("no occurrences of '2'...")
                list_of_Loci.append(Locus(0, rsIDs[j], 0, snpPositions[j], alleles[1], alleles[0]))
            elif 0 in allele_counts[0] and 2 in allele_counts[0]:
                # print("no occrrences of '1'...")
                if allele_counts[1][0] < allele_counts[1][1]:
                    list_of_Loci.append(Locus(0, rsIDs[j], 0, snpPositions[j], alleles[0], alleles[1]))
                else:
                    list_of_Loci.append(Locus(0, rsIDs[j], 0, snpPositions[j], alleles[1], alleles[0]))
            else:
                # print("no occurrences of '0'...")
                list_of_Loci.append(Locus(0, rsIDs[j], 0, snpPositions[j], alleles[0], alleles[1]))
            # sys.exit(0)

    file_prefix = 'sim_plinkfiles/cluster_'+model_flag+'_'+str(random.randint(0, 9999999))
    print('File prefix: '+file_prefix)

    # Create PLINK files corresponding to the data
    plink_obj = WritablePlinkFile(file_prefix,list_of_Samples)
    # plink_obj = WritablePlinkFile(file_prefix,[Sample('0', '0', '0', '0', 0, 0)])


    for i in range(0,X.shape[0]):
        # print(X[i,:])
        # print(X[i,:].shape)
        plink_obj.write_row(list_of_Loci[i], X[i,:])
    # plink_obj.loci = list_of_Loci

    # print("YAY")
    # sys.exit(0)

    return plink_obj

if len(sys.argv) == 1:
    model_flag = "BN"
else:
    model_flag = str(sys.argv[1])

print("Model type: ", model_flag)

plot_flag = 1
sketch_flag = 1

start = time.time()

if model_flag == "BN":
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%%%%%% Load HapMap Info
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HM_inf = pd.read_csv('CEUASWMEX_fst_frq.txt',sep=' ')
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %get allele freq and Fst for each SNP

    # % allele freq: frequency (proportion) of that SNP in the data
    # % Fst: % of contribution to the total genetic variation that each SNP has
    frq = HM_inf['FRQ'].values #cell2mat(HM_inf(:,4));
    Fst = HM_inf['FST'].values #cell2mat(HM_inf(:,3));

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # % REMEMBER: 'n' here is INDIVIDUALS not SNPs
    # % HapMap3 Data ~1000 individuals and ~1000000 SNPs
    m = int(1e4) #number of SNPs
    n = int(1e3) #number of individuals
    pvalue = 0.0025
    d = 3 #number of populations

    # % each row of Gamma will be populated with 3 i.i.d. draws from BN model
    # % thanks to the nature of the HapMap data (sampled from 3 populations)
    #
    # % allele freq of population: allele freq of each SNP described by that
    # % population
    G = np.zeros((m,d)); #allele freq for each population

    # % columns of S will be populated with indicator vectors s.t. each
    # % individual assigned to one of the 3 subpopulations i.e. admixture
    S = np.zeros((d,n)); #individual population admixture

    # %X = zeros(m,n); %genotype matrix
    # % random seeding here...

    random.seed(datetime.now())

    # %populate the allele freq matrix from BN with (p,F) from HapMap
    # % for each SNP...
    for i in range(0,m):
        # each row will generate 'd' variates drawing from this distribution
        G[i,:] = np.random.beta(frq[i]*(1-Fst[i])/Fst[i], (1-frq[i])*(1-Fst[i])/Fst[i], size=d)

    # print('Mean of allele freq matrix: ', np.mean(G, axis=0))

    # %populate the population admixture matrix
    # %this part is tailor made for 3 populations as described in
    # %Song et al. Nat. Genet. (2015)
    # % Treating the probabilities as ranges
    # % 1: <60/210, 2: bet. 60/210 and 120/210, 3: >=120/210
    popidx = np.zeros((n,1));
    for i in range(0,n):
        p = random.uniform(0, 1);
        if p < (60.0/210):
            pick = 1;
            popidx[i]= 1;
        elif p < (2*(60.0/210)):
            pick = 2;
            popidx[i] = 2;
        else:
            pick = 3;
            popidx[i] = 3;

        S[pick-1,i] = 1;
        # Leaving all other values at zero (indiv only assigned to one subpop)

elif model_flag == "PSD":
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%%%%%% Load HapMap Info
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HM_inf = pd.read_csv('CEUASWMEX_fst_frq.txt',sep=' ')
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %get allele freq and Fst for each SNP

    # % allele freq: frequency (proportion) of that SNP in the data
    # % Fst: % of contribution to the total genetic variation that each SNP has
    frq = HM_inf['FRQ'].values #cell2mat(HM_inf(:,4));
    Fst = HM_inf['FST'].values #cell2mat(HM_inf(:,3));

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # % REMEMBER: 'n' here is INDIVIDUALS not SNPs
    # % HapMap3 Data ~1000 individuals and ~1000000 SNPs
    m = int(1e4) #number of SNPs
    n = int(1e3) #number of individuals
    pvalue = 250.0/m
    d = 3 #number of populations

    # % each row of Gamma will be populated with 3 i.i.d. draws from BN model
    # % thanks to the nature of the HapMap data (sampled from 3 populations)
    #
    # % allele freq of population: allele freq of each SNP described by that
    # % population
    G = np.zeros((m,d)); #allele freq for each population

    # % columns of S will be populated with indicator vectors s.t. each
    # % individual assigned to one of the 3 subpopulations i.e. admixture
    S = np.zeros((d,n)); #individual population admixture

    # %X = zeros(m,n); %genotype matrix
    # % random seeding here...

    random.seed(datetime.now())

    # %populate the allele freq matrix from BN with (p,F) from HapMap
    # % for each SNP...
    for i in range(0,m):
        # each row will generate 'd' variates drawing from this distribution
        G[i,:] = np.random.beta(frq[i]*(1-Fst[i])/Fst[i], (1-frq[i])*(1-Fst[i])/Fst[i], size=d)

    # print('Mean of allele freq matrix: ', np.mean(G, axis=0))

    # %populate the population admixture matrix
    # %this part is tailor made for 3 populations as described in
    # %Song et al. Nat. Genet. (2015)
    # % Treating the probabilities as ranges
    # % 1: <60/210, 2: bet. 60/210 and 120/210, 3: >=120/210

    alpha = 0.1*np.ones((d,1))
    popidx = np.zeros((n,1));
    for i in range(0,n):
        for j in range(0,d):
            S[j,i] = np.random.gamma(alpha[j],1)

        S[:,i] = S[:,i]/np.sum(S[:,i])
        I = np.argmax(S[:,i])
        popidx[i] = I+1
elif model_flag == "HGDP":
    # REMEMBER: 'n' here is INDIVIDUALS not SNPs
    # Downsampling for simulation (computationally easier)
    m = int(1e4) #number of SNPs
    n = int(305) #no. of individuals
    pvalue = 0.0025;
    flag = 0 #plot flag
    d = int(10) #number of populations (see log of --fst result)

    # allele freq of population: allele freq of each SNP described by that
    # population
    G = np.zeros((m,d)) #allele freq for each population

    # columns of S will be populated with indicator vectors s.t. each
    # individual assigned to one of the 51 subpopulations i.e. admixture
    S = np.zeros((d,n)) #individual population admixture

    # random seeding here...
    random.seed(datetime.now())

    #populate the allele freq matrix from BN with (p,F) from HapMap
    # for each SNP...
    for i in range(0,m):
        # each row will generate 'd' variates drawing from this distribution
        G[i,:] = 0.9*np.random.uniform(0, 0.5, size=d)

    # set last column to 0.05 per Song et al. 2015
    G[:,d-1] = 0.05;

    HGDP_PCs = pd.read_csv('pruned_HGDP_topPops_singVecs.txt',sep=' ',header=None)
    topPCs = HGDP_PCs.values

    for i in range(0,d):
       S[i,:] = (topPCs[:,i]-np.min(topPCs[:,i]))/(np.max(topPCs[:,i])-np.min(topPCs[:,i]))
    S[d-1,:] = 1

    popidx = np.zeros((n,1));

    HGDP_subpops = pd.read_csv('subpops_pruned_HGDP.txt',sep=' ',header=None)

    for i in range(0,HGDP_subpops.values.shape[0]):
        if HGDP_subpops.values[i] == "Biaka_Pygmies":
            popidx[i] = 1;
        elif HGDP_subpops.values[i] == "French":
            popidx[i] = 2;
        elif HGDP_subpops.values[i] == "Han":
            popidx[i] = 3;
        elif HGDP_subpops.values[i] == "Japanese":
            popidx[i] = 4;
        elif HGDP_subpops.values[i] == "Palestinian":
            popidx[i] = 5;
        elif HGDP_subpops.values[i] == "Papuan":
            popidx[i] = 6;
        elif HGDP_subpops.values[i] == "Pima":
            popidx[i] = 7;
        elif HGDP_subpops.values[i] == "Russian":
            popidx[i] = 8;
        elif HGDP_subpops.values[i] == "Sardinian":
            popidx[i] = 9;
        else:
            # Sindhi
            popidx[i] = 10;
elif model_flag == "TGP":
    # REMEMBER: 'n' here is INDIVIDUALS not SNPs
    # Downsampling for simulation (computationally easier)
    m = int(1e4) #number of SNPs
    n = int(1056) #no. of individuals
    pvalue = 0.0025;
    flag = 0 #plot flag
    d = int(10) #number of populations (see log of --fst result)

    # allele freq of population: allele freq of each SNP described by that
    # population
    G = np.zeros((m,d)) #allele freq for each population

    # columns of S will be populated with indicator vectors s.t. each
    # individual assigned to one of the 51 subpopulations i.e. admixture
    S = np.zeros((d,n)) #individual population admixture

    # random seeding here...
    random.seed(datetime.now())

    #populate the allele freq matrix from BN with (p,F) from HapMap
    # for each SNP...
    for i in range(0,m):
        # each row will generate 'd' variates drawing from this distribution
        G[i,:] = 0.9*np.random.uniform(0, 0.5, size=d)

    # set last column to 0.05 per Song et al. 2015
    G[:,d-1] = 0.05;

    TGP_PCs = pd.read_csv('pruned_TGP_topPops_singVecs.txt',sep=' ',header=None)
    topPCs = TGP_PCs.values

    for i in range(0,d):
       S[i,:] = (topPCs[:,i]-np.min(topPCs[:,i]))/(np.max(topPCs[:,i])-np.min(topPCs[:,i]))
    S[d-1,:] = 1

    popidx = np.zeros((n,1));

    TGP_subpops = pd.read_csv('subpops_pruned_TGP.txt',sep=' ',header=None)

    for i in range(0,TGP_subpops.values.shape[0]):
        if TGP_subpops.values[i] == "CHB":
            popidx[i] = 1;
        elif TGP_subpops.values[i] == "CHS":
            popidx[i] = 2;
        elif TGP_subpops.values[i] == "GIH":
            popidx[i] = 3;
        elif TGP_subpops.values[i] == "GWD":
            popidx[i] = 4;
        elif TGP_subpops.values[i] == "IBS":
            popidx[i] = 5;
        elif TGP_subpops.values[i] == "JPT":
            popidx[i] = 6;
        elif TGP_subpops.values[i] == "PUR":
            popidx[i] = 7;
        elif TGP_subpops.values[i] == "STU":
            popidx[i] = 8;
        elif TGP_subpops.values[i] == "TSI":
            popidx[i] = 9;
        else:
            # YRI
            popidx[i] = 10;
else:
    print("Not a valid model type! Options are BN, PSD, HGDP or TGP.")
    sys.exit(0)

# %Get the allele frequency matrix for each individual (SNP -by- indiv)
# % GOAL: efficiently estimate the individual-specific frequencies, F
# % Hao 2015 paper
F = np.matmul(G,S)

if model_flag == "HGDP" or "TGP":
    #####################################
    # Normalize F by column (making sure each column i.e. individual is bet. [0 1])
    F = F/F.max(axis=0)
    #####################################

# % simulating X using binomial distribution of estimated F
X = np.random.binomial(2, F)
print(X)
print(X.shape)
print(" ")

################################################################################
if sketch_flag == 1:
    X = matrix_sketch(X.T,500)
    print(X)
    print(X.shape)
    print(" ")
    m = X.shape[0]
    n = X.shape[1]
################################################################################

# # % if A is a matrix, then sum(A,2) is a column vector containing the sum of each row.
idxzer = np.where(~X.any(axis=1))[0]
# print(idxzer)
# # % randomly fill those rows with 0/1 in a random column
X[idxzer,random.randint(0,n-1)] = 1;
# print('Mean of simulated matrix: ', np.mean(X, axis=0))


# # %the second arg is related to whether we just want to mean center X or
# # %standardize by the 2*p*(1-p)
# # % same as normalize(X,2)
normX = normalize.norm(X,0);
# print(normX)
# print(normX.shape)
# sys.exit(0)

# % proportions of genetic, environmental and noise contributions
# % respectively for simulated traits (one of the 3 configs from paper)
v = [10, 0, 90]

# % simulate the traits
# % the strategy simulates non-genetic effects and random variation
traits, status = traitSim.simulate(normX,S,v,m,n,d)
Y = status

# ###############################################################################
# # Write data in cluster to plink formatted file
# if not os.path.exists('sim_plinkfiles'):
#     os.makedirs('sim_plinkfiles')
#
# new_plinkfile = convert2plink(X,Y,model_flag)
# ###############################################################################


# % Get p-values based on correlation between each individual and the status
# % (for every SNP)
# % essentially getting a measurement of how well each SNP is correlated to
# % the status
pvstat = np.zeros((m,1))
for i in range(0,m):
    pvstat[i] = ArmitChisq.chisq(normX[i,:],Y.T,n)


# Find a value that exceeds (1-0.0025)% of the samples from chi-sq with 1
# degree of freedom
# Only going to observe values greater than 'chisqst' 0.0025% of the time
chisqst = chi2.ppf(1-pvalue,df=1);
# Grab pvals that are larger than 'chisqst'
sigset = pvstat[np.where(pvstat > chisqst)[0]];
# number of candidate SNPs (1)
candSNP1 = sigset.shape[0]
# certain measurement to calculate new 'pvstat'
lam = np.median(pvstat)/0.456;
newchsq = pvstat/lam;
# number of candidate SNPs (4)
candSNP4 = pvstat[np.where(newchsq > chisqst)[0]].shape[0]

# PCA step on indivs by indivs matrix
temp = np.matmul(normX.T,normX)
U, Sig, _ = np.linalg.svd(temp);
corrPC =  np.corrcoef(U[:,0],popidx.T);
print('Correlation of top axis of variation with populations membership: ',
        corrPC[0,1]**2);


# project out the top 'K' PCs (in this case 10)
K = 10;
adjR, adjstat = EigStrat.project(normX.T,U[:,0:K],Y);

# Get p-values based on correlation between each individual and the status
# (for every SNP) *now on the adjusted data
adjpvstat = np.zeros((m,1))
for i in range(0,m):
    adjpvstat[i] = ArmitChisq.chisq(adjR[:,i],adjstat.T,n-K-1);

# Grab pvals that are larger than 'chisqst'
adjsigset = adjpvstat[np.where(adjpvstat > chisqst)[0]]
# number of candidate SNPs (2)
candSNP2 = adjsigset.shape[0]

# 1000 choose 2 pairs (calculating distance between each pair)
D = pdist(normX.T)
dele = [0, 1, 2, 3, 4]
normR = normX.T

# hierarchical clustering of individuals
# pdist(X) returns the Euclidean distance between pairs of observations in X.
clustering_time = time.time()
CS, clustcount, allidx = CluStrat.cluster(normX.T, U, D, d, Y, chisqst, dele)
end_of_clustering = time.time()
print('Time elapsed running CluStrat (in seconds): ', end_of_clustering-clustering_time)
# basically take the max "top PCs" from the clusterings and that is how
# many more candidate SNPs we have (also record for which cluster gave the max)
candSNP3 = np.max(CS)
I = np.where(CS == candSNP3)[0][0]


# if the plot flag is active...
if(plot_flag == 1):
    if model_flag == "BN" or model_flag == "PSD":
        # grab indices for cases and controls
        idxcase = np.where(status == 1)[0]
        idxcontr = np.where(status == 0)[0]

        idxA = np.where(popidx == 1)[0]
        idxcaseA = np.intersect1d(idxA,idxcase)
        idxcontrA = np.intersect1d(idxA,idxcontr)

        idxB = np.where(popidx == 2)[0]
        idxcaseB = np.intersect1d(idxB,idxcase)
        idxcontrB = np.intersect1d(idxB,idxcontr)

        idxC = np.where(popidx == 3)[0]
        idxcaseC = np.intersect1d(idxC,idxcase)
        idxcontrC = np.intersect1d(idxC,idxcontr)

        # plot top 2 PCs grouped by population membership for cases
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # print(len(idxcaseA))
        ax.scatter(U[idxcaseA, 0], U[idxcaseA, 1], marker='o', color='red', s=8, label="CaseA")
        ax.scatter(U[idxcaseB, 0], U[idxcaseB, 1], marker='o', color='blue', s=8, label="CaseB")
        ax.scatter(U[idxcaseC, 0], U[idxcaseC, 1], marker='o', color='lawngreen', s=8, label="CaseC")

        ax.scatter(U[idxcontrA, 0], U[idxcontrA, 1], marker='*', color='red', s=8, label="ControlA")
        ax.scatter(U[idxcontrB, 0], U[idxcontrB, 1], marker='*', color='blue', s=8, label="ControlB")
        ax.scatter(U[idxcontrC, 0], U[idxcontrC, 1], marker='*', color='lawngreen', s=8, label="ControlC")

        ax.legend()
        ax.set_xlabel('PC 1')
        ax.set_ylabel('PC 2')
        ax.set_title('Random SNPs');

        # Documentation to get 'plt.show()' to work:
        # https://stackoverflow.com/questions/43397162/show-matplotlib-plots-in-ubuntu-windows-subsystem-for-linux

        plt.savefig('top2PCs.png')
        # plt.show()
    # else:


end = time.time()
print('Total time elapsed (in seconds): ', end - start)
