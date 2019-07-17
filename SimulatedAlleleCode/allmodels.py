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

def convert2plink(X, status, model_flag, model_directory):

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

    file_prefix = model_directory+'/simdata_'+model_flag+'_'+str(random.randint(0, 9999999))
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

# Create 5 simulated datasets for each model and set of proportions (60 total datasets)
model_flags = ["BN","PSD","HGDP","TGP"]
# % proportions of genetic, environmental and noise contributions
# % respectively for simulated traits (one of the 3 configs from paper)
v_set = [[10, 0, 90],[20, 10, 70],[5, 5, 90]]

nums = [1,2,3,4,5]

for model_flag, v, ctr in list(itertools.product(model_flags,v_set,nums)):
    print(model_flag)
    print(v)
    print(ctr)
    print(" ")


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
    print(normX)
    # print(normX.shape)
    # sys.exit(0)

    # % simulate the traits
    # % the strategy simulates non-genetic effects and random variation
    traits, status = traitSim.simulate(normX,S,v,m,n,d)
    Y = status

    ###############################################################################
    # Write data in cluster to plink formatted file
    if v == [10, 0, 90]:
        prop = 1
    elif v == [20, 10, 70]:
        prop = 2
    else:
        prop = 3

    model_directory = 'sim_plinkfiles/'+model_flag+'/proportion'+str(prop)
    if not os.path.exists(model_directory):
        os.makedirs(model_directory)

    new_plinkfile = convert2plink(X,Y,model_flag, model_directory)
    ###############################################################################
