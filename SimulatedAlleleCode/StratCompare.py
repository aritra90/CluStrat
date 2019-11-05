#######################################################################
## This generates datasets for 4 methods: BN, PSD, HGDP, TGP as 
## shown in Song et al. Nature (2015). 
## The code then performs a comparison of stratification correction
## techniques. Namely, uncorrected stratification as baseline with 
## Armitage trend chi-sq. Then, the two popular techniques: 
## PCA based and LMM based (GEMMA, EMMAX). Thereafter, it runs 
## CluStrat, a structure informed clustering based stratification 
## technique which outperform all the other techniques in identifying 
## the maximum number of causal SNPs while allowing for some spurious 
## associations. 

### Run: python or python3 StratCompare.py 1 
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

from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pdb
import sys, csv, os, math, subprocess, itertools, time, random
from operator import add
from datetime import datetime
from subprocess import Popen, PIPE
from sklearn.cluster import KMeans
from scipy.stats.distributions import chi2
from scipy.spatial.distance import pdist
from sklearn.model_selection import train_test_split
import normalize, traitSim , ArmitChisq, EigStrat, getMH, CluStrat
from plinkio.plinkfile import WritablePlinkFile, Sample, Locus

#######################################################################
###################   DEFINE VARIABLES      ###########################
#######################################################################
PATH_TO_GEMMA = '/depot/pdrineas/data/gwas_software/gemma/gemma-0.98.1' 
PATH_TO_EMMAX_KIN = '/depot/pdrineas/data/gwas_software/emmax/emmax-kin'
PATH_TO_EMMAX = '/depot/pdrineas/data/gwas_software/emmax/emmax'
PRE_PATH = subprocess.call(['pwd']) 
NUMRUN = 100
# Create 5 simulated datasets for each model and set of proportions (60 total datasets)
# model_flags = ["BN"]
model_flags = ["BN","PSD1","PSD01","PSD5","TGP"]
# % proportions of genetic, environmental and noise contributions
# % respectively for simulated traits (one of the 3 configs from paper)
v_set = [[5,5,90],[10, 0, 90],[20, 10, 70]]
nums = np.arange(1)
#nums = [1,2,3,4,5]
trait_flag = 0
start = time.time()


#######################################################################

def convert2plink(X, status, model_flag, model_directory, continuous_trait=None):

    #create fake sample IDs
    # sampleIDs = [unicode(model_flag)+u"000"+unicode(s) for s in xrange(1,X.shape[1]+1)]
    #create fake rsIDs
    rsIDs = ["rs000"+str(s) for s in range(1,X.shape[1]+1)]
    #create fake positions increasing randomly between 200 and 20k centimorgans
    snpPositions = [float(s) for s in range(70000,70000+X.shape[1])]
    for k in range(1,X.shape[1]):
        snpPositions[k] = snpPositions[k-1] + random.randint(200,1000)

    # print(snpPositions)
    # print(' ')
    # print(sampleIDs)
    # print(" ")
    # print(rsIDs)

    if continuous_trait is None:
        list_of_Samples = []
        # Create samples
        for i in range(X.shape[0]):
            # print sampleIDs[i]
            num = i+1
            fid = '%d'%num
            iid = fid
            father_iid = fid
            mother_iid = fid
            list_of_Samples.append(Sample(fid, iid, father_iid, mother_iid, -9, status[i]))

        # sys.exit(0)
        # print(list_of_Samples)
        # print(len(list_of_Samples))
    else:
        list_of_Samples = []
        # Create samples
        for i in range(X.shape[0]):
            # print sampleIDs[i]
            num = i+1
            fid = '%d'%num
            iid = fid
            father_iid = fid
            mother_iid = fid
            list_of_Samples.append(Sample(fid, iid, father_iid, mother_iid, -9, 2, continuous_trait[i]))

        # sys.exit(0)
        # print(list_of_Samples)
        # print(len(list_of_Samples))

    list_of_Loci = []
    # Create loci
    for j in range(0,X.shape[1]):
        rng = random.randint(0,1)
        if rng == 0:
            alleles = ['A','T']
        else:
            alleles = ['C','G']

        # Choosing to encode rng = 0 as allele 'A'/'C' and 2 as allele 'T'/'G' (1 as heterozygous occurrences)
        # Get the allele frequencies
        allele_counts = np.unique(X[:,j], return_counts = True)
        # print(allele_counts)
        # if 0,1,2 all occur in the SNP...
        if allele_counts[0].shape[0] == 3:
            # Determine which is the minor allele and place that allele as the first allele argument
            if allele_counts[1][0] < allele_counts[1][2]:
                list_of_Loci.append(Locus(1, rsIDs[j], 0, snpPositions[j], alleles[0], alleles[1]))
            else:
                list_of_Loci.append(Locus(1, rsIDs[j], 0, snpPositions[j], alleles[1], alleles[0]))
        else:
            if 0 in allele_counts[0] and 1 in allele_counts[0]:
                # print("no occurrences of '2'...")
                list_of_Loci.append(Locus(1, rsIDs[j], 0, snpPositions[j], alleles[1], alleles[0]))
            elif 0 in allele_counts[0] and 2 in allele_counts[0]:
                # print("no occrrences of '1'...")
                if allele_counts[1][0] < allele_counts[1][1]:
                    list_of_Loci.append(Locus(1, rsIDs[j], 0, snpPositions[j], alleles[0], alleles[1]))
                else:
                    list_of_Loci.append(Locus(1, rsIDs[j], 0, snpPositions[j], alleles[1], alleles[0]))
            else:
                # print("no occurrences of '0'...") then its automatically the minor allele cause it doesnt show up
                list_of_Loci.append(Locus(1, rsIDs[j], 0, snpPositions[j], alleles[0], alleles[1]))
            # sys.exit(0)

    #file_prefix = model_directory+'/simdata_'+model_flag+'_'+str(random.randint(0, 9999999))
        if continuous_trait is None:
                file_prefix = model_directory + '/simdata_'+model_flag+'_'+'bin'
        else: 
                file_prefix = model_directory + '/simdata_'+model_flag+'_'+'cont'
    print('File prefix: '+file_prefix)

    # Create PLINK files corresponding to the data
    # pdb.set_trace()
    # Sample(fid, iid, iid, iid, sex, affection, phenotype = 0.0)
    plink_obj = WritablePlinkFile(file_prefix,list_of_Samples)
    # plink_obj = WritablePlinkFile(file_prefix,[Sample('0', '0', '0', '0', 0, 0)])

    for i in range(0,X.shape[1]):
        # print(X[i,:])
        # print(X[i,:].shape)
        plink_obj.write_row(list_of_Loci[i], X[:,i])
    # plink_obj.loci = list_of_Loci
    print("Files created")
    # print("YAY")
    # sys.exit(0)
    plink_obj.close()
    del plink_obj
    del list_of_Samples
    del list_of_Loci
    #subprocess.call(['rm', 'file_prefix*'])
    #print("Deleting files as created") 
    return file_prefix
    # return plink_obj

### check here 
trait_flag = sys.argv[1]

for iter in range(NUMRUN):
    start_epoch = time.time()
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

        elif model_flag == "PSD1":
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
        elif model_flag == "PSD01":
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

            alpha = 0.01*np.ones((d,1))
            popidx = np.zeros((n,1));
            for i in range(0,n):
                for j in range(0,d):
                    S[j,i] = np.random.gamma(alpha[j],1)

                S[:,i] = S[:,i]/np.sum(S[:,i])
                I = np.argmax(S[:,i])
                popidx[i] = I+1
        elif model_flag == "PSD5":
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

            alpha = 0.5*np.ones((d,1))
            popidx = np.zeros((n,1));
            for i in range(0,n):
                for j in range(0,d):
                    S[j,i] = np.random.gamma(alpha[j],1)

                S[:,i] = S[:,i]/np.sum(S[:,i])
                I = np.argmax(S[:,i])
                popidx[i] = I+1
		
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
        # will show num_snps x num_indivs
        print(X.shape)
        X = X.T
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
        normX,_ = normalize.norm(X,0);
        print(normX)
        # print(normX.shape)
        # sys.exit(0)

        # % simulate the traits
        # % the strategy simulates non-genetic effects and random variation
        traits, status = traitSim.simulate(normX.T,S,v,m,n,d)
        if trait_flag == 1:
                Y = status
        else:
                Y = traits
        # print status
        # print traits
        # print type(traits)
        # sys.exit(0)

        ###############################################################################
        # Write data in cluster to plink formatted file
        if v == [10, 0, 90]:
            prop = 1
        elif v == [20, 10, 70]:
            prop = 2
        else:
            prop = 3

        if trait_flag == 1:
            model_directory = 'sim_plinkfiles/'+model_flag+'/proportion'+str(prop)+'/binary'
            if not os.path.exists(model_directory):
                os.makedirs(model_directory)

            # model_directory = u" " #os.path.abspath(model_directory)

            file_prefix = convert2plink(X,status,model_flag, model_directory)
        else:
            model_directory = 'sim_plinkfiles/'+model_flag+'/proportion'+str(prop)+'/continuous'
            if not os.path.exists(model_directory):
                os.makedirs(model_directory)

            # model_directory = u" " #os.path.abspath(model_directory)

            file_prefix = convert2plink(X,status,model_flag, model_directory, continuous_trait=traits)
        ###############################################################################
        ofSP = open(model_flag + "_" + str(prop) + "_" + str(trait_flag) + "_SP_out.txt", "a+") 
        ofCS = open(model_flag + "_" + str(prop) + "_" + str(trait_flag) +  "_CS_out.txt", "a+")
        # ###############################################################################
        # ####################      RUN ARMITAGE TREND CHISQ          ###################
        # ###############################################################################
        st2 = time.time()
        # % Get p-values based on correlation between each individual and the status
        # % (for every SNP)
        # % essentially getting a measurement of how well each SNP is correlated to
        # % the status
        pvstat = np.zeros((m,1))
        for i in range(0,m):
            pvstat[i] = ArmitChisq.chisq(normX[:,i],Y.T,n)
        # Find a value that exceeds (1-0.0025)% of the samples from chi-sq with 1
        # degree of freedom
        # Only going to observe values greater than 'chisqst' 0.0025% of the time
        chisqst = chi2.ppf(1-pvalue,df=1);
        # Grab pvals that are larger than 'chisqst'
        sigset = np.where(pvstat > chisqst)[0];
        #get the significant SNPs in test data     
        # number of candidate SNPs (1)    
        candSNP1 = sigset.shape[0]
        #print("Candidate SNP after Armitage trend Chisq: " + str(candSNP1))
        sig1 = [i for i in sigset if i <= 10]
        print("\n ================= \n ")
        print("number of causal SNPs in Armitage trend chisq = " + str(len(sig1)))
        print("Number of spurious associations " + str(candSNP1 - len(sig1)))

        
        SP1 = candSNP1 - len(sig1)
        CS1 = len(sig1)
        en2 = time.time()
       
        print('Time to ArmitChisq : ', en2 - st2)
        print("\n ================= \n ")
        ###############################################################################
        ####################           RUN EIGSTRAT                 ###################
        ###############################################################################
        ##PCA step on indivs by indivs matrix
        st3 = time.time()
        temp = np.matmul(normX,normX.T)
        U, Sig, _ = np.linalg.svd(temp);
        corrPC =  np.corrcoef(U[:,0],popidx.T);
        print('Correlation of top axis of variation with populations membership: ',
                corrPC[0,1]**2);
        print("\n ================= \n ")
        # project out the top 'K' PCs (in this case 10)
        K = 10;
        adjR, adjstat = EigStrat.project(normX,U[:,0:K],Y);
        # Get p-values based on correlation between each individual and the status
        # (for every SNP) *now on the adjusted data
        adjpvstat = np.zeros((m,1))
        for i in range(0,m):
            adjpvstat[i] = ArmitChisq.chisq(adjR[:,i],adjstat.T,n-K-1);
        # Grab pvals that are larger than 'chisqst'
        adjsigset = np.where(adjpvstat > chisqst)[0]
        # number of candidate SNPs (2)
        candSNP2 = adjsigset.shape[0]
        sig2 = [i for i in adjsigset if i <= 10]
        print("\n ================= \n ")
        print("number of causal SNPs in EigStrat = " + str(len(sig2)))
        print("Number of spurious associations " + str(candSNP2 - len(sig2)))
        SP2 = candSNP2 - len(sig2) 
        CS2 = len(sig2)
        
        en3 = time.time()

        print('Time to EIGENTRAT : ', en3 - st3) 
        print("\n ================= \n ")
        ##############################################################################
        ####################           RUN GEMMA                    ###################
        ###############################################################################
        st4 = time.time()
        GRM_PATH = str(os.path.basename(file_prefix))+'_'+str(prop)+'_GRM'
        ASSOC_PATH = str(os.path.basename(file_prefix))+'_'+str(prop)+'_ASSOC'
        subprocess.call([PATH_TO_GEMMA, '-bfile', file_prefix, '-gk', '1', '-o', GRM_PATH])
        GRM_PATH = 'output/'+str(GRM_PATH)+'.cXX.txt'
        subprocess.call([PATH_TO_GEMMA, '-bfile', file_prefix, '-k', \
            GRM_PATH, '-lmm', '2', '-o', ASSOC_PATH])    
        ASSOC_PATH = 'output/'+str(ASSOC_PATH)+'.assoc.txt'
        SIGSNP_PATH = 'output/'+str(os.path.basename(file_prefix))+ '_sigsnp.txt'
        f1 = open(SIGSNP_PATH, "w")
        subprocess.call(['awk', '$NF < 0.0025{print NR}', ASSOC_PATH], stdout=f1)
        f1.close()
        p1 = Popen(['wc', '-l', SIGSNP_PATH], stdout=subprocess.PIPE)
        candSNP3 = p1.communicate()[0].split(' ')  
        candSNP3 = int(candSNP3[0])
        p1.stdout.close() 
        p2 = Popen(['awk', '$1 < 11{print $1}', SIGSNP_PATH], stdout=subprocess.PIPE)
        p3 = Popen(['wc', '-l'], stdin=p2.stdout, stdout=subprocess.PIPE)
        p2.stdout.close()
        sig3 = p3.communicate()[0].split(' ')
        sig3 = int(sig3[0])
        SP3 = candSNP3 - sig3  
        CS3 = sig3
        p3.stdout.close()
        print("\n ================= \n ")
        print("number of causal SNPs in GEMMA = " + str(sig3))
        print("Number of spurious associations =  " + str(candSNP3 - sig3))
        en4 = time.time() 
        print('Time to GEMMA : ', en4 - st4) 
        print("\n ================= \n ")
        ###############################################################################
        ####################           RUN GEMMA                    ###################
        ###############################################################################
        st5 = time.time()
        subprocess.call(['/depot/pdrineas/data/gwas_software/plink2/plink', '--bfile', file_prefix, '--recode', '12', '--output-missing-genotype', '0', '--transpose','--out', file_prefix])
        TFAM_PATH = file_prefix + '.tfam'
        TMP_TFAM_PATH = file_prefix + '_tmp.tfam'
        f2 = open(TMP_TFAM_PATH, "w")
        subprocess.call(['awk', '{print $1, $2, $NF}', TFAM_PATH], stdout=f2)
        f2.close()
        subprocess.call([PATH_TO_EMMAX_KIN, '-v', '-d', '10', file_prefix])
        KIN_PATH = str(file_prefix)+'.BN.kinf'
        EMMAXOUT_PATH = 'output/'+str(os.path.basename(file_prefix))+'_'+str(prop)+'_EMMAX'
        subprocess.call([PATH_TO_EMMAX, '-v', '-d', '10', '-t', file_prefix, '-p', TMP_TFAM_PATH, '-k', KIN_PATH, '-o', EMMAXOUT_PATH])
        EMMAXOUT_PATH = EMMAXOUT_PATH + '.ps'
        
        p4 = Popen(['awk', '$NF < 0.0025{print NR}', EMMAXOUT_PATH], stdout=subprocess.PIPE)
        p7 = Popen(['wc', '-l'], stdin = p4.stdout, stdout=subprocess.PIPE)     
        candSNP4 = p7.communicate()[0].split(' ')  
        candSNP4 = int(candSNP4[0])
        p7.stdout.close()
        p5 = Popen(['awk', '$1 < 11{print $1}'], stdin=p4.stdout, stdout=subprocess.PIPE)
        p6 = Popen(['wc', '-l'], stdin=p5.stdout, stdout=subprocess.PIPE)
        p4.stdout.close()
        p5.stdout.close()
        sig4 = p6.communicate()[0].split(' ')
        sig4 = int(sig4[0])
        p6.stdout.close()
        SP4 = candSNP4 - sig4  
        CS4 = sig4
        print("\n ================= \n ")
        print("number of causal SNPs in EMMAX = " + str(sig4))
        print("Number of spurious associations =  " + str(candSNP4 - sig4))
        en5 = time.time() 
        print('Time to EMMAX : ', en5 - st5) 
        print("\n ================= \n ")
        ###############################################################################
        ####################           RUN CLUSTRAT                   ###################
        ###############################################################################
        
        # 1000 choose 2 pairs (calculating distance between each pair)
        #D = pdist(normX.T)	
        st6 = time.time() 
        D = getMH.MH(normX)
        print(D)
        sketch_flag = 0	
        if model_flag == "TGP":
            dele = [0]
        else: 
            dele = [3, 5]
        # hierarchical clustering of individuals
        # pdist(X) returns the Euclidean distance between pairs of observations in X.
        clustering_time = time.time()
        CS, clustcount, SP = CluStrat.cluster(X, D, d, Y, pvalue, dele, sketch_flag)
        CS = list(CS)
        SP = list(SP)
        end_of_clustering = time.time()
        print('Time elapsed running CluStrat (in seconds): ', end_of_clustering-clustering_time)
        # # basically take the max "top PCs" from the clusterings and that is how
        # # many more candidate SNPs we have (also record for which cluster gave the max)
        maxidx = CS.index(max(CS))
        SPmax = SP[maxidx]
        # print("Candidate SNP after CluStrat: " + str(SPmax))
        # # SP3.append(SPmax)
        # # CS3.append(CS[maxidx])
        SP5 = SPmax
        CS5 = CS[maxidx]
        en6 = time.time()
        print('Time to CluStrat : ', en6 - st6)
        ofSP.write("%d \t %d \t %d \t %d \t %d \n" %(SP1, SP2, SP3, SP4, SP5))
        ofCS.write("%d \t %d \t %d \t %d \t %d \n" %(CS1, CS2, CS3, CS4 ,CS5))
        #ofSP.write("%d \n" %(SP5))
        #ofCS.write("%d \n" %(CS5))
        print("One epoch took (mins): ", (time.time() - start_epoch)/60.0)
# Putting headers
fline = "CHISQ\tEIGENSTRAT\tGEMMA\tEMMAX\tCluStrat\n"    
oline = ofCS.readlines()
oline.insert(0,fline) 
oline1 = ofSP.readlines()
oline1.insert(0,fline)
ofCS.close()
ofSP.close() 

end = time.time()
print('Total time elapsed (in seconds): ', end - start)

