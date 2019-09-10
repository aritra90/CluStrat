from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random
from datetime import datetime
import statsmodels.api as sm
from scipy.stats.distributions import chi2
from scipy.spatial.distance import pdist
from sklearn.model_selection import train_test_split
import normalize, traitSim, ArmitChisq, EigStrat, CluStrat, getMH

if len(sys.argv) == 1:
    model_flag = "BN"
else:
    model_flag = str(sys.argv[1])
# % proportions of genetic, environmental and noise contributions
# % respectively for simulated traits (one of the 3 configs from paper)
varflag = sys.argv[2]

if int(varflag) == 1:
	v = [5, 5, 90]
elif int(varflag) == 2:
	v = [10, 0, 90]
else:
	v = [10, 20, 70]

prop = ''.join(str(e) for e in v)

print("Model type: ", model_flag)
sketch_flag = 1 #useless
plot_flag = 0
#change number of iteration here
NUMRUN = 1

pdl = []
#Open Output files
of = open(model_flag + "_" + prop+"_SP_out.txt","a+")
of.write("CHISQ\tEIGENSTRAT\tCluStrat\n")

of1 = open(model_flag + "_" + prop+"_CS_out.txt","a+")
of1.write("CHISQ\tEIGENSTRAT\tCluStrat\n")

################################################################################
# Variable SNPs for timing purposes
# snp_list = [5000, 10000, 15000]
# snp_list = [5e3, 10e3, 15e3]
snp_list = [5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000]
# snp_list = list(np.linspace(5000,500000,19).astype(int))
################################################################################

start = time.time()

# Timing in seconds...
time_D = []
time_dele4 = []
time_dele5 = []
time_dele6 = []
time_clustrat = []

for num_snps in snp_list:

    for iter in range(NUMRUN):
    	st1 = time.time()
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
    		m = num_snps #number of SNPs
    		n = int(1e3) #number of individuals
    		pvalue = 25.0/m
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
    		m = num_snps #number of SNPs
    		n = int(1e3) #number of individuals
    		pvalue = 25.0/m
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
    		m = num_snps #number of SNPs
    		n = int(5e2) #no. of individuals
    		pvalue = 25.0/m;
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
    		m = num_snps #number of SNPs
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
    	# print(X)


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

    	# % simulate the traits
    	# % the strategy simulates non-genetic effects and random variation
    	traits, status = traitSim.simulate(normX,S,v,m,n,d)
    	Y = traits

    	en1 = time.time()
    	print('Time to get data : ', en1 - st1)

    	st2 = time.time()
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
    	sigset = np.where(pvstat > chisqst)[0];
    	#get the significant SNPs in test data


    	# number of candidate SNPs (1)
    	candSNP1 = sigset.shape[0]
    	#print("Candidate SNP after Armitage trend Chisq: " + str(candSNP1))
    	sig1 = [i for i in sigset if i <= 10]
    	print("\n ================= \n ")
    	print("number of causal SNPs in Armitage trend chisq = " + str(len(sig1)))
    	print("Number of spurious associations " + str(candSNP1 - len(sig1)))
    	print("\n ================= \n ")
    	# SP1.append(candSNP1 - len(sig1))
    	# CS1.append(len(sig1))
    	SP1 = candSNP1 - len(sig1)
    	CS1 = len(sig1)
    	pdl.append({'Spurious': SP1, 'Causal': CS1,
    				'Method': 'Armitage trend $\chi^2$'})
    	en2 = time.time()

    	print('Time to ArmitChisq : ', en2 - st2)

    	# # certain measurement to calculate new 'pvstat'
    	# lam = np.median(pvstat)/0.456;
    	# newchsq = pvstat/lam
    	# # number of candidate SNPs (4)
    	# candSNP4 = pvstat[np.where(newchsq > chisqst)[0]].shape[0]

    	# PCA step on indivs by indivs matrix
    	st3 = time.time()
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

    	# Get p-values based on correlation between each individual and the status
    	# (for every SNP) *now on the adjusted data
    	#Radjpvstat = CluStrat1.ridge_pvals(adjR, adjstat)
    	#Radjsigset = np.where(Radjpvstat < pvalue)[0]
    	#RcandSNP2 = Radjsigset.shape[0]
    	# Rsig2 = [i for i in Radjsigset if i <= 10]
    	# print("number of causal SNPs in Eigen Ridge = " + str(len(Rsig2)))
    	# print("Number of spurious associations " + str(RcandSNP2 - len(Rsig2)))
    	# print("\n ================= \n ")
    	# SP2.append(RcandSNP2 - len(Rsig2))
    	# CS2.append(len(Rsig2))
    	# SP2 = RcandSNP2 - len(Rsig2)
    	# CS2 = len(Rsig2)

    	pdl.append({'Spurious': SP2, 'Causal': CS2,
    				'Method': 'EIGENSTRAT'})
    	en3 = time.time()
    	#print('Time to EIGENSTRAT : ', en1 - st1)

    	# 1000 choose 2 pairs (calculating distance between each pair)
    	#D = pdist(normX.T)
    	st4 = time.time()

        ########################################################################
        # TIMING CALCULATION OF D
    	stD = time.time()
    	D = getMH.MH(normX.T)
    	enD = time.time()
    	time_D.append(enD - stD)
    	print(D)
    	dele_list = [4, 5, 6]
    	normR = normX.T
    	# hierarchical clustering of individuals
    	# pdist(X) returns the Euclidean distance between pairs of observations in X.
    	clustering_time = time.time()


        ########################################################################
        # TIMING CLUSTRAT ACROSS VALUES OF 'dele'
    	CS = []
    	SP = []
    	for dele in dele_list:
            stdele = time.time()
            CS_temp, clustcount, SP_temp = CluStrat.cluster(X.T, U, D, d, Y, pvalue, [dele])
            CS.append(CS_temp)
            SP.append(SP_temp)
            endele = time.time()
            if dele == 4:
                time_dele4.append(endele - stdele)
            elif dele == 5:
                time_dele5.append(endele - stdele)
            else:
                time_dele6.append(endele - stdele)

    	print(CS)

    	end_of_clustering = time.time()
    	print('Time elapsed running CluStrat (in seconds): ', end_of_clustering-clustering_time)
    	time_clustrat.append(end_of_clustering-clustering_time)
    	# # basically take the max "top PCs" from the clusterings and that is how
    	# # many more candidate SNPs we have (also record for which cluster gave the max)
    	maxidx = CS.index(max(CS))
    	print(max(CS))
    	print(maxidx)
    	SPmax = SP[maxidx]
    	# print("Candidate SNP after CluStrat: " + str(SPmax))
    	# # SP3.append(SPmax)
    	# # CS3.append(CS[maxidx])
    	SP3 = SPmax
    	CS3 = CS[maxidx]
    	en4 = time.time()

    	print('Time to CluStrat : ', en4 - st4)

    	#pdl contains a direct way of input for box plots if we want to integrate the code later.
    	#might be of use later.
    	pdl.append({'Spurious': SP3[0], 'Causal': CS3[0],
    				'Method': 'CluStrat'})
    	#Write to output
    	of.write("%d \t %d \t %d \n" %(SP1,SP2,SP3))
    	of1.write("%d \t %d \t %d \n" %(CS1,CS2,CS3))


of.close()
of1.close()

df = pd.DataFrame(pdl)
#print(df)
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
        plt.style.use('seaborn-darkgrid')
        # plot top 2 PCs grouped by population membership for cases
        fig = plt.figure()
        ax = fig.add_subplot(111)
        scale = 200.0 * np.random.rand(750)

        ax.scatter(U[idxA, 0], U[idxA, 1], s=scale, marker='*', color='red', alpha = 0.6, label="A")
        ax.scatter(U[idxB, 0], U[idxB, 1], s=scale, marker='*', color='blue', alpha = 0.6, label="B")
        ax.scatter(U[idxC, 0], U[idxC, 1], s=scale, marker='*', color='lawngreen', alpha = 0.6, label="C")

        # print(len(idxcaseA))
        # ax.scatter(U[idxcaseA, 0], U[idxcaseA, 1], marker='o', color='red', s=10, label="A")
        # ax.scatter(U[idxcaseB, 0], U[idxcaseB, 1], marker='o', color='blue', s=10, label="B")
        # ax.scatter(U[idxcaseC, 0], U[idxcaseC, 1], marker='o', color='lawngreen', s=10, label="C")

        # ax.scatter(U[idxcontrA, 0], U[idxcontrA, 1], marker='*', color='red', s=10, label="A")
        # ax.scatter(U[idxcontrB, 0], U[idxcontrB, 1], marker='*', color='blue', s=10, label="B")
        # ax.scatter(U[idxcontrC, 0], U[idxcontrC, 1], marker='*', color='lawngreen', s=10, label="C")

        ax.legend()
        ax.set_xlabel('PC 1')
        ax.set_ylabel('PC 2')
        #ax.set_title('Random SNPs');

        # Documentation to get 'plt.show()' to work:
        # https://stackoverflow.com/questions/43397162/show-matplotlib-plots-in-ubuntu-windows-subsystem-for-linux

        plt.savefig(model_flag + '_top2PCs.png', bbox_inches='tight', dpi=600)
        # plt.show()
    # else:


# PLOT THE TIMING FIGURES
plt.hlines(time_D, 0, snp_list, linestyle="dashed")
plt.plot(snp_list, time_D,'g')
plt.xlabel('# of SNPs')
plt.ylabel('Time (in seconds)')
plt.title('Time to compute D')
plt.savefig(str(model_flag)+'_timeD_prop'+str(varflag)+'.png')
plt.clf()
plt.close()

plt.hlines(time_dele4, 0, snp_list, linestyle="dashed")
plt.plot(snp_list, time_dele4,'r')
plt.xlabel('# of SNPs')
plt.ylabel('Time (in seconds)')
plt.title('CluStrat time (dele = 4)')
plt.savefig(str(model_flag)+'_timedele4_prop'+str(varflag)+'.png')
plt.clf()
plt.close()

plt.hlines(time_dele5, 0, snp_list, linestyle="dashed")
plt.plot(snp_list, time_dele5,'c')
plt.xlabel('# of SNPs')
plt.ylabel('Time (in seconds)')
plt.title('CluStrat time (dele = 5)')
plt.savefig(str(model_flag)+'_timedele5_prop'+str(varflag)+'.png')
plt.clf()
plt.close()

plt.hlines(time_dele6, 0, snp_list, linestyle="dashed")
plt.plot(snp_list, time_dele6,'k')
plt.xlabel('# of SNPs')
plt.ylabel('Time (in seconds)')
plt.title('CluStrat time (dele = 6)')
plt.savefig(str(model_flag)+'_timedele6_prop'+str(varflag)+'.png')
plt.clf()
plt.close()

plt.hlines(time_clustrat, 0, snp_list, linestyle="dashed")
plt.plot(snp_list, time_clustrat,'m')
plt.xlabel('# of SNPs')
plt.ylabel('Time (in seconds)')
plt.title('Total CluStrat time')
plt.savefig(str(model_flag)+'_timeclustrat_prop'+str(varflag)+'.png')
plt.clf()
plt.close()


end = time.time()
print('Total time elapsed (in seconds): ', end - start)
