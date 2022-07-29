from __future__ import division
from __future__ import absolute_import
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use(u"Agg")
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random, argparse
from operator import add
from datetime import datetime

from sklearn.cluster import KMeans
from scipy.stats.distributions import chi2
from scipy.spatial.distance import pdist
import normalize, traitSim, traitSim_risk #, ArmitChisq, EigStrat, getMH, CluStrat
from plinkio.plinkfile import WritablePlinkFile, Sample, Locus

import pdb

def msg(name=None):
    return '''data_simulate.py
         >> python data_simulate.py --model BN --prop 10,20,70 --pheno 1 --size 1000,10000
        '''

def parse_arguments():
    parser = argparse.ArgumentParser(usage=msg())

    parser.add_argument("-m", "--model", dest='model', action='store', help="Indicate which model to use for simulation (--model PSD/BN/TGP).",
                        metavar="MODEL")

    parser.add_argument("-p", "--prop", dest='prop', action='store', help="Indicate which proportion to use for simulation (--prop 1/2/3).",
                        metavar="PROP")

    parser.add_argument("-ph", "--pheno", dest='phenotype', action='store', help="Indicate binary or continuous phenotype type to use for simulation (--pheno 0/1).",
                        metavar="PHENO")

    parser.add_argument("-sz", "--size", dest='size', action='store', help="Enter the desired simulation matrix dimensions (individuals by SNPs)",
                        metavar="SIZE")

    parser.add_argument("-r", "--risk", dest='risk_flag', action='store', help="Set Risk flag if you want risk for one population to be greater than the others",
                        metavar="RISK")
    args = parser.parse_args()

    return args

def convert2plink(X, status, model_flag, plink_filenm, continuous_trait=None):

    #create fake sample IDs
    # sampleIDs = [unicode(model_flag)+u"000"+unicode(s) for s in xrange(1,X.shape[1]+1)]
    #create fake rsIDs
    rsIDs = [u"rs000"+unicode(s) for s in xrange(1,X.shape[0]+1)]
    #create fake positions increasing randomly between 200 and 20k centimorgans
    snpPositions = [float(s) for s in xrange(70000,70000+X.shape[0])]
    for k in xrange(1,X.shape[0]):
        snpPositions[k] = snpPositions[k-1] + random.randint(200,1000)

    # print(snpPositions)
    # print(' ')
    # print(sampleIDs)
    # print(" ")
    # print(rsIDs)

    if continuous_trait is None:
        list_of_Samples = []
        # Create samples
        for i in range(X.shape[1]):
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
        for i in range(X.shape[1]):
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
    for j in xrange(0,X.shape[0]):
        rng = random.randint(0,1)
        if rng == 0:
            alleles = [u'A',u'T']
        else:
            alleles = [u'C',u'G']

        # Choosing to encode rng = 0 as allele 'A'/'C' and 2 as allele 'T'/'G' (1 as heterozygous occurrences)
        # Get the allele frequencies
        allele_counts = np.unique(X[j,:], return_counts = True)
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

    print u'File prefix: '+plink_filenm

    # Create PLINK files corresponding to the data
    # pdb.set_trace()
    # Sample(fid, iid, iid, iid, sex, affection, phenotype = 0.0)
    plink_obj = WritablePlinkFile(plink_filenm,list_of_Samples)
    # plink_obj = WritablePlinkFile(file_prefix,[Sample('0', '0', '0', '0', 0, 0)])


    for i in xrange(0,X.shape[0]):
        # print(X[i,:])
        # print(X[i,:].shape)
        plink_obj.write_row(list_of_Loci[i], X[i,:])
    # plink_obj.loci = list_of_Loci

    # print("YAY")
    # sys.exit(0)
    plink_obj.close()
    del plink_obj
    del list_of_Samples
    del list_of_Loci

    # return plink_obj
if __name__ == '__main__':

    args = parse_arguments()

    if args.model:
        if args.model == "TGP":
            model_flags = [str(args.model)]
        elif args.model == "PSD":
            model_flags = [str(args.model)]
        else:
            model_flags = [u"BN"]
    else:
        print("No model argument given. Setting model to 'BN'.\n")
        model_flags = [u"BN"]

    if args.prop:
        v = args.prop.split(',')
        try:
            v_set = [int(i) for i in v]

            if len(v) != 3 or np.sum(v_set) != 100:
                print("Usage: -pr 10,20,70 means 10% genetic, 20% environmental and 70% noise variance.")
                sys.exit(1)

        except ValueError:
            print("Usage: -pr 10,20,70 means 10% genetic, 20% environmental and 70% noise variance.")
            sys.exit(1)
    else:
        print("Usage: -pr 10,20,70 means 10% genetic, 20% environmental and 70% noise variance.")
        sys.exit(1)

    if args.size:
        dims = args.size.split(',')
        try:
            sim_shape = [int(i) for i in dims]

            if len(dims) != 2:
                print("Usage: -sz 1000,10000 means 1k individuals and 10k SNPs")
                sys.exit(1)

        except ValueError:
            print("Usage: -sz 1000,10000 means 1k individuals and 10k SNPs")
            sys.exit(1)
    else:
        print("Usage: -sz 1000,10000 means 1k individuals and 10k SNPs")
        sys.exit(1)


    if args.phenotype:
        if args.phenotype == "0":
            trait_flag = [0]
        else:
            trait_flag = [1]
    else:
        print("No phenotype argument given. Setting pheno to 'continuous (1)'.\n")
        trait_flag = [0]

    if args.risk_flag:
        if args.risk_flag == "1":
            risk_flag = 1
        else:
            risk_flag = 0
    else:
        print("No risk flag argument given. Setting risk to 0.\n")
        risk_flag = 0  

    # Generates 5 simulations
    nums = [1]

    for model_flag, v, ctr, trait_flag, risk_flag in list(itertools.product(model_flags,[v_set],nums,trait_flag,risk_flag)):
        print model_flag
        print v
        print ctr

        n = int(sim_shape[0])
        m = int(sim_shape[1])


        if model_flag == u"BN":
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # %%%%%%% Load HapMap Info
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            HM_inf = pd.read_csv(u'CEUASWMEX_fst_frq.txt',sep=u' ')
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # %get allele freq and Fst for each SNP

            # % allele freq: frequency (proportion) of that SNP in the data
            # % Fst: % of contribution to the total genetic variation that each SNP has
            frq = HM_inf[u'FRQ'].values #cell2mat(HM_inf(:,4));
            Fst = HM_inf[u'FST'].values #cell2mat(HM_inf(:,3));

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            #
            # % REMEMBER: 'n' here is INDIVIDUALS not SNPs
            # % HapMap3 Data ~1000 individuals and ~1000000 SNPs
            # m = int(1e4) #number of SNPs
            # n = int(1e3) #number of individuals
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
            for i in xrange(0,m):
                # each row will generate 'd' variates drawing from this distribution
                G[i,:] = np.random.beta(frq[i]*(1-Fst[i])/Fst[i], (1-frq[i])*(1-Fst[i])/Fst[i], size=d)

            # print('Mean of allele freq matrix: ', np.mean(G, axis=0))

            # %populate the population admixture matrix
            # %this part is tailor made for 3 populations as described in
            # %Song et al. Nat. Genet. (2015)
            # % Treating the probabilities as ranges
            # % 1: <60/210, 2: bet. 60/210 and 120/210, 3: >=120/210
            popidx = np.zeros((n,1));
            for i in xrange(0,n):
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

        elif model_flag == u"PSD":
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # %%%%%%% Load HapMap Info
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            HM_inf = pd.read_csv(u'CEUASWMEX_fst_frq.txt',sep=u' ')
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # %get allele freq and Fst for each SNP

            # % allele freq: frequency (proportion) of that SNP in the data
            # % Fst: % of contribution to the total genetic variation that each SNP has
            frq = HM_inf[u'FRQ'].values #cell2mat(HM_inf(:,4));
            Fst = HM_inf[u'FST'].values #cell2mat(HM_inf(:,3));

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            #
            # % REMEMBER: 'n' here is INDIVIDUALS not SNPs
            # % HapMap3 Data ~1000 individuals and ~1000000 SNPs
            # m = int(1e4) #number of SNPs
            # n = int(1e3) #number of individuals
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
            for i in xrange(0,m):
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
            for i in xrange(0,n):
                for j in xrange(0,d):
                    S[j,i] = np.random.gamma(alpha[j],1)

                S[:,i] = S[:,i]/np.sum(S[:,i])
                I = np.argmax(S[:,i])
                popidx[i] = I+1
        elif model_flag == u"HGDP":
            # REMEMBER: 'n' here is INDIVIDUALS not SNPs
            # Downsampling for simulation (computationally easier)
            # m = int(1e4) #number of SNPs
            # n = int(305) #no. of individuals
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
            for i in xrange(0,m):
                # each row will generate 'd' variates drawing from this distribution
                G[i,:] = 0.9*np.random.uniform(0, 0.5, size=d)

            # set last column to 0.05 per Song et al. 2015
            G[:,d-1] = 0.05;

            HGDP_PCs = pd.read_csv(u'pruned_HGDP_topPops_singVecs.txt',sep=u' ',header=None)
            topPCs = HGDP_PCs.values

            for i in xrange(0,d):
               S[i,:] = (topPCs[:,i]-np.min(topPCs[:,i]))/(np.max(topPCs[:,i])-np.min(topPCs[:,i]))
            S[d-1,:] = 1

            popidx = np.zeros((n,1));

            HGDP_subpops = pd.read_csv(u'subpops_pruned_HGDP.txt',sep=u' ',header=None)

            for i in xrange(0,HGDP_subpops.values.shape[0]):
                if HGDP_subpops.values[i] == u"Biaka_Pygmies":
                    popidx[i] = 1;
                elif HGDP_subpops.values[i] == u"French":
                    popidx[i] = 2;
                elif HGDP_subpops.values[i] == u"Han":
                    popidx[i] = 3;
                elif HGDP_subpops.values[i] == u"Japanese":
                    popidx[i] = 4;
                elif HGDP_subpops.values[i] == u"Palestinian":
                    popidx[i] = 5;
                elif HGDP_subpops.values[i] == u"Papuan":
                    popidx[i] = 6;
                elif HGDP_subpops.values[i] == u"Pima":
                    popidx[i] = 7;
                elif HGDP_subpops.values[i] == u"Russian":
                    popidx[i] = 8;
                elif HGDP_subpops.values[i] == u"Sardinian":
                    popidx[i] = 9;
                else:
                    # Sindhi
                    popidx[i] = 10;
        elif model_flag == u"TGP":
            # REMEMBER: 'n' here is INDIVIDUALS not SNPs
            # Downsampling for simulation (computationally easier)
            # m = int(1e4) #number of SNPs
            # n = int(1056) #no. of individuals
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
            for i in xrange(0,m):
                # each row will generate 'd' variates drawing from this distribution
                G[i,:] = 0.9*np.random.uniform(0, 0.5, size=d)

            # set last column to 0.05 per Song et al. 2015
            G[:,d-1] = 0.05;

            TGP_PCs = pd.read_csv(u'pruned_TGP_topPops_singVecs.txt',sep=u' ',header=None)
            topPCs = TGP_PCs.values

            for i in xrange(0,d):
               S[i,:] = (topPCs[:,i]-np.min(topPCs[:,i]))/(np.max(topPCs[:,i])-np.min(topPCs[:,i]))
            S[d-1,:] = 1

            popidx = np.zeros((n,1));

            TGP_subpops = pd.read_csv(u'subpops_pruned_TGP.txt',sep=u' ',header=None)

            for i in xrange(0,TGP_subpops.values.shape[0]):
                if TGP_subpops.values[i] == u"CHB":
                    popidx[i] = 1;
                elif TGP_subpops.values[i] == u"CHS":
                    popidx[i] = 2;
                elif TGP_subpops.values[i] == u"GIH":
                    popidx[i] = 3;
                elif TGP_subpops.values[i] == u"GWD":
                    popidx[i] = 4;
                elif TGP_subpops.values[i] == u"IBS":
                    popidx[i] = 5;
                elif TGP_subpops.values[i] == u"JPT":
                    popidx[i] = 6;
                elif TGP_subpops.values[i] == u"PUR":
                    popidx[i] = 7;
                elif TGP_subpops.values[i] == u"STU":
                    popidx[i] = 8;
                elif TGP_subpops.values[i] == u"TSI":
                    popidx[i] = 9;
                else:
                    # YRI
                    popidx[i] = 10;
        else:
            print u"Not a valid model type! Options are BN, PSD, HGDP or TGP."
            sys.exit(0)

        # %Get the allele frequency matrix for each individual (SNP -by- indiv)
        # % GOAL: efficiently estimate the individual-specific frequencies, F
        # % Hao 2015 paper
        F = np.matmul(G,S)

        if model_flag == u"HGDP" or u"TGP":
            #####################################
            # Normalize F by column (making sure each column i.e. individual is bet. [0 1])
            F = F/F.max(axis=0)
            d = 10
        else:
            d = 3
            #####################################

        # % simulating X using binomial distribution of estimated F
        X = np.random.binomial(2, F)

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

        risk_pop = random.randint(d)
        # % simulate the traits
        # % the strategy simulates non-genetic effects and random variation
        if risk_flag == 0:
            traits, status = traitSim.simulate(X,S,v,m,n,d)
        else:
            traits, status = traitSim_risk.simulate(X,S,v,m,n,d,risk_pop)
        Y = traits
        # print status
        # print traits
        # print type(traits)
        # sys.exit(0)

        ###############################################################################
        # Write data in cluster to plink formatted file

        print("Converting to PLINK format...")
        if trait_flag == 0:
            plink_filenm = "simfile_"+str(args.model)+"_"+str(v[0])+"_"+str(v[1])+"_"+str(v[2])+"_"+str(args.phenotype)
            # plink_filenm = u" " #os.path.abspath(plink_filenm)

            convert2plink(X,status,model_flag, plink_filenm)
        else:
            plink_filenm = "simfile_"+str(args.model)+"_"+str(v[0])+"_"+str(v[1])+"_"+str(v[2])+"_"+str(args.phenotype)
            # plink_filenm = u" " #os.path.abspath(plink_filenm)

            convert2plink(X,status,model_flag, plink_filenm, continuous_trait=traits)

        ###############################################################################
