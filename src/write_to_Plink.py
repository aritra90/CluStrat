import random, time, itertools, subprocess, math, os, csv, sys, pdb
from plinkio.plinkfile import WritablePlinkFile, Sample, Locus
from datetime import datetime
from operator import add
import numpy as np
import pandas as pd

def convert2plink(X, model_directory):

    #create fake sample IDs
    # sampleIDs = [unicode(model_flag)+u"000"+unicode(s) for s in xrange(1,X.shape[1]+1)]
    #create fake rsIDs
    rsIDs = ["rs000"+str(s) for s in range(1, X.shape[0]+1)]
    #create fake positions increasing randomly between 200 and 20k centimorgans
    snpPositions = [float(s) for s in range(70000, 70000+X.shape[0])]
    for k in range(1, X.shape[0]):
        snpPositions[k] = snpPositions[k-1] + random.randint(200, 1000)

    list_of_Samples = []
    # Create samples
    for i in range(X.shape[1]):
        # print sampleIDs[i]
        num = i+1
        fid = 'id%d' % num
        iid = fid
        father_iid = '%d' % 0
        mother_iid = '%d' % 0
        list_of_Samples.append(
            Sample(fid, iid, father_iid, mother_iid, -9, -9))

    list_of_Loci = []
    # chromosome numbers based on how many SNPs being simulated
    num_chroms = 22
    chromosomes = list(range(1,num_chroms+1)) # chrs 1 thru 22
    segments = np.repeat(int(X.shape[0]/num_chroms),num_chroms)
    if m%num_chroms != 0:
        segments[0] = int(m/num_chroms) + int(m%num_chroms) # add remainder to the first chromosome
    chr_labels = np.repeat(chromosomes, segments)

    # Create loci
    for j in range(0, X.shape[0]):
        rng = random.randint(0, 1)
        if rng == 0:
            alleles = ['A', 'T']
        else:
            alleles = ['C', 'G']

        # Choosing to encode rng = 0 as allele 'A'/'C' and 2 as allele 'T'/'G' (1 as heterozygous occurrences)
        # Get the allele frequencies
        allele_counts = np.unique(X[j, :], return_counts=True)
        # print(allele_counts)
        # if 0,1,2 all occur in the SNP...
        if allele_counts[0].shape[0] == 3:
            # Determine which is the minor allele and place that allele as the first allele argument
            if allele_counts[1][0] < allele_counts[1][2]:
                list_of_Loci.append(
                    Locus(chr_labels[j], rsIDs[j], 0, snpPositions[j], alleles[0], alleles[1]))
            else:
                list_of_Loci.append(
                    Locus(chr_labels[j], rsIDs[j], 0, snpPositions[j], alleles[1], alleles[0]))
        else:
            if 0 in allele_counts[0] and 1 in allele_counts[0]:
                list_of_Loci.append(
                    Locus(chr_labels[j], rsIDs[j], 0, snpPositions[j], alleles[1], alleles[0]))
            elif 0 in allele_counts[0] and 2 in allele_counts[0]:
                # print("no occrrences of '1'...")
                if allele_counts[1][0] < allele_counts[1][1]:
                    list_of_Loci.append(
                        Locus(chr_labels[j], rsIDs[j], 0, snpPositions[j], alleles[0], alleles[1]))
                else:
                    list_of_Loci.append(
                        Locus(chr_labels[j], rsIDs[j], 0, snpPositions[j], alleles[1], alleles[0]))
            else:
                # print("no occurrences of '0'...") then its automatically the minor allele cause it doesnt show up
                list_of_Loci.append(
                    Locus(chr_labels[j], rsIDs[j], 0, snpPositions[j], alleles[0], alleles[1]))

    file_prefix = str(model_directory)+'/simdata'
    print('File prefix: '+file_prefix)

    # Create PLINK files corresponding to the data
    # Sample(fid, iid, iid, iid, sex, affection, phenotype = 0.0)
    plink_obj = WritablePlinkFile(file_prefix, list_of_Samples)
    # plink_obj = WritablePlinkFile(file_prefix,[Sample('0', '0', '0', '0', 0, 0)])

    for i in range(0, X.shape[0]):
        plink_obj.write_row(list_of_Loci[i], X[i, :])
    print("Files created")
    plink_obj.close()
    del plink_obj
    del list_of_Samples
    del list_of_Loci
    return file_prefix

# TRYING R PACKAGE BN SIMULATOR
# df = pd.read_csv("/home/mcburch/MaSk-LMM-main/simulator/BN_simulator/bn_matrix")
# df = pd.read_csv("/home/mcburch/mask_lmm_main/simulator/BN_SIM_R/bn_matrix")
df = pd.read_csv(str(sys.argv[1]))
df = df.iloc[: , 1:]
X = df.values
m, n = X.shape

# df = pd.read_csv("/home/mcburch/MaSk-LMM-main/simulator/BN_simulator/bn_trait")
# df = df.iloc[: , 1:]
# trait = df.values

# model_directory = 'bn_sim_files'
# model_directory = '1000g_sim_files'
model_directory = str(sys.argv[2])
if not os.path.exists(model_directory):
    os.makedirs(model_directory)

file_prefix = convert2plink(X, model_directory)
