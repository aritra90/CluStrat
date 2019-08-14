import numpy as np
import numpy.matlib
import math, sys, argparse, subprocess
from plinkio import plinkfile
import parseplink as pp

def get_weights(data):

    num_indivs = len(data)
    num_SNPs = data.shape[1]

    # Calculate the weights based on MAF
    # Minor allele = the allele that occurs the least (not necessarily dominant/recessive)
    # So, count the number of dominant/recessive alleles. Determine which is the minor and then do the frequency
    weights = np.arange(0,num_SNPs)
    weights = weights.astype(float)
    frequencies = np.arange(0,num_SNPs)
    frequencies = frequencies.astype(float)

    for elem in np.arange(0,num_SNPs):

        dominant_allele = 0
        recessive_allele = 0
        total_alleles = 0

        unique, counts = np.unique(data[:,elem], return_counts=True)
        allele_counts = dict(zip(unique,counts))
        # print(allele_counts)

        for x in allele_counts:
            if x == 1.0:
                dominant_allele = dominant_allele + 2*allele_counts[x]
                # total_alleles = total_alleles + 2*allele_counts[x]
            elif x == 2.0:
                dominant_allele = dominant_allele + allele_counts[x]
                recessive_allele = recessive_allele + allele_counts[x]
                # total_alleles = total_alleles + 2*allele_counts[x]
            elif x == 3.0:
                recessive_allele = recessive_allele + 2*allele_counts[x]
                # total_alleles = total_alleles + 2*allele_counts[x]
        total_alleles = dominant_allele + recessive_allele

        if dominant_allele < recessive_allele:
            # pk
            maf = dominant_allele/float(total_alleles)
            if maf == 0.0:
                SNP_weight = 0.0
            else:
                SNP_weight = 1/(math.sqrt(maf*(1-maf)))
        else:
            maf = recessive_allele/float(total_alleles)
            if maf == 0.0:
                SNP_weight = 0.0
            else:
                SNP_weight = 1/(math.sqrt(maf*(1-maf)))

        frequencies[elem] = maf
        weights[elem] = SNP_weight

    return weights

def kernelize(data1, data2, weights):

    # Size of first matrix
    num_indivs1 = len(data1)
    num_SNPs1 = data1.shape[1]

    # Mean center the data about the columns
    means1 = np.nanmean(data1,axis=0)
    data1 = data1 - means1

    # Size of second matrix
    num_indivs2 = len(data2)
    num_SNPs2 = data2.shape[1]

    # Mean center the data about the columns
    means2 = np.nanmean(data2,axis=0)
    data2 = data2 - means2

    # Make the kernel (dot product between rows)
    # Resultant matrix should be num_indivs x num_indivs
    # Genetic relationship matrix (covariance of individuals)
    k_data = np.inner(data1*weights, data2)

    # Quadratic kernel
    # K_quadratic = np.square(1 + data1*weights)

    # Gaussian kernel (source of the delta value???)
    # K_gaussian = SUBTRACTION LOOP???

    # Alike-instate kernel
    # K_alike = SUBTRACTION LOOP???

    return k_data

if __name__ == '__main__':

    # Generate the dummy training data
    num_indivs = 3
    num_SNPs = 10
    train_data = np.random.randint(3,size=(num_indivs,num_SNPs))
    train_data = train_data + 1
    train_data = train_data.astype(float)

    # Generate the dummy testing data
    num_indivs = 5
    num_SNPs = 10
    test_data = np.random.randint(3,size=(num_indivs,num_SNPs))
    test_data = test_data + 1
    test_data = test_data.astype(float)


    the_weights = get_weights(train_data)
    K = kernelize(test_data,train_data,the_weights)
