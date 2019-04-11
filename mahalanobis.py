import numpy as np
import numpy.matlib
import math, sys, argparse, subprocess
from plinkio import plinkfile
import parseplink as pp

def kernelize(data1, S_inv, data2):

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

    k_data = np.dot(np.dot(data1,S_inv),data2.T)

    return k_data

def MDist(data1):

    # Size of first matrix
    num_indivs1 = len(data1)
    num_SNPs1 = data1.shape[1]

    # Mean center the data about the columns
    means1 = np.nanmean(data1,axis=0)
    data1 = data1 - means1

    # Covariance of SNPs (encoding the LD correlation and structure into the kernel)
    # Treats 'S' as our weights
    S = np.inner(data1.T, data1.T)

    print(S)
    print(" ")
    # Zero determinant --> singular matrix so add noise to S? Regularization?
    print("Determinant of S: ", np.linalg.det(S))
    # Large condition number corresponds to a matrix being difficult to invert
    # (takes a while to compute so use determinant as red flag)
    # print("Condition No. of S: ", np.linalg.cond(S))

    if np.linalg.det(S) == 0.0:
        noise = np.random.normal(0,1,num_SNPs1*num_SNPs1)
        noise = noise.reshape(num_SNPs1, num_SNPs1)
        # print(noise)
        S = S + noise

    S_inv = np.linalg.inv(S)

    return S_inv

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

    S_inv = MDist(train_data)
    K = kernelize(test_data, S_inv, train_data)

    print(K)
