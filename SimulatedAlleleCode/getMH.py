import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random
from scipy.sparse.linalg import svds


def MH(X):

    n = X.shape[0]
    m = X.shape[1]
    U, _, _ = svds(np.matmul(X,X.T), k=10)
    ALS = np.matmul(U,U.T)
    #get the Leverage Scores
    d = np.diagonal(ALS)
    s = d.repeat(n)
    s.shape = (n,n)
    #calculate H_i + H_j - 2H_ij
    for i in range(len(d)):
        s[i] += d
        s[i] -= 2*ALS[i]
    #scale it
    s = (n-1)*s
    d = (n-1)*(d - 1/n)
    np.fill_diagonal(s,d)


    return s


if __name__ == '__main__':
    print("HEY HEY HEY")
