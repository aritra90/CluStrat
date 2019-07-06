import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random
from scipy.sparse.linalg import svds

# Transposing downstream will be an issue here
def MH(X):

    n = X.shape[1]
    U, _, _ = svds(np.matmul(X,X.T), k=5)
    ALS = np.matmul(U,U.T)
    MH = (n-1)*(ALS - 1/n)

    return MH


if __name__ == '__main__':
    print("HEY HEY HEY")
