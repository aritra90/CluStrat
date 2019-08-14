import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random

# Transposing downstream will be an issue here
def norm(X, flag):
    R = X.T.astype(float)
    # print(R)
    if flag == 1:
        for i in range(0, R.shape[1]):
            R[:,i] = R[:,i] - np.mean(R[:,i])
            fr = (1+sum(R[:,i]))/(2+2*R.shape[0])
            R[:,i] = R[:,i]/math.sqrt(fr*(1-fr))
    else:
        for i in range(0, R.shape[1]):
            R[:,i] = R[:,i] - np.mean(R[:,i])

    R = R.T

    return R


if __name__ == '__main__':
    print("HEY HEY HEY")
