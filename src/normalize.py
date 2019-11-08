import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random

# Transposing downstream will be an issue here
def norm(X, flag):
    R = X.astype(float)
    fr = np.zeros(R.shape[1])
    # print(R)
    #if flag == 1:
    for i in range(0, R.shape[1]):
        fr[i] = (1+sum(R[:,i]))/(2+2*R.shape[0])
        if flag == 1:    
            R[:,i] = R[:,i] - np.mean(R[:,i])
            R[:,i] = R[:,i]/math.sqrt(fr[i]*(1-fr[i]))
        else:
            R[:,i] = R[:,i] - np.mean(R[:,i])

    return R, fr


if __name__ == '__main__':
    print("HEY HEY HEY")
