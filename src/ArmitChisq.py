import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random


def chisq(snp, status, fact):

    # 'corr' returns a correlation matrix so the pval is on the
    # off-diagonal of 'corr'
    corr = np.corrcoef(snp,status);
    pval = fact*(corr[0,1]**2);

    return pval

if __name__ == '__main__':
    print("HEY HEY HEY")
