import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random


# function [adjR,adjstat] = EigStrat(R, U, status)
#     status = status - mean(status);
#     proj = U*U'*R;
#     adjR = R - proj;
#     adjstat = status - (U*U'*status);
# end

def project(R, U, status):

    status = status - np.mean(status);
    proj = np.matmul(np.matmul(U,U.T),R)
    adjR = R - proj;
    adjstat = status - np.matmul(np.matmul(U,U.T),status)

    return adjR, adjstat

if __name__ == '__main__':
    print("HEY HEY HEY")
