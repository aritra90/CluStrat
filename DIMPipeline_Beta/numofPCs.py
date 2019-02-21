################################
#Author: Aritra Bose. 11/08/2018
################################

import numpy as np
import argparse
import os
import subprocess
import sys

def determine_k(S):
    #Determine K by observing the gaps between the proportion of variance
    #explained by the singular values.
    #S is the singular value diagonal matrix.
    thres = 0.001
    #We can play with the threshold
    #print(S)
    S = np.diag(S)
    #print(S)
    pve = S/np.sum(S)
    #print(pve)
    # plot(1:numel(pve), pve,'r-', 'linewidth',8)
    d_pve = abs(np.diff(pve))
    #print(d_pve)
    [pve_thres,] = np.where(d_pve < thres)
    #this is a corner case when none of the differences are below the
    #thresholdo

    if (len(pve_thres) == 0):
        #print("HEY")
        pve_thres = np.argsort(d_pve)

    pve_thres = pve_thres + 1
    #print(pve_thres)
    d_pve_thres = abs(np.diff(pve_thres))
    #print(d_pve_thres)
    [det_k,] = np.where(d_pve_thres > 1)
    #print(det_k)
    #another corner case when there is convergence in pve, we select the
    #first index
    if (len(det_k) == 0):
        tmpk = 1
    else:
        tmpk = det_k[-1]+1

    #print(tmpk)

    k = pve_thres[tmpk-1]

    return k


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-S", "--singval", dest='singvals_file', action='store', help="Put the file to the matrix of singular values.",
                        metavar="eigenvalFile")
    args = parser.parse_args()
    return args.singvals_file


if __name__ == '__main__':

    # pass in the pca.covar or pca.evec file
    singvals_file = parse_arguments()
    COMMAND = "awk 'NR==1 {for (i=2; i<=NF; i++) print $i}' "+str(singvals_file)
    x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())
    singvals = np.diag(x)

    # print(singvals)
    # matrix-fy the singular values
    # below should get k = 3
    # singvals = np.array([[90.5728, 0, 0,0],[0,23.0567,0,0],[0,0,6.3494,0],[0,0,0,1.2777]])
    # below should get k = 1
    # singvals = np.array([[2.4605,0,0],[0,1.6996,0],[0,0,0.2391]])

    k = determine_k(singvals)
    print(k)
