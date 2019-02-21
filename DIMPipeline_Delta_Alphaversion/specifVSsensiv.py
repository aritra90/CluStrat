import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.interpolate import spline
import os
import argparse
import subprocess
from scipy import interpolate
import numpy as np
import sys

# def parse_arguments():
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-dt", "--data", dest='data', action='store', help="Put the dataset of accuracy file.",
#                         metavar="DATA")
#     args = parser.parse_args()
#     return args.data


if __name__ == "__main__":
    # err_flag = 0
    # cvFile = open("crossvalType.txt","r")
    # cv_flag = str(cvFile.read()).strip()
    # print(cv_flag)
    # data = parse_arguments()
    pvals = [100, 500, 1000]
    ratios = ["1:1", "1:1.25", "1:1.5", "1:2", "1:2.25", "1:2.5", "1:2.75", "1:3"]
    xvals = [1,2,3,4,5,6,7,8]
    # alg_options = ['LDA','QDA','PolyReg','Ridge','Lasso','Elastic','SVM','kSVM','RidgeSVM','RFESVM','RandomForest']


    thefile = 'firstSpecificity.txt'
    thenum = 2*len(ratios)
    COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > firstStats.txt"
    subprocess.call(COMMAND, shell=True)
    COMMAND = "awk 'NR%2==0 {print $NF}' firstStats.txt"
    x1 = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

    thefile = 'firstSensitivity.txt'
    thenum = 2*len(ratios)
    COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > firstStats.txt"
    subprocess.call(COMMAND, shell=True)
    COMMAND = "awk 'NR%2==0 {print $NF}' firstStats.txt"
    x2 = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

    plt.figure(1)
    plt.plot(xvals,x1)

    plt.figure(4)
    plt.plot(xvals,x2)
#------------------------------------------------------------------------------------------------------------------------------
    thefile = 'secondSpecificity.txt'
    thenum = 2*len(ratios)
    COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > secondStats.txt"
    subprocess.call(COMMAND, shell=True)
    COMMAND = "awk 'NR%2==0 {print $NF}' secondStats.txt"
    y1 = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

    thefile = 'secondSensitivity.txt'
    thenum = 2*len(ratios)
    COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > secondStats.txt"
    subprocess.call(COMMAND, shell=True)
    COMMAND = "awk 'NR%2==0 {print $NF}' secondStats.txt"
    y2 = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

    plt.figure(2)
    plt.plot(xvals,y1)

    plt.figure(5)
    plt.plot(xvals,y2)
#------------------------------------------------------------------------------------------------------------------------------
    thefile = 'thirdSpecificity.txt'
    thenum = 2*len(ratios)
    COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > thirdStats.txt"
    subprocess.call(COMMAND, shell=True)
    COMMAND = "awk 'NR%2==0 {print $NF}' thirdStats.txt"
    z1 = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

    thefile = 'thirdSensitivity.txt'
    thenum = 2*len(ratios)
    COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > thirdStats.txt"
    subprocess.call(COMMAND, shell=True)
    COMMAND = "awk 'NR%2==0 {print $NF}' thirdStats.txt"
    z2 = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

    plt.figure(3)
    plt.plot(xvals,z1)

    plt.figure(6)
    plt.plot(xvals,z2)

#------------------------------------------------------------------------------------------------------------------------------

    plt.figure(1)
    plt.xticks(xvals,ratios)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Case:Control Ratios")
    plt.ylabel("Specificity")
    plt.title('Case:Control Ratios vs Specificity with 100 SNPs Extracted')
    plt.grid()
    plt.tight_layout()
    plt.savefig('RatiosVsSpecificity_100.png')

    plt.figure(2)
    plt.xticks(xvals,ratios)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Case:Control Ratios")
    plt.ylabel("Specificity")
    plt.title('Case:Control Ratios vs Specificity with 500 SNPs Extracted')
    plt.grid()
    plt.tight_layout()
    plt.savefig('RatiosVsSpecificity_500.png')

    plt.figure(3)
    plt.xticks(xvals,ratios)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Case:Control Ratios")
    plt.ylabel("Specificity")
    plt.title('Case:Control Ratios vs Specificity with 1000 SNPs Extracted')
    plt.grid()
    plt.tight_layout()
    plt.savefig('RatiosVsSpecificity_1000.png')

#------------------------------------------------------------------------------------------------------------------------------

    plt.figure(4)
    plt.xticks(xvals,ratios)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Case:Control Ratios")
    plt.ylabel("Sensitivity")
    plt.title('Case:Control Ratios vs Sensitivity with 100 SNPs Extracted')
    plt.grid()
    plt.tight_layout()
    plt.savefig('RatiosVsSensitivity_100.png')

    plt.figure(5)
    plt.xticks(xvals,ratios)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Case:Control Ratios")
    plt.ylabel("Sensitivity")
    plt.title('Case:Control Ratios vs Sensitivity with 500 SNPs Extracted')
    plt.grid()
    plt.tight_layout()
    plt.savefig('RatiosVsSensitivity_500.png')

    plt.figure(6)
    plt.xticks(xvals,ratios)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Case:Control Ratios")
    plt.ylabel("Sensitivity")
    plt.title('Case:Control Ratios vs Sensitivity with 1000 SNPs Extracted')
    plt.grid()
    plt.tight_layout()
    plt.savefig('RatiosVsSensitivity_1000.png')


###############################################################################################################################
