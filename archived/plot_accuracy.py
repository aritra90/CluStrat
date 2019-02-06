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

# Class Options = [LDA,QDA,PolyReg,Ridge,Lasso,Elastic,SVM,kSVM,RidgeSVM,RFESVM,RandomForest]

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dt", "--data", dest='data', action='store', help="Put the dataset of accuracy file.",
                        metavar="DATA")
    args = parser.parse_args()
    return args.data

# fc = interpolate.interp1d(pvals, x, kind='quadratic')
# pvalsn = np.linspace(50,2000,10000)
# ync = fc(pvalsn)
# plt.plot(pvalsn,ync,label='SVM')

if __name__ == "__main__":
    err_flag = 0
    cvFile = open("crossvalType.txt","r")
    cv_flag = str(cvFile.read()).strip()
    print(cv_flag)
    data = parse_arguments()
    pvals = [50, 100, 200, 300, 500, 1000, 2000]
    alg_options = ['LDA','QDA','PolyReg','Ridge','Lasso','Elastic','SVM','kSVM','RidgeSVM','RFESVM','RandomForest']
    for alg in alg_options:
        print(alg)
        if (os.path.isfile(str(data)+'_'+str(alg)+'_bestaccuracy.txt')):
            # Awk command to grab the last so many lines of the file so that I plot the
            # most recent run of accuracies
            thefile = str(data)+'_'+str(alg)+'_bestaccuracy.txt'
            # Times 2 because of formatting that is how many lines there are per run
            thenum = 2*len(pvals)
            COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
            subprocess.call(COMMAND, shell=True)

            # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
            # also grab the last column i.e. the accuracy
            COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

            # Get the accuracy using the command
            x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

            plt.figure(121)
            fc = interpolate.interp1d(pvals, x, kind='quadratic')
            pvalsn = np.linspace(50,2000,10000)
            ync = fc(pvalsn)
            plt.plot(pvalsn,ync,label=str(alg))
            plt.legend(loc=4)

            plt.figure(111)
            plt.plot(pvals,x,label=str(alg))


            #########################################################################
            # Make naive assumption that if the above file exists then so do the pos/neg files
            # Positive accuracy plots

            # Awk command to grab the last so many lines of the file so that I plot the
            # most recent run of accuracies
            thefile = str(data)+'_'+str(alg)+'_positiveaccuracy.txt'
            # Times 2 because of formatting that is how many lines there are per run
            thenum = 2*len(pvals)
            COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
            subprocess.call(COMMAND, shell=True)

            # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
            # also grab the last column i.e. the accuracy
            COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

            x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())
            plt.figure(0)
            fc = interpolate.interp1d(pvals, x, kind='quadratic')
            pvalsn = np.linspace(50,2000,10000)
            ync = fc(pvalsn)
            plt.plot(pvalsn,ync,label=str(alg))
            plt.legend(loc=4)

            plt.figure(1)
            plt.plot(pvals,x,label=str(alg))


            # Negative accuracy plots

            # Awk command to grab the last so many lines of the file so that I plot the
            # most recent run of accuracies
            thefile = str(data)+'_'+str(alg)+'_negativeaccuracy.txt'
            # Times 2 because of formatting that is how many lines there are per run
            thenum = 2*len(pvals)
            COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
            subprocess.call(COMMAND, shell=True)

            # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
            # also grab the last column i.e. the accuracy
            COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

            x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())
            plt.figure(2)
            fc = interpolate.interp1d(pvals, x, kind='quadratic')
            pvalsn = np.linspace(50,2000,10000)
            ync = fc(pvalsn)
            plt.plot(pvalsn,ync,label=str(alg))
            plt.legend(loc=4)

            plt.figure(3)
            plt.plot(pvals,x,label=str(alg))

            # F1 score plots

            # Awk command to grab the last so many lines of the file so that I plot the
            # most recent run of accuracies
            thefile = str(data)+'_'+str(alg)+'_f1score.txt'
            # Times 2 because of formatting that is how many lines there are per run
            thenum = 2*len(pvals)
            COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
            subprocess.call(COMMAND, shell=True)

            # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
            # also grab the last column i.e. the accuracy
            COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

            x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())
            plt.figure(15)
            fc = interpolate.interp1d(pvals, x, kind='quadratic')
            pvalsn = np.linspace(50,2000,10000)
            ync = fc(pvalsn)
            plt.plot(pvalsn,ync,label=str(alg))
            plt.legend(loc=4)

            plt.figure(16)
            plt.plot(pvals,x,label=str(alg))

            err_flag = 1

    if err_flag == 0:
        # print("YES!")
        print("Could not find an accuracy file corresponding to dataset prefix: "+str(data))
        print("Make sure file is in the following format: prefix_classifier_bestaccuracy.txt")
        sys.exit(1)

    ################### Total Accuracy Plots ###################
    plt.figure(111)
    plt.xticks(pvals,rotation=-45)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Number of Top SNPs Extracted")
    plt.ylabel("Accuracy")
    plt.title('Classifier Accuracy')
    plt.grid()
    plt.legend(loc=4)
    plt.tight_layout()
    plt.savefig(str(data)+'_classifier_accuracy_'+str(cv_flag)+'.png')

    plt.figure(121)
    plt.xticks(pvals,rotation=-45)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Number of Top SNPs Extracted")
    plt.ylabel("Accuracy")
    plt.title('Classifier Accuracy')
    plt.grid()
    plt.tight_layout()
    # plt.legend(loc=4)
    plt.savefig(str(data)+'_classifier_accuracy_smooth_'+str(cv_flag)+'.png')

    ################### Case Accuracy Plots (true positives) ###################
    plt.figure(1)
    plt.xticks(pvals,rotation=-45)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Number of Top SNPs Extracted")
    plt.ylabel("Accuracy")
    plt.title('Case Classifier Accuracy')
    plt.grid()
    plt.tight_layout()
    plt.legend(loc=4)
    plt.savefig(str(data)+'_classifier_positiveaccuracy_'+str(cv_flag)+'.png')

    plt.figure(0)
    plt.xticks(pvals,rotation=-45)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Number of Top SNPs Extracted")
    plt.ylabel("Accuracy")
    plt.title('Case Classifier Accuracy')
    plt.grid()
    plt.tight_layout()
    # plt.legend(loc=4)
    plt.savefig(str(data)+'_classifier_positiveaccuracy_smooth_'+str(cv_flag)+'.png')


    ################### Control Accuracy Plots (true negatives) ###################
    plt.figure(3)
    plt.xticks(pvals,rotation=-45)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Number of Top SNPs Extracted")
    plt.ylabel("Accuracy")
    plt.title('Control Classifier Accuracy')
    plt.grid()
    plt.tight_layout()
    plt.legend(loc=4)
    plt.savefig(str(data)+'_classifier_negativeaccuracy_'+str(cv_flag)+'.png')

    plt.figure(2)
    plt.xticks(pvals,rotation=-45)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Number of Top SNPs Extracted")
    plt.ylabel("Accuracy")
    plt.title('Control Classifier Accuracy')
    plt.grid()
    plt.tight_layout()
    # plt.legend(loc=4)
    plt.savefig(str(data)+'_classifier_negativeaccuracy_smooth_'+str(cv_flag)+'.png')


    ################### F1 Score Plots ###################
    plt.figure(16)
    plt.xticks(pvals,rotation=-45)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Number of Top SNPs Extracted")
    plt.ylabel("F1 Score")
    plt.title('Classifier F1 Score')
    plt.grid()
    plt.tight_layout()
    plt.legend(loc=4)
    plt.savefig(str(data)+'_classifier_f1score_'+str(cv_flag)+'.png')

    plt.figure(15)
    plt.xticks(pvals,rotation=-45)
    plt.ylim(ymin=0,ymax=1.2)
    plt.xlabel("Number of Top SNPs Extracted")
    plt.ylabel("F1 Score")
    plt.title('Classifier F1 Score')
    plt.grid()
    plt.tight_layout()
    # plt.legend(loc=4)
    plt.savefig(str(data)+'_classifier_f1score_smooth_'+str(cv_flag)+'.png')
    plt.close('all')


###############################################################################################################################
