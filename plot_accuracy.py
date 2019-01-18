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
    if (os.path.isfile(str(data)+'_SVM_bestaccuracy.txt')):
        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_SVM_bestaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='SVM')
        plt.legend(loc=4)

        plt.figure(111)
        plt.plot(pvals,x,label='SVM')


        #########################################################################
        # Make naive assumption that if the above file exists then so do the pos/neg files
        # Positive accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_SVM_positiveaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='SVM')
        plt.legend(loc=4)

        plt.figure(1)
        plt.plot(pvals,x,label='SVM')


        # Negative accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_SVM_negativeaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='SVM')
        plt.legend(loc=4)

        plt.figure(3)
        plt.plot(pvals,x,label='SVM')

        # F1 score plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_SVM_f1score.txt'
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
        plt.plot(pvalsn,ync,label='SVM')
        plt.legend(loc=4)

        plt.figure(16)
        plt.plot(pvals,x,label='SVM')

        err_flag = 1

    if (os.path.isfile(str(data)+'_QDA_bestaccuracy.txt')):
        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_QDA_bestaccuracy.txt'
        # Times 2 because of formatting that is how many lines there are per run
        thenum = 2*len(pvals)
        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
        subprocess.call(COMMAND, shell=True)

        # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
        # also grab the last column i.e. the accuracy
        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

        x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

        plt.figure(121)
        fc = interpolate.interp1d(pvals, x, kind='quadratic')
        pvalsn = np.linspace(50,2000,10000)
        ync = fc(pvalsn)
        plt.plot(pvalsn,ync,label='QDA')
        plt.legend(loc=4)

        plt.figure(111)
        plt.plot(pvals,x,label='QDA')

        #########################################################################
        # Make naive assumption that if the above file exists then so do the pos/neg files
        # Positive accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_QDA_positiveaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='QDA')
        plt.legend(loc=4)

        plt.figure(1)
        plt.plot(pvals,x,label='QDA')


        # Negative accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_QDA_negativeaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='QDA')
        plt.legend(loc=4)

        plt.figure(3)
        plt.plot(pvals,x,label='QDA')

        # F1 score plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_QDA_f1score.txt'
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
        plt.plot(pvalsn,ync,label='QDA')
        plt.legend(loc=4)

        plt.figure(16)
        plt.plot(pvals,x,label='QDA')

        err_flag = 1

    if (os.path.isfile(str(data)+'_LDA_bestaccuracy.txt')):
        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_LDA_bestaccuracy.txt'
        # Times 2 because of formatting that is how many lines there are per run
        thenum = 2*len(pvals)
        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
        subprocess.call(COMMAND, shell=True)

        # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
        # also grab the last column i.e. the accuracy
        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

        x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

        plt.figure(121)
        fc = interpolate.interp1d(pvals, x, kind='quadratic')
        pvalsn = np.linspace(50,2000,10000)
        ync = fc(pvalsn)
        plt.plot(pvalsn,ync,label='LDA')
        plt.legend(loc=4)

        plt.figure(111)
        plt.plot(pvals,x,label='LDA')

        #########################################################################
        # Make naive assumption that if the above file exists then so do the pos/neg files
        # Positive accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_LDA_positiveaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='LDA')
        plt.legend(loc=4)

        plt.figure(1)
        plt.plot(pvals,x,label='LDA')


        # Negative accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_LDA_negativeaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='LDA')
        plt.legend(loc=4)

        plt.figure(3)
        plt.plot(pvals,x,label='LDA')

        # F1 score plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_LDA_f1score.txt'
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
        plt.plot(pvalsn,ync,label='LDA')
        plt.legend(loc=4)

        plt.figure(16)
        plt.plot(pvals,x,label='LDA')

        err_flag = 1

    if (os.path.isfile(str(data)+'_PolyReg_bestaccuracy.txt')):
        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_PolyReg_bestaccuracy.txt'
        # Times 2 because of formatting that is how many lines there are per run
        thenum = 2*len(pvals)
        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
        subprocess.call(COMMAND, shell=True)

        # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
        # also grab the last column i.e. the accuracy
        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

        x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

        plt.figure(121)
        fc = interpolate.interp1d(pvals, x, kind='quadratic')
        pvalsn = np.linspace(50,2000,10000)
        ync = fc(pvalsn)
        plt.plot(pvalsn,ync,label='PolyReg')
        plt.legend(loc=4)

        plt.figure(111)
        plt.plot(pvals,x,label='PolyReg')

        #########################################################################
        # Make naive assumption that if the above file exists then so do the pos/neg files
        # Positive accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_PolyReg_positiveaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='PolyReg')
        plt.legend(loc=4)

        plt.figure(1)
        plt.plot(pvals,x,label='PolyReg')


        # Negative accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_PolyReg_negativeaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='PolyReg')
        plt.legend(loc=4)

        plt.figure(3)
        plt.plot(pvals,x,label='PolyReg')

        # F1 score plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_PolyReg_f1score.txt'
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
        plt.plot(pvalsn,ync,label='PolyReg')
        plt.legend(loc=4)

        plt.figure(16)
        plt.plot(pvals,x,label='PolyReg')

        err_flag = 1

    if (os.path.isfile(str(data)+'_Ridge_bestaccuracy.txt')):
        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_Ridge_bestaccuracy.txt'
        # Times 2 because of formatting that is how many lines there are per run
        thenum = 2*len(pvals)
        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
        subprocess.call(COMMAND, shell=True)

        # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
        # also grab the last column i.e. the accuracy
        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

        x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

        plt.figure(121)
        fc = interpolate.interp1d(pvals, x, kind='quadratic')
        pvalsn = np.linspace(50,2000,10000)
        ync = fc(pvalsn)
        plt.plot(pvalsn,ync,label='Ridge')
        plt.legend(loc=4)

        plt.figure(111)
        plt.plot(pvals,x,label='Ridge')

        #########################################################################
        # Make naive assumption that if the above file exists then so do the pos/neg files
        # Positive accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_Ridge_positiveaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='Ridge')
        plt.legend(loc=4)

        plt.figure(1)
        plt.plot(pvals,x,label='Ridge')


        # Negative accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_Ridge_negativeaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='Ridge')
        plt.legend(loc=4)

        plt.figure(3)
        plt.plot(pvals,x,label='Ridge')

        # F1 score plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_Ridge_f1score.txt'
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
        plt.plot(pvalsn,ync,label='Ridge')
        plt.legend(loc=4)

        plt.figure(16)
        plt.plot(pvals,x,label='Ridge')

        err_flag = 1

    if (os.path.isfile(str(data)+'_Lasso_bestaccuracy.txt')):
        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_Lasso_bestaccuracy.txt'
        # Times 2 because of formatting that is how many lines there are per run
        thenum = 2*len(pvals)
        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
        subprocess.call(COMMAND, shell=True)

        # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
        # also grab the last column i.e. the accuracy
        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

        x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

        plt.figure(121)
        fc = interpolate.interp1d(pvals, x, kind='quadratic')
        pvalsn = np.linspace(50,2000,10000)
        ync = fc(pvalsn)
        plt.plot(pvalsn,ync,label='Lasso')
        plt.legend(loc=4)

        plt.figure(111)
        plt.plot(pvals,x,label='Lasso')

        #########################################################################
        # Make naive assumption that if the above file exists then so do the pos/neg files
        # Positive accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_Lasso_positiveaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='Lasso')
        plt.legend(loc=4)

        plt.figure(1)
        plt.plot(pvals,x,label='Lasso')


        # Negative accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_Lasso_negativeaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='Lasso')
        plt.legend(loc=4)

        plt.figure(3)
        plt.plot(pvals,x,label='Lasso')

        # F1 score plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_Lasso_f1score.txt'
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
        plt.plot(pvalsn,ync,label='Lasso')
        plt.legend(loc=4)

        plt.figure(16)
        plt.plot(pvals,x,label='Lasso')

        err_flag = 1

    if (os.path.isfile(str(data)+'_Elastic_bestaccuracy.txt')):
        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_Elastic_bestaccuracy.txt'
        # Times 2 because of formatting that is how many lines there are per run
        thenum = 2*len(pvals)
        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
        subprocess.call(COMMAND, shell=True)

        # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
        # also grab the last column i.e. the accuracy
        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

        x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

        plt.figure(121)
        fc = interpolate.interp1d(pvals, x, kind='quadratic')
        pvalsn = np.linspace(50,2000,10000)
        ync = fc(pvalsn)
        plt.plot(pvalsn,ync,label='Elastic')
        plt.legend(loc=4)

        plt.figure(111)
        plt.plot(pvals,x,label='Elastic')

        #########################################################################
        # Make naive assumption that if the above file exists then so do the pos/neg files
        # Positive accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_Elastic_positiveaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='Elastic')
        plt.legend(loc=4)

        plt.figure(1)
        plt.plot(pvals,x,label='Elastic')


        # Negative accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_Elastic_negativeaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='Elastic')
        plt.legend(loc=4)

        plt.figure(3)
        plt.plot(pvals,x,label='Elastic')

        # F1 score plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_Elastic_f1score.txt'
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
        plt.plot(pvalsn,ync,label='Elastic')
        plt.legend(loc=4)

        plt.figure(16)
        plt.plot(pvals,x,label='Elastic')

        err_flag = 1

    if (os.path.isfile(str(data)+'_kSVM_bestaccuracy.txt')):
        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_kSVM_bestaccuracy.txt'
        # Times 2 because of formatting that is how many lines there are per run
        thenum = 2*len(pvals)
        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
        subprocess.call(COMMAND, shell=True)

        # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
        # also grab the last column i.e. the accuracy
        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

        x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

        plt.figure(121)
        fc = interpolate.interp1d(pvals, x, kind='quadratic')
        pvalsn = np.linspace(50,2000,10000)
        ync = fc(pvalsn)
        plt.plot(pvalsn,ync,label='kSVM')
        plt.legend(loc=4)

        plt.figure(111)
        plt.plot(pvals,x,label='kSVM')

        #########################################################################
        # Make naive assumption that if the above file exists then so do the pos/neg files
        # Positive accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_kSVM_positiveaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='kSVM')
        plt.legend(loc=4)

        plt.figure(1)
        plt.plot(pvals,x,label='kSVM')


        # Negative accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_kSVM_negativeaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='kSVM')
        plt.legend(loc=4)

        plt.figure(3)
        plt.plot(pvals,x,label='kSVM')

        # F1 score plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_kSVM_f1score.txt'
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
        plt.plot(pvalsn,ync,label='kSVM')
        plt.legend(loc=4)

        plt.figure(16)
        plt.plot(pvals,x,label='kSVM')

        err_flag = 1

    if (os.path.isfile(str(data)+'_RidgeSVM_bestaccuracy.txt')):
        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_RidgeSVM_bestaccuracy.txt'
        # Times 2 because of formatting that is how many lines there are per run
        thenum = 2*len(pvals)
        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
        subprocess.call(COMMAND, shell=True)

        # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
        # also grab the last column i.e. the accuracy
        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

        x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

        plt.figure(121)
        fc = interpolate.interp1d(pvals, x, kind='quadratic')
        pvalsn = np.linspace(50,2000,10000)
        ync = fc(pvalsn)
        plt.plot(pvalsn,ync,label='RidgeSVM')
        plt.legend(loc=4)

        plt.figure(111)
        plt.plot(pvals,x,label='RidgeSVM')

        #########################################################################
        # Make naive assumption that if the above file exists then so do the pos/neg files
        # Positive accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_RidgeSVM_positiveaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='RidgeSVM')
        plt.legend(loc=4)

        plt.figure(1)
        plt.plot(pvals,x,label='RidgeSVM')


        # Negative accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_RidgeSVM_negativeaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='RidgeSVM')
        plt.legend(loc=4)

        plt.figure(3)
        plt.plot(pvals,x,label='RidgeSVM')

        # F1 score plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_RidgeSVM_f1score.txt'
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
        plt.plot(pvalsn,ync,label='RidgeSVM')
        plt.legend(loc=4)

        plt.figure(16)
        plt.plot(pvals,x,label='RidgeSVM')

        err_flag = 1

    if (os.path.isfile(str(data)+'_RFESVM_bestaccuracy.txt')):
        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_RFESVM_bestaccuracy.txt'
        # Times 2 because of formatting that is how many lines there are per run
        thenum = 2*len(pvals)
        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
        subprocess.call(COMMAND, shell=True)

        # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
        # also grab the last column i.e. the accuracy
        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

        x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

        plt.figure(121)
        fc = interpolate.interp1d(pvals, x, kind='quadratic')
        pvalsn = np.linspace(50,2000,10000)
        ync = fc(pvalsn)
        plt.plot(pvalsn,ync,label='RFESVM')
        plt.legend(loc=4)

        plt.figure(111)
        plt.plot(pvals,x,label='RFESVM')

        #########################################################################
        # Make naive assumption that if the above file exists then so do the pos/neg files
        # Positive accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_RFESVM_positiveaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='RFESVM')
        plt.legend(loc=4)

        plt.figure(1)
        plt.plot(pvals,x,label='RFESVM')


        # Negative accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_RFESVM_negativeaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='RFESVM')
        plt.legend(loc=4)

        plt.figure(3)
        plt.plot(pvals,x,label='RFESVM')

        # F1 score plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_RFESVM_f1score.txt'
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
        plt.plot(pvalsn,ync,label='RFESVM')
        plt.legend(loc=4)

        plt.figure(16)
        plt.plot(pvals,x,label='RFESVM')

        err_flag = 1

    if (os.path.isfile(str(data)+'_RandomForest_bestaccuracy.txt')):
        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_RandomForest_bestaccuracy.txt'
        # Times 2 because of formatting that is how many lines there are per run
        thenum = 2*len(pvals)
        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
        subprocess.call(COMMAND, shell=True)

        # Awk command to grab every even line (odd lines are just formatting even lines have the accuracies)
        # also grab the last column i.e. the accuracy
        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"

        x = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

        plt.figure(121)
        fc = interpolate.interp1d(pvals, x, kind='quadratic')
        pvalsn = np.linspace(50,2000,10000)
        ync = fc(pvalsn)
        plt.plot(pvalsn,ync,label='RandomForest')
        plt.legend(loc=4)

        plt.figure(111)
        plt.plot(pvals,x,label='RandomForest')

        #########################################################################
        # Make naive assumption that if the above file exists then so do the pos/neg files
        # Positive accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_RandomForest_positiveaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='RandomForest')
        plt.legend(loc=4)

        plt.figure(1)
        plt.plot(pvals,x,label='RandomForest')


        # Negative accuracy plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_RandomForest_negativeaccuracy.txt'
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
        plt.plot(pvalsn,ync,label='RandomForest')
        plt.legend(loc=4)

        plt.figure(3)
        plt.plot(pvals,x,label='RandomForest')

        # F1 score plots

        # Awk command to grab the last so many lines of the file so that I plot the
        # most recent run of accuracies
        thefile = str(data)+'_RandomForest_f1score.txt'
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
        plt.plot(pvalsn,ync,label='RandomForest')
        plt.legend(loc=4)

        plt.figure(16)
        plt.plot(pvals,x,label='RandomForest')

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
