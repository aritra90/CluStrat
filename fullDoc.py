from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, LongTable,\
    Plot, Figure, SubFigure, NoEscape, Matrix, Alignat, NewPage, TextBlock, MediumText, HugeText,\
    SmallText
from pylatex.utils import italic
from pylatex.utils import bold
import pylatex
import os

import argparse
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.interpolate import spline
from scipy import interpolate
import numpy as np

# Class Options = [LDA,QDA,PolyReg,Ridge,Lasso,Elastic,SVM,kSVM,RidgeSVM,RFESVM,RandomForest]

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dt", "--data", dest='data', action='store', help="Put the dataset of accuracy file.",
                        metavar="DATA")
    args = parser.parse_args()
    return args.data

if __name__ == '__main__':
    data = parse_arguments()

    if (str(data) == 'PRK'):
        str_dataset = "PKHCGSRRS"
    elif (str(data) == 'SCZ'):
        str_dataset = "SZHCGSRRS"
    else:
        sys.exit(1)

    # image_filename = os.path.join(os.path.dirname(__file__), 'kitten.jpg')
    geometry_options = {"tmargin": "1cm", "lmargin": "1cm"}
    doc = Document('MLGWAS_Results_'+str_dataset,geometry_options=geometry_options)

    doc.append(pylatex.Command('fontsize', arguments = ['10', '12']))
    doc.append(pylatex.Command('selectfont'))
    doc.change_document_style("empty")

    with doc.create(Section('Dataset Description')):
        # doc.append('Some regular text and some')
        # doc.append(italic('italic text. '))
        # doc.append('\nAlso some crazy characters: $&#{}')
        doc.append('The datasets used for these experiments are the Parkinsons dataset (PKHCGSRRS) and the Schizophreina '+
        'dataset (SZHCGSRRS). For the Schizophreina dataset, there are 5,921 individuals. 5,440 are control and 481 are cases. '+
        '40 individuals (20 cases and 20 controls) are used for testing and the rest are used for training (5,881 individuals, 461 cases and 5,420 controls). '+
        'For the Parkinsons dataset, there are 4,706 individuals. 2,837 are control and 1,869 are cases. '+
        '40 individuals (20 cases and 20 controls) are used for testing and the rest are used for training (4,666 individuals, 1,849 cases and 2,817 controls). '+
        'We will always take this approach to split train and test sets unless we use balanced subsampling for leave-one-out cross-validation '+
        'This way we can train on a small class-balanced set and this make leave-one-out cross-validation much more efficient. '+
        '\n\nEach dataset is ran through the GWAS pipeline to find the top ~2K SNPs that we believe are '+
        'associated with the corresponding disease. Different p-value thresholds are used to get ~2K top SNPs '+
        '(0.02 for Schizophreina dataset and 0.003 for Parkinsons dataset). The other default parameters for the GWAS pipeline '+
        'are used. ')

    # doc.append(NewPage())

    err_flag = 0
    cvFile = open("crossvalType.txt","r")
    cv_flag = str(cvFile.read()).strip()
    # print(cv_flag)

    pvals = [50, 100, 200, 300, 500, 1000, 2000]
    alg_options = ['LDA','QDA','PolyReg','Ridge','Lasso','Elastic','SVM','SVR','kSVM','RidgeSVM','RidgeCV','RFESVM','RandomForest']

    with doc.create(Section('Results')):
        if (str(data) == 'PRK'):
            doc.append('Dataset:'+SmallText(bold(" PKHCGSRRS"))+'\nCross-validation type:'+SmallText(bold(' ' +str(cv_flag))))
        elif (str(data) == 'SCZ'):
            doc.append('Dataset:'+SmallText(bold(" SZHCGSRRS"))+'\nCross-validation:'+SmallText(bold(' '+str(cv_flag))))
        else:
            sys.exit(1)
        with doc.create(Subsection('Accuracy Tables')):
            with doc.create(LongTable('c|c|c|c|c|c|c|c')) as table:
                table.add_row(('Algorithm','# of SNPs','Training Acc','Confusion Matrix','Testing Acc','F1 Score','Case Acc','Control Acc'))
                table.add_hline()
                table.add_hline()
                for alg in alg_options:
                    if (os.path.isfile(str(data)+'_'+str(alg)+'_bestaccuracy.txt')):

                        thefile = str(data)+'_'+str(alg)+'_trainingaccuracy.txt'
                        thenum = 2*len(pvals)
                        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
                        subprocess.call(COMMAND, shell=True)
                        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"
                        a = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

                        plt.figure(18)
                        fc = interpolate.interp1d(pvals, a, kind='quadratic')
                        pvalsn = np.linspace(50,2000,10000)
                        ync = fc(pvalsn)
                        plt.plot(pvalsn,ync,label=str(alg))
                        plt.legend(loc=4)
                        # debugging
                        # plt.savefig('dummy1.png')
                        plt.figure(17)
                        plt.plot(pvals,a,label=str(alg))
                        plt.legend(loc=4)
                        # debugging
                        # plt.savefig('dummy2.png')




                        thefile = str(data)+'_'+str(alg)+'_bestaccuracy.txt'
                        thenum = 2*len(pvals)
                        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
                        subprocess.call(COMMAND, shell=True)
                        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"
                        b = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

                        plt.figure(121)
                        fc = interpolate.interp1d(pvals, b, kind='quadratic')
                        pvalsn = np.linspace(50,2000,10000)
                        ync = fc(pvalsn)
                        plt.plot(pvalsn,ync,label=str(alg))
                        plt.legend(loc=4)
                        plt.figure(111)
                        plt.plot(pvals,b,label=str(alg))




                        thefile = str(data)+'_'+str(alg)+'_caseaccuracy.txt'
                        thenum = 2*len(pvals)
                        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
                        subprocess.call(COMMAND, shell=True)
                        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"
                        c = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

                        plt.figure(0)
                        fc = interpolate.interp1d(pvals, c, kind='quadratic')
                        pvalsn = np.linspace(50,2000,10000)
                        ync = fc(pvalsn)
                        plt.plot(pvalsn,ync,label=str(alg))
                        plt.legend(loc=4)
                        plt.figure(1)
                        plt.plot(pvals,c,label=str(alg))




                        thefile = str(data)+'_'+str(alg)+'_controlaccuracy.txt'
                        thenum = 2*len(pvals)
                        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
                        subprocess.call(COMMAND, shell=True)
                        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"
                        d = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

                        plt.figure(2)
                        fc = interpolate.interp1d(pvals, d, kind='quadratic')
                        pvalsn = np.linspace(50,2000,10000)
                        ync = fc(pvalsn)
                        plt.plot(pvalsn,ync,label=str(alg))
                        plt.legend(loc=4)
                        plt.figure(3)
                        plt.plot(pvals,d,label=str(alg))




                        thefile = str(data)+'_'+str(alg)+'_f1score.txt'
                        thenum = 2*len(pvals)
                        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
                        subprocess.call(COMMAND, shell=True)
                        COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"
                        e = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

                        plt.figure(15)
                        fc = interpolate.interp1d(pvals, e, kind='quadratic')
                        pvalsn = np.linspace(50,2000,10000)
                        ync = fc(pvalsn)
                        plt.plot(pvalsn,ync,label=str(alg))
                        plt.legend(loc=4)
                        plt.figure(16)
                        plt.plot(pvals,e,label=str(alg))


                        err_flag = 1

                        thefile = str(data)+'_'+str(alg)+'_confusionmatrices.txt'
                        thenum = 4*len(pvals)
                        COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
                        subprocess.call(COMMAND, shell=True)
                        COMMAND = "awk 'NR%4==0||(NR+1)%4==0' tempStats.txt > tempCnfmat.txt"
                        subprocess.call(COMMAND, shell=True)
                        thefile = open( "tempCnfmat.txt", "r" )
                        f = []
                        for line in thefile:
                            f.append(line)

                        max_case_accuracy = max(c)
                        max_case_index = c.index(max_case_accuracy)

                        for ind in range(len(pvals)):
                            table.add_empty_row()
                            if ind == 3:
                                table.add_row((str(alg), pvals[ind], a[ind], f[2*ind], b[ind], e[ind], c[ind], d[ind]))
                                table.add_row(('', '', '', f[2*ind+1], '', '', '', ''))
                            elif ind == max_case_index:
                                table.add_row(('', pvals[ind], a[ind], f[2*ind], b[ind], e[ind], c[ind], d[ind]),
                                    mapper=bold, color="yellow")
                                table.add_row(('', '', '', f[2*ind+1], '', '', '', ''),
                                    mapper=bold, color="yellow")
                            else:
                                table.add_row(('', pvals[ind], a[ind], f[2*ind], b[ind], e[ind], c[ind], d[ind]))
                                table.add_row(('', '', '', f[2*ind+1], '', '', '', ''))

                            table.add_empty_row()
                            if ind != len(pvals)-1:
                                table.add_hline(2,3)
                                table.add_hline(4,5)
                                table.add_hline(6,7)
                                table.add_hline(8)
                            else:
                                table.add_hline()
                                table.add_hline()

                if err_flag == 0:
                    # print("YES!")
                    print("Could not find an accuracy file corresponding to dataset prefix: "+str(data))
                    print("Make sure file is in the following format: prefix_classifier_bestaccuracy.txt")
                    sys.exit(1)
        doc.append(NewPage())
        with doc.create(Subsection('Accuracy Plots')):

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
            plt.savefig(str(data)+'_classifier_caseaccuracy_'+str(cv_flag)+'.png')

            plt.figure(0)
            plt.xticks(pvals,rotation=-45)
            plt.ylim(ymin=0,ymax=1.2)
            plt.xlabel("Number of Top SNPs Extracted")
            plt.ylabel("Accuracy")
            plt.title('Case Classifier Accuracy')
            plt.grid()
            plt.tight_layout()
            # plt.legend(loc=4)
            plt.savefig(str(data)+'_classifier_caseaccuracy_smooth_'+str(cv_flag)+'.png')


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
            plt.savefig(str(data)+'_classifier_controlaccuracy_'+str(cv_flag)+'.png')

            plt.figure(2)
            plt.xticks(pvals,rotation=-45)
            plt.ylim(ymin=0,ymax=1.2)
            plt.xlabel("Number of Top SNPs Extracted")
            plt.ylabel("Accuracy")
            plt.title('Control Classifier Accuracy')
            plt.grid()
            plt.tight_layout()
            # plt.legend(loc=4)
            plt.savefig(str(data)+'_classifier_controlaccuracy_smooth_'+str(cv_flag)+'.png')


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


            ################### Training Accuracy Plots ###################
            plt.figure(17)
            plt.xticks(pvals,rotation=-45)
            plt.ylim(ymin=0,ymax=1.2)
            plt.xlabel("Number of Top SNPs Extracted")
            plt.ylabel("Training Accuracy")
            plt.title('Classifier Accuracy')
            plt.grid()
            # plt.legend(loc=4)
            plt.tight_layout()
            plt.savefig(str(data)+'_classifier_trainaccuracy_'+str(cv_flag)+'.png')

            plt.figure(18)
            plt.xticks(pvals,rotation=-45)
            plt.ylim(ymin=0,ymax=1.2)
            plt.xlabel("Number of Top SNPs Extracted")
            plt.ylabel("Training Accuracy")
            plt.title('Classifier Accuracy')
            plt.grid()
            plt.tight_layout()
            # plt.legend(loc=4)
            plt.savefig(str(data)+'_classifier_trainaccuracy_smooth_'+str(cv_flag)+'.png')
            plt.close('all')


            with doc.create(Figure(position='h!')) as accuracies:
                # with doc.create(SubFigure(position='b',width=NoEscape(r'0.75\linewidth'))) as left_plot:
                #     left_plot.add_image(str(data)+'_classifier_accuracy_'+str(cv_flag)+'.png',
                #                           width=NoEscape(r'\linewidth'))
                #     left_plot.add_caption('Kitten on the left')
                # with doc.create(SubFigure(position='b',width=NoEscape(r'0.75\linewidth'))) as right_plot:
                #     right_plot.add_image(str(data)+'_classifier_accuracy_smooth_'+str(cv_flag)+'.png',
                #                            width=NoEscape(r'\linewidth'))
                #     right_plot.add_caption('Kitten on the right')
                accuracies.add_image(str(data)+'_classifier_accuracy_smooth_'+str(cv_flag)+'.png',
                                       width=NoEscape(r'0.7\linewidth'))
                accuracies.add_caption('Test accuracy')

            with doc.create(Figure(position='h!')) as accuracies:
                accuracies.add_image(str(data)+'_classifier_caseaccuracy_smooth_'+str(cv_flag)+'.png',
                                       width=NoEscape(r'0.7\linewidth'))
                accuracies.add_caption('Case accuracy')

            with doc.create(Figure(position='h!')) as accuracies:
                accuracies.add_image(str(data)+'_classifier_controlaccuracy_smooth_'+str(cv_flag)+'.png',
                                       width=NoEscape(r'0.7\linewidth'))
                accuracies.add_caption('Control accuracy')

            with doc.create(Figure(position='h!')) as accuracies:
                accuracies.add_image(str(data)+'_classifier_trainaccuracy_smooth_'+str(cv_flag)+'.png',
                                       width=NoEscape(r'0.7\linewidth'))
                accuracies.add_caption('Train Accuracy')

            with doc.create(Figure(position='h!')) as accuracies:
                accuracies.add_image(str(data)+'_classifier_f1score_smooth_'+str(cv_flag)+'.png',
                                       width=NoEscape(r'0.7\linewidth'))
                accuracies.add_caption('F1 Score')

    # VERY WEIRD POSITIONING AND ODD BEHAVIOR WITH FIGURES ABOVE
    # doc.append(NewPage())
    # with doc.create(Section('Discussion')):
    #     doc.append('INSERT DISCUSSION RIGHT HERE...')

    doc.generate_pdf(clean_tex=True)
    doc.generate_tex()
