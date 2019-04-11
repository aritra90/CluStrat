from sklearn.metrics import accuracy_score, f1_score, confusion_matrix
from sklearn.metrics import precision_recall_fscore_support as prfs

from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

import pickle

# def plot_roc_curve(fpr, tpr):
#     plt.plot(fpr, tpr, color='orange', label='ROC')
#     plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
#     plt.xlabel('False Positive Rate')
#     plt.ylabel('True Positive Rate')
#     plt.title('Receiver Operating Characteristic (ROC) Curve')
#     plt.legend()
#     plt.show()

def generate_output(test_pheno, ytestpred, classifiertype, prefix, value, best_training_accuracy, best_params):
    # metrics = ['accuracy', 'f1_micro', 'f1_macro', 'f1_weighted']
    # # headers for output files
    dataset_name = str(prefix)
    file_prefix = dataset_name+'_'+classifiertype
    # # time_string = datetime.datetime.now()
    # rn = random.randint(1,200000)
    #
    # outdir = file_prefix+'_'+str(value)+'_Results_'+str(rn)
    # os.mkdir(outdir)
    #
    # doc = rp.generate_pdf_with_name(file_prefix,rn)
    #
    #
    # # build the confusion matrix
    stats_var = prfs(test_pheno, ytestpred, average='macro')
    cnfmat = confusion_matrix(test_pheno, ytestpred)
    # rp.fill_document(doc, [rp_model], metrics, cnfmat, dataset_name, num_indivs, n_snps, G_train.shape, G_test.shape,5)

    print('\n****RESULT STATISTICS****')
    print('\n===============================================\n')
    print("\n Precision, Recall and F-score: Macro\n")
    print("\n Precision, Recall and F-score: Micro\n")
    print(prfs(test_pheno, ytestpred, average='micro'))
    print("\n Precision, Recall and F-score: Weighted\n")
    print(prfs(test_pheno, ytestpred, average='weighted'))
    print('\n===============================================\n')
    print('\nConfusion Matrix:')
    print('\n===============================================\n')
    print(cnfmat)
    print('\n===============================================\n')
    # print(classification_report(test_pheno, ytestpred))

    # # Move the files to their own directory
    # absolute_prefix = file_prefix+'_'+str(rn)
    #
    # shutil.move(absolute_prefix+'.tex', outdir+'/'+absolute_prefix+'.tex')
    # shutil.move(absolute_prefix+'.pdf', outdir+'/'+absolute_prefix+'.pdf')

    best_accuracy = (float(cnfmat[0][0]+cnfmat[1][1]))/(cnfmat[0][0]+cnfmat[0][1]+cnfmat[1][0]+cnfmat[1][1])
    best_cnfmat = cnfmat
    # Also store the corresponding case(postive) and control(negative) accuracies
    # TN/(TN+FP) i.e. specificity
    control_accuracy = (float(cnfmat[0][0]))/(cnfmat[0][0]+cnfmat[0][1])
    # TP/(TP+FN) i.e. sensitivity
    case_accuracy = (float(cnfmat[1][1]))/(cnfmat[1][0]+cnfmat[1][1])
    # Store corresponding f1 score
    f1score = round(stats_var[1],3)

    # ROC Curve stuff
    fpr, tpr, thresholds = roc_curve(test_pheno, ytestpred)

    with open(file_prefix+'_ROC_stats.txt', 'wb') as f:
        pickle.dump([fpr, tpr, thresholds], f)

    f = open(file_prefix+'_bestparameters.txt','a+')
    f.write('#--------------------------------------------------------------------#\n')
    f.write('Best parameters for top '+str(value)+' SNPs extracted: '+'\n')
    for key in best_params:
        f.write('{'+str(key)+','+str(best_params[key])+'}\n')
    f.close()

    f = open(file_prefix+'_bestaccuracy.txt','a+')
    f.write('#--------------------------------------------------------------------#\n')
    f.write('Accuracy for top '+str(value)+' SNPs extracted: '+str(best_accuracy)+'\n')
    f.close()

    f = open(file_prefix+'_f1score.txt','a+')
    f.write('#--------------------------------------------------------------------#\n')
    f.write('F1 Score for top '+str(value)+' SNPs extracted: '+str(f1score)+'\n')
    f.close()

    f = open(file_prefix+'_caseaccuracy.txt','a+')
    f.write('#--------------------------------------------------------------------#\n')
    f.write('Positive Accuracy for top '+str(value)+' SNPs extracted: '+str(case_accuracy)+'\n')
    f.close()

    f = open(file_prefix+'_controlaccuracy.txt','a+')
    f.write('#--------------------------------------------------------------------#\n')
    f.write('Negative Accuracy for top '+str(value)+' SNPs extracted: '+str(control_accuracy)+'\n')
    f.close()

    f = open(file_prefix+'_trainingaccuracy.txt','a+')
    f.write('#--------------------------------------------------------------------#\n')
    f.write('Training Accuracy for top '+str(value)+' SNPs extracted: '+str(round(best_training_accuracy,3))+'\n')
    f.close()

    f = open(file_prefix+'_confusionmatrices.txt','a+')
    f.write('#--------------------------------------------------------------------#\n')
    f.write('Confusion matrix for top '+str(value)+' SNPs extracted: \n'+str(best_cnfmat)+'\n')
    f.close()


if __name__ == '__main__':

    print("HI")
