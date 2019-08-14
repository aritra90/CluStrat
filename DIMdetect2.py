import matplotlib
matplotlib.use('Agg')
from plinkio import plinkfile
import parseplink as pp
import numpy as np

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.linear_model import Ridge, Lasso, ElasticNet, LinearRegression, RidgeCV
from sklearn.svm import SVC, SVR
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import PolynomialFeatures

from sklearn.metrics import accuracy_score, f1_score, confusion_matrix
from sklearn.metrics import precision_recall_fscore_support as prfs
from sklearn.model_selection import train_test_split, KFold

import os, datetime, random, sys, shutil, itertools, math, argparse, time
import matplotlib.pyplot as plt

import kernelize as K
import mahalanobis as MH
import crossvalidation as CV
import classificationGWAS as cy
import outputgenerator

# AB: Put the plot in a PDF along with other results.
def get_conf(model, X, y, vis=False):
    preds = model.predict(X)
    conf = confusion_matrix(y, preds.round())
    if (vis):
        labels = [0, 1]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        cax = ax.matshow(conf)
        plt.title('Confusion matrix of the classifier')
        fig.colorbar(cax)
        ax.set_xticklabels([''] + labels)
        ax.set_yticklabels([''] + labels)
        plt.xlabel('Predicted')
        plt.ylabel('True')
        plt.show()

    return conf

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-tr", "--train", dest='trainhandle', action='store', help="Put the path to training set file.",
                        metavar="TRAIN")
    parser.add_argument("-te", "--test", dest='testhandle', action='store', help="Put the path to test set file.",
                        metavar="TEST")
    parser.add_argument("-pf", "--prefix", dest='prefix', action='store', help="Put the data prefix.",
                        metavar="PREFIX")
    parser.add_argument("-cl", "--classifier", dest='classifiertype', action='store', help="Put the classifier type. See readme for types.",
                        metavar="CLASS")
    # parser.add_argument("-dir", "--directory", dest='outdir', action='store', help="Put the name of the output directory.",
    #                     metavar="DIR")
    parser.add_argument("-cv", "--crossval", dest='crossval', action='store', help="Put the cross-validation type (kfold or loo).",
                        metavar="CV")
    parser.add_argument("-pval", "--value", dest='value', action='store', help="Put the number of p values kept.",
                        metavar="PVAL")
    parser.add_argument("-knel", "--kernel", dest='kernel', action='store', help="Put 1 or 0 for kernelizing or not.",
                        metavar="KNEL")
    args = parser.parse_args()
    return args.trainhandle, args.testhandle, args.prefix, args.classifiertype, args.crossval, args.value, args.kernel

def read_handlers(trainhandle, testhandle):
    print(trainhandle)
    print(testhandle)
    plink_file = plinkfile.open(trainhandle)
    G_train, train_pheno = pp.read_n_parse(plink_file)
    traincases = np.where(train_pheno == 1)[0]
    traincontrols = np.where(train_pheno == 0)[0]
    print('Training Case to Control Ratio: ')
    print(str(len(traincases))+":"+str(len(traincontrols)))

    plink_file = plinkfile.open(testhandle)
    G_test, test_pheno = pp.read_n_parse(plink_file)
    testcases = np.where(test_pheno == 1)[0]
    testcontrols = np.where(test_pheno == 0)[0]
    print('Testing Case to Control Ratio: ')
    print(str(len(testcases))+":"+str(len(testcontrols)))
    return G_train, G_test, train_pheno, test_pheno

def set_parameters(classifiertype):

    #--------------------------------------------------------------------#
    if classifiertype == 'LDA':
        # Parameters
        # Solver to use (options: 'svd', 'lsqr', and 'eigen')
        solver_LDA = ['svd']#,'lsqr','eigen']
        # Shrinkage parameter (options: None, 'auto' and 0 < float < 1)
        shrinkage = [None]
        # Class priors
        priors = [None]
        # Number of components for dimensionality reduction
        n_components = [None]
        # Class covariance matrix (only in 'svd' solver)
        store_covariance = [False]
        # Threshold for rank estimation (only in 'svd' solver)
        tol_LDA = [1.05e-4]
        iterations = list(itertools.product(solver_LDA,shrinkage,priors,n_components,store_covariance,tol_LDA))
    #--------------------------------------------------------------------#
    elif classifiertype == 'QDA':
        # Parameters
        # Class priors
        priors = [None]
        # Regularization
        reg_param = [1.0]
        # reg_param = [0.0,0.5,1.0,10.0,100.0]
        # Class covariance matrix
        store_covariance = [False]
        # Threshold for rank estimation
        tol_QDA = [1.05e-3]
        # tol_QDA = [1.05e-4,1.05e-3,1.05e-5]
        iterations = list(itertools.product(priors,reg_param,store_covariance,tol_QDA))
    #--------------------------------------------------------------------#
    elif classifiertype == 'PolyReg':
        # Parameters
        # Intercept of the model
        fit_intercept = [True]
        # Normalize (subtract the mean and divide by L2 norm)
        normalize = [False]
        # Degree of polynomial features
        degree = [3]
        # Interaction of features
        interaction_only = [True]
        # Include bias column (feature with all zeros)
        include_bias = [True]
        iterations = list(itertools.product(fit_intercept,normalize,degree,interaction_only,include_bias))
    #--------------------------------------------------------------------#
    elif classifiertype == 'Ridge':
        # Parameters
        # Regularization strength
        alpha = [0.001,0.01,0.1,0.5,1.0]
        # Intercept of the model
        fit_intercept = [True]
        # Normalize (subtract the mean and divide by L2 norm)
        normalize = [False]
        # Max iterations for the solver
        max_iter = [None]
        # Precision of solution
        tol_Ridge = [0.001,0.01,0.1]
        # Solver to use (options: 'auto', 'svd', 'cholesky', 'lsqr','sparse_cg', 'sag', and 'saga')
        solver_Ridge = ['auto']
        # Random seed
        random_state = [None]
        iterations = list(itertools.product(alpha,fit_intercept,normalize,max_iter,tol_Ridge,solver_Ridge,random_state))
    #--------------------------------------------------------------------#
    elif classifiertype == 'Lasso':
        # Parameters
        # Regularization strength
        alpha = [0.001,0.01,0.1,1.0]
        # Intercept of the model
        fit_intercept = [True]
        # Normalize (subtract the mean and divide by L2 norm)
        normalize = [False]
        # Max iterations for the solver
        max_iter = [100,500,1000,2000]
        # Precision of solution
        tol_Lasso = [0.001,0.01,0.1]
        # Set coefficients to be positive
        positive = [False]
        # Random seed
        random_state = [None]
        # Convergence method (options: 'cyclic' or 'random')
        selection = ['cyclic','random']
        iterations = list(itertools.product(alpha,fit_intercept,normalize,max_iter,tol_Lasso,positive,random_state,selection))
    #--------------------------------------------------------------------#
    elif classifiertype == 'Elastic':
        # Parameters
        # Regularization strength
        alpha = [0.001,0.01,0.1,1.0]
        # Mixing parameter
        l1_ratio = [0.1,0.5,0.7]
        # Intercept of the model
        fit_intercept = [True]
        # Normalize (subtract the mean and divide by L2 norm)
        normalize = [False]
        # Max iterations for the solver
        max_iter = [100,500,1000,2000]
        # Precision of solution
        tol_Elastic = [0.001,0.01,0.1]
        # Set coefficients to be positive
        positive = [False]
        # Random seed
        random_state = [None]
        # Convergence method (options: 'cyclic' or 'random')
        selection = ['cyclic','random']
        iterations = list(itertools.product(alpha,l1_ratio,fit_intercept,normalize,max_iter,tol_Elastic,positive,random_state,selection))
    #--------------------------------------------------------------------#
    elif classifiertype == 'SVM':
        # Parameters
        # Penalty parameter of error term
        C = [0.001,0.1,0.5,1.0]
        # Type of kernel to use (options: 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed')
        kernel = ['linear']
        # Degree of the polynomial kernel
        degree = [3]
        # Kernel coefficient for 'rbf', 'poly', and 'sigmoid'
        gamma = ['auto']
        # Independent term in kernel function for 'poly' and 'sigmoid'
        coef0 = [0.0]
        # Shrinking heuristic
        shrinking = [True]
        # Enable probability estimates
        probability = [False]
        # Stopping criterion
        tol_SVM = [1e-3,1e-2,1e-1]
        # Verbose output
        verbose = [False]
        # Limit on iterations within solver
        max_iter = [-1]
        # One versus rest decision boundary or one versus one (options: 'ovr' and 'ovo')
        decision_function_shape = ['ovr','ovo']
        # Pseudo random number generator seed
        random_state = [None]
        iterations = list(itertools.product(C,kernel,degree,gamma,coef0,shrinking,probability,tol_SVM,verbose,max_iter,decision_function_shape,random_state))
    #--------------------------------------------------------------------#
    elif classifiertype == 'SVR':
        # Parameters
        # Penalty parameter of error term
        C = [0.001,0.5,1.0]
        # Type of kernel to use (options: 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed')
        kernel = ['linear']
        # Degree of the polynomial kernel
        degree = [3]
        # Kernel coefficient for 'rbf', 'poly', and 'sigmoid'
        gamma = ['auto']
        # Independent term in kernel function for 'poly' and 'sigmoid'
        coef0 = [0.0,0.5,1.0]
        # Shrinking heuristic
        shrinking = [True]
        # Stopping criterion
        tol_SVR = [1e-3,1e-2]
        # Verbose output
        verbose = [False]
        # Limit on iterations within solver
        max_iter = [-1]
        # Epsilon for the margin
        epsilon = [0.1,0.5]
        iterations = list(itertools.product(C,kernel,degree,gamma,coef0,shrinking,tol_SVR,verbose,max_iter,epsilon))
    #--------------------------------------------------------------------#
    elif classifiertype == 'kSVM':
        # Parameters
        # Penalty parameter of error term
        C = [0.001,0.5,1.0]
        # C = [1]
        # Type of kernel to use (options: 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed')
        # kernel = ['poly','rbf','sigmoid']
        kernel = ['poly']
        # Degree of the polynomial kernel
        degree = [8]
        # Kernel coefficient for 'rbf', 'poly', and 'sigmoid'
        # gamma = ['auto','scale']
        gamma = ['auto']
        # Independent term in kernel function for 'poly' and 'sigmoid'
        # coef0 = [0.0,0.5,1.0]
        coef0 = [0.0]
        # Shrinking heuristic
        shrinking = [True]
        # Enable probability estimates
        probability = [False]
        # Stopping criterion
        # tol_kSVM = [1e-3,1e-2]
        tol_kSVM = [1e-3]
        # Verbose output
        verbose = [False]
        # Limit on iterations within solver
        max_iter = [-1]
        # One versus rest decision boundary or one versus one (options: 'ovr' and 'ovo')
        decision_function_shape = ['ovr']
        # Pseudo random number generator seed
        random_state = [None]
        iterations = list(itertools.product(C,kernel,degree,gamma,coef0,shrinking,probability,tol_kSVM,verbose,max_iter,decision_function_shape,random_state))
    #--------------------------------------------------------------------#
    elif classifiertype == 'RidgeSVM':
        # Parameters
        # Regularization strength
        alpha = [1.0]
        # Intercept of the model
        fit_intercept = [True]
        # Normalize (subtract the mean and divide by L2 norm)
        normalize = [False]
        # Max iterations for the solver
        max_iter = [None]
        # Precision of solution
        tol_RidgeSVM = [1e-3]
        # Solver to use (options: 'auto', 'svd', 'cholesky', 'lsqr','sparse_cg', 'sag', and 'saga')
        solver_RidgeSVM = ['svd']
        # Random seed
        random_state = [None]
        iterations = list(itertools.product(alpha,fit_intercept,normalize,max_iter,tol_RidgeSVM,solver_RidgeSVM,random_state))
    #--------------------------------------------------------------------#
    elif classifiertype == 'RFESVM':
        # Parameters
        # Penalty parameter of error term
        # C = [0.001,0.1,0.5,1.0]
        C = [0.001]
        # Number of features to select (None is actually half)
        n_features_to_select = [None]
        # Number of features to remove at each iteration
        step = [1]
        # Verbosity of output
        verbose = [0]
        iterations = list(itertools.product(C,n_features_to_select,step,verbose))
    #--------------------------------------------------------------------#
    elif classifiertype == 'RidgeCV':
        # Parameters
        # Array of alpha values to try with built-in cross-validation
        alphas = [[0.0001,0.001,0.1,1.0,10.0,100.0]]
        # Calculate intercept or not (not when data is centered)
        fit_intercept = [True]
        # Ignored when fit_intercept is False
        normalize = [False]
        # Scorer callable object
        scoring = [None]
        # Cross-Validation scheme (None = LOOCV,integer for number of folds, others...)
        cv = [None]
        # Generalized CV (see documentation)
        gcv_mode = [None]
        # Cross validation values corresponding to each alpha
        store_cv_values = [True]
        iterations = list(itertools.product(alphas,fit_intercept,normalize,scoring,cv,gcv_mode,store_cv_values))
    #--------------------------------------------------------------------#
    elif classifiertype == 'RandomForest':
        # Parameters
        # Number of trees in the forest
        n_estimators = [10,50,100]
        # Function to measure quality of the split (options: 'gini' or 'entropy')
        criterion = ['gini','entropy']
        # Max depth of the trees
        max_depth = [None]
        # Minimum number of samples required to split a node
        min_samples_split = [2]
        # Minimum number of samples required to be at leaf node
        min_samples_leaf = [1]
        # Minimum weighted fraction of the sum total of weights to be at leaf node
        min_weight_fraction_leaf = [0.0]
        # Number of features when looking at each split (options: int, float, 'auto', 'sqrt', 'log2', and None)
        max_features = ['auto']
        # Grow trees with this many 'best' nodes
        max_leaf_nodes = [None]
        # A node will split if it induces a decrease of impurity greater than this value
        min_impurity_decrease = [0.0]
        # Threshold for stopping in tree growth
        min_impurity_split = [None]
        # Bootstrap samples when building tree
        bootstrap = [True]
        # Out-of-Bag samples to estimate accuracy
        oob_score = [False]
        # Number of jobs to run in parallel for fit and predict
        n_jobs = [None]
        # Random seed
        random_state = [None]
        # Verbosity of the output
        verbose = [0]
        # Weights for the classes
        class_weight = [None]
        iterations = list(itertools.product(n_estimators,criterion,max_depth,min_samples_split,min_samples_leaf,min_weight_fraction_leaf,max_features,
                                            max_leaf_nodes,min_impurity_decrease,min_impurity_split,bootstrap,oob_score,n_jobs,random_state,verbose,class_weight))
    #--------------------------------------------------------------------#
    else:
        print('Not a classifier covered in this script. See MLGWAS.rmd for the options...')
        sys.exit(0)
    #--------------------------------------------------------------------#

    return iterations

if __name__ == '__main__':
    total_time = time.time()

    parse_time = time.time()
    ########################### Parsing Input Files ###########################
    trainhandle, testhandle, prefix, classifiertype, crossval, value, kernel = parse_arguments()
    assert os.path.isdir(trainhandle) == os.path.isdir(testhandle), 'Either both are directories or both are not'
    train_handles = []
    test_handles = []
    if os.path.isdir(trainhandle):
        files_train = os.listdir(trainhandle)
        files_test = os.listdir(testhandle)
        train_handles = list(set(map(lambda x: ''.join(x.split('.')[0:-1]), files_train)))
        test_handles = list(set(map(lambda x: ''.join(x.split('.')[0:-1]), files_test)))
        assert len(train_handles) == len(test_handles)
    else:
        train_handles = [trainhandle]
        test_handles = [testhandle]
    ########################## End of Parsing Files ###########################
    t0= "Parse time: "+str(time.time()-parse_time)+" seconds (wall time)"

    # mlstuff_time = time.time()
    ################################ ML Stuff #################################
    for trainhandle, testhandle in zip(train_handles, test_handles):
        dataread_time = time.time()
        G_train, G_test, train_pheno, test_pheno = read_handlers(trainhandle, testhandle)
        t_load = "Data Load Computation Time: "+ str(time.time()-dataread_time)+" seconds (wall time)"
        num_indivs = len(G_test) + len(G_train)
        n_snps = G_test.shape[1]

        # CODE KERNEL OPTION AS A FLAG IN PARSE ARGS LATER
        # flag = 1

        # No Kernels
        if int(kernel) == 0:
            # Mean center the data about the columns
            means = np.nanmean(G_train,axis=0)
            G_train = G_train - means

            # Mean center the data about the columns
            means = np.nanmean(G_test,axis=0)
            G_test = G_test - means

            trainingdata = G_train
            testingdata = G_test
            t1 = "No Kernel Used..."
            t_inv = "No Inverse Used..."
        # Basic Kernel (just inner product) --> weighted linear kernel (Wang et al. 2014)
        elif int(kernel) == 1:
            kernel_time = time.time()
            # Get the weights of the training set
            the_weights = K.get_weights(G_train)
            # Kernelize the training and testing data (w.r.t. the training data)
            # Mean centering is done inside this function
            kG_train = K.kernelize(G_train,G_train,the_weights)
            kG_test = K.kernelize(G_test,G_train,the_weights)
            t1 = "Kernel (Basic) Computation Time: "+ str(time.time()-kernel_time)+" seconds (wall time)"
            t_inv = "No Inverse Used..."
            trainingdata = kG_train
            testingdata = kG_test
        # Kernel using Mahalanobis Distance as weights
        else:
            kernel_time = time.time()
            # Get the Mahalanobis Distance SNP covariance matrix
            # Mean centering is done inside this function
            inverse_time = time.time()
            S_inv = MH.MDist(G_train)
            t_inv = "Inverse (Mahalanobis) Computation Time: "+ str(time.time()-inverse_time)+" seconds (wall time)"
            # Kernelize the training and testing data (w.r.t. the training data)
            kG_train = MH.kernelize(G_train, S_inv, G_train)
            kG_test = MH.kernelize(G_test, S_inv, G_train)
            t1 = "Kernel (Mahalanobis) Computation Time: "+ str(time.time()-kernel_time)+" seconds (wall time)"

            trainingdata = kG_train
            testingdata = kG_test

        # Set parameters first so we can loop through various lists as well
        iterations = set_parameters(classifiertype)

        # The format of the list 'elem' will depend on the algorithm...
        # Iterate through different options of parameters (cross validation step)
        crossval_time = time.time()
        for elem in iterations:
            # Keep track of accuracy for plotting (keeping the best across all parameters)
            best_training_accuracy = 0
            params_dict = {}
            best_params = {}

            scores, params_dict = CV.crossvalidation(trainingdata, train_pheno, classifiertype, crossval, elem)

            if scores[0] > best_training_accuracy:
                best_training_accuracy = scores[0]
                # Store corresponding parameters
                best_params = params_dict
                best_param_arr = elem
        t2 = "Cross Validation Time: "+str(time.time()-crossval_time)+" seconds (wall time)"

        classification_time = time.time()
        # Model fit and classify
        ytestpred = cy.classify(trainingdata, train_pheno, testingdata, best_param_arr, classifiertype)
        t3 = "Classification Time: "+str(time.time()-classification_time)+" seconds (wall time)"

        # t4 = "Learning+Classification Time: "+ str(time.time()-mlstuff_time)+" seconds (wall time)"
    ################################ ML Stuff #################################
        outputgenerator.generate_output(test_pheno, ytestpred, classifiertype, prefix, value, best_training_accuracy, best_params)
        t5 = "Total time: "+str(time.time()-total_time)+" seconds (wall time)"

        f = open("timing_Info.txt",'a+')
        f.write('#--------------------------------------------------------------------#\n')
        f.write('Dataset dimensions: '+str(num_indivs)+'x'+str(n_snps)+'\n')
        f.write(t0+'\n')
        f.write(t_load+'\n')
        f.write(t1+'\n')
        f.write(t_inv+'\n')
        f.write(t2+'\n')
        f.write(t3+'\n')
        # f.write(t4+'\n')
        f.write(t5+'\n\n')

############################## End of File ################################
