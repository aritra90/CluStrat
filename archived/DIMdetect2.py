import matplotlib
matplotlib.use('Agg')
from plinkio import plinkfile
import parseplink as pp
import numpy as np
import argparse

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.linear_model import Ridge, Lasso, ElasticNet, LinearRegression, RidgeCV
from sklearn.svm import SVC, SVR
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import PolynomialFeatures

from sklearn.metrics import accuracy_score, f1_score, confusion_matrix
from sklearn.metrics import precision_recall_fscore_support as prfs
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import classification_report

import os
import report_generator as rp
import matplotlib.pyplot as plt
import datetime
import random
import sys
import shutil
import itertools

def cross_validation_custom(X, y, model_init, crossval, verbose=False):
    if crossval == 'kfold':
        print('Cross-validation type: K-Fold')
        # k fold cross validation (k = 5)
        NUM_FOLDS = 5
    else:
        print('Cross-validation type: Leave-one-out')
        # Leave one out cross validation (k = n)
        NUM_FOLDS = X.shape[0]
    kf = KFold(n_splits=NUM_FOLDS, shuffle=True)
    score = 0
    i = 1
    print(X.shape, y.shape)
    f_1_score_micro = 0
    f_1_score_macro = 0
    f_1_score_wt = 0
    conf = None
    for train_index, test_index in kf.split(X):
        model = model_init
        # print(i) if verbose else False
        i += 1
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        if classifiertype == 'PolyReg':
            poly = PolynomialFeatures(degree=2,interaction_only=True)
            X_train_ply = poly.fit_transform(X_train)
            X_test_ply = poly.fit_transform(X_test)
            model.fit(X_train_ply, y_train)
            preds = model.predict(X_test_ply)
            # discretization of the continuous predictions
            preds = np.asarray([0 if counter < 1 else 1 for counter in preds])
            preds = np.where(preds.round() >= 1, 1, 0)
        else:
            model.fit(X_train, y_train)
            preds = model.predict(X_test)
            preds = np.where(preds.round() >= 1, 1, 0)

        # print(preds.round())

        score += accuracy_score(y_test, preds)
        f_1_score_micro += f1_score(y_test, preds, average='micro')
        f_1_score_macro += f1_score(y_test, preds, average='macro')
        f_1_score_wt += f1_score(y_test, preds, average='weighted')
        # if conf is None:
        #     conf = get_conf(model, X_test, y_test)
        # else:
        #     print(conf)
        #     conf += get_conf(model, X_test, y_test)
    print('Cross Validation accuracy score:', score / NUM_FOLDS)
    print('Cross Validation micro f-1 score:', f_1_score_micro / NUM_FOLDS)
    print('Cross Validation macro f-1 score:', f_1_score_macro / NUM_FOLDS)
    print('Cross Validation weighted f-1 score:', f_1_score_wt / NUM_FOLDS)
    scores = [score / NUM_FOLDS, f_1_score_micro / NUM_FOLDS, f_1_score_macro / NUM_FOLDS, f_1_score_wt / NUM_FOLDS]
    return scores, model

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
    args = parser.parse_args()
    return args.trainhandle, args.testhandle, args.prefix, args.classifiertype, args.crossval, args.value

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

if __name__ == '__main__':
    ########################### Parsing Input Files ###########################
    trainhandle, testhandle, prefix, classifiertype, crossval, value = parse_arguments()
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


    ################################ ML Stuff #################################
    for trainhandle, testhandle in zip(train_handles, test_handles):
        G_train, G_test, train_pheno, test_pheno = read_handlers(trainhandle, testhandle)

        # Set parameters first so we can loop through various lists as well
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
            reg_param = [0.0,0.5,1.0,10.0,100.0]
            # Class covariance matrix
            store_covariance = [False]
            # Threshold for rank estimation
            tol_QDA = [1.05e-4,1.05e-3,1.05e-5]

            iterations = list(itertools.product(priors,reg_param,store_covariance,tol_QDA))

        #--------------------------------------------------------------------#
        elif classifiertype == 'PolyReg':
            # Parameters
            # Intercept of the model
            fit_intercept = [True]
            # Normalize (subtract the mean and divide by L2 norm)
            normalize = [False]
            # Degree of polynomial features
            degree = [2]
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
            alpha = [0.001,0.01,0.1,1.0]
            # Intercept of the model
            fit_intercept = [True]
            # Normalize (subtract the mean and divide by L2 norm)
            normalize = [False]
            # Max iterations for the solver
            max_iter = [None]
            # Precision of solution
            tol_RidgeSVM = [1e-3,1e-2]
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

        # This for loop will wrap around all the code below...
        # The format of the list 'elem' will depend on the algorithm...
        # Iterate through different options of parameters
        for elem in iterations:
            # Keep track of accuracy for plotting (keeping the best across all parameters)
            best_training_accuracy = 0
            params_dict = {}
            best_params = {}
            #--------------------------------------------------------------------#
            if classifiertype == 'LDA':
                print('Model Type: LDA')
                # Parameters
                # Solver to use (options: 'svd', 'lsqr', and 'eigen')
                solver_LDA = elem[0]
                # Shrinkage parameter (options: None, 'auto' and 0 < float < 1)
                shrinkage = elem[1]
                # Class priors
                priors = elem[2]
                # Number of components for dimensionality reduction
                n_components = elem[3]
                # Class covariance matrix (only in 'svd' solver)
                store_covariance = elem[4]
                # Threshold for rank estimation (only in 'svd' solver)
                tol_LDA = elem[5]

                params_dict = {'solver':solver_LDA,'shrinkage':shrinkage,'tolerance':tol_LDA}

                model_init = LinearDiscriminantAnalysis(solver=solver_LDA,shrinkage=shrinkage,tol=tol_LDA)

                scores, model = cross_validation_custom(G_train, train_pheno,model_init, crossval)

                if scores[0] > best_training_accuracy:
                    best_training_accuracy = scores[0]
                    # Store corresponding parameters
                    best_params = params_dict
                    best_param_arr = elem
                #     # Store corresponding training accuracy and confusion matrix
                #     best_training_accuracy = round(scores[1],3)
                #     best_cnfmat = cnfmat
                #     # Also store the corresponding case(postive) and control(negative) accuracies
                #     # TN/(TN+FP)
                #     control_accuracy = (float(cnfmat[0][0]))/(cnfmat[0][0]+cnfmat[0][1])
                #     # TP/(TP+FN)
                #     case_accuracy = (float(cnfmat[1][1]))/(cnfmat[1][0]+cnfmat[1][1])
                #     # Store corresponding f1 score
                #     f1score = round(stats_var[1],3)
                # # print(best_accuracy)

            #--------------------------------------------------------------------#
            elif classifiertype == 'QDA':
                print('Model Type: QDA')
                # Parameters
                # Class priors
                priors = elem[0]
                # Regularization
                reg_param = elem[1]
                # Class covariance matrix
                store_covariance = elem[2]
                # Threshold for rank estimation
                tol_QDA = elem[3]

                model_init = QuadraticDiscriminantAnalysis(reg_param=reg_param,tol=tol_QDA)

                scores, model = cross_validation_custom(G_train, train_pheno,model_init, crossval)

                model = QuadraticDiscriminantAnalysis(reg_param=reg_param,tol=tol_QDA)

                params_dict = {'regularization':reg_param,'tolerance':tol_QDA}

                if scores[0] > best_training_accuracy:
                    best_training_accuracy = scores[0]
                    # Store corresponding parameters
                    best_params = params_dict
                    best_param_arr = elem

            #--------------------------------------------------------------------#
            elif classifiertype == 'PolyReg':
                print('Model Type: Polynomial Regression')
                # Parameters
                # Intercept of the model
                fit_intercept = elem[0]
                # Normalize (subtract the mean and divide by L2 norm)
                normalize = elem[1]
                # Degree of polynomial features
                degree = elem[2]
                # Interaction of features
                interaction_only = elem[3]
                # Include bias column (feature with all zeros)
                include_bias = elem[4]

                model_init = LinearRegression(fit_intercept=fit_intercept)

                scores, model = cross_validation_custom(G_train, train_pheno,model_init, crossval)

                params_dict = {'fit_intercept':fit_intercept,'degree':degree,'interaction_only':interaction_only,'include_bias':include_bias}

                if scores[0] > best_training_accuracy:
                    best_training_accuracy = scores[0]
                    # Store corresponding parameters
                    best_params = params_dict
                    best_param_arr = elem
            #--------------------------------------------------------------------#
            elif classifiertype == 'Ridge':
                print('Model Type: Ridge Regression')
                # Parameters
                # Regularization strength
                alpha = elem[0]
                # Intercept of the model
                fit_intercept = elem[1]
                # Normalize (subtract the mean and divide by L2 norm)
                normalize = elem[2]
                # Max iterations for the solver
                max_iter = elem[3]
                # Precision of solution
                tol_Ridge = elem[4]
                # Solver to use (options: 'auto', 'svd', 'cholesky', 'lsqr','sparse_cg', 'sag', and 'saga')
                solver_Ridge = elem[5]
                # Random seed
                random_state = elem[6]


                model_init = Ridge(alpha=alpha,fit_intercept=fit_intercept,max_iter=max_iter,tol=tol_Ridge,solver=solver_Ridge)

                scores, model = cross_validation_custom(G_train, train_pheno,model_init, crossval)

                params_dict = {'alpha':alpha,'fit_intercept':fit_intercept,'max_iter':max_iter,'tolerance':tol_Ridge,'solver':solver_Ridge}

                if scores[0] > best_training_accuracy:
                    best_training_accuracy = scores[0]
                    # Store corresponding parameters
                    best_params = params_dict
                    best_param_arr = elem
            #--------------------------------------------------------------------#
            elif classifiertype == 'Lasso':
                print('Model Type: Lasso Regression')
                # Parameters
                # Regularization strength
                alpha = elem[0]
                # Intercept of the model
                fit_intercept = elem[1]
                # Normalize (subtract the mean and divide by L2 norm)
                normalize = elem[2]
                # Max iterations for the solver
                max_iter = elem[3]
                # Precision of solution
                tol_Lasso = elem[4]
                # Set coefficients to be positive
                positive = elem[5]
                # Random seed
                random_state = elem[6]
                # Convergence method (options: 'cyclic' or 'random')
                selection = elem[7]

                model_init = Lasso(alpha=alpha,fit_intercept=fit_intercept,max_iter=max_iter,tol=tol_Lasso,positive=positive,selection=selection)

                scores, model = cross_validation_custom(G_train, train_pheno,model_init, crossval)

                params_dict = {'alpha':alpha,'fit_intercept':fit_intercept,'max_iter':max_iter,'tolerance':tol_Lasso,'positive':positive,'selection':selection}

                if scores[0] > best_training_accuracy:
                    best_training_accuracy = scores[0]
                    # Store corresponding parameters
                    best_params = params_dict
                    best_param_arr = elem
            #--------------------------------------------------------------------#
            elif classifiertype == 'Elastic':
                print('Model Type: Elastic Net')
                # Parameters
                # Regularization strength
                alpha = elem[0]
                # Mixing parameter
                l1_ratio = elem[1]
                # Intercept of the model
                fit_intercept = elem[2]
                # Normalize (subtract the mean and divide by L2 norm)
                normalize = elem[3]
                # Max iterations for the solver
                max_iter = elem[4]
                # Precision of solution
                tol_Elastic = elem[5]
                # Set coefficients to be positive
                positive = elem[6]
                # Random seed
                random_state = elem[7]
                # Convergence method (options: 'cyclic' or 'random')
                selection = elem[8]


                model_init = ElasticNet(alpha=alpha,l1_ratio=l1_ratio,fit_intercept=fit_intercept,max_iter=max_iter,tol=tol_Elastic,positive=positive,selection=selection)

                scores, model = cross_validation_custom(G_train, train_pheno,model_init, crossval)

                params_dict = {'alpha':alpha,'l1_ratio':l1_ratio,'fit_intercept':fit_intercept,'max_iter':max_iter,'tolerance':tol_Elastic,'positive':positive,'selection':selection}

                if scores[0] > best_training_accuracy:
                    best_training_accuracy = scores[0]
                    # Store corresponding parameters
                    best_params = params_dict
                    best_param_arr = elem
            #--------------------------------------------------------------------#
            elif classifiertype == 'SVM':
                print('Model Type: SVM')
                # Parameters
                # Penalty parameter of error term
                C = elem[0]
                # Type of kernel to use (options: 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed')
                kernel = elem[1]
                # Degree of the polynomial kernel
                degree = elem[2]
                # Kernel coefficient for 'rbf', 'poly', and 'sigmoid'
                gamma = elem[3]
                # Independent term in kernel function for 'poly' and 'sigmoid'
                coef0 = elem[4]
                # Shrinking heuristic
                shrinking = elem[5]
                # Enable probability estimates
                probability = elem[6]
                # Stopping criterion
                tol = elem[7]
                # Verbose output
                verbose = elem[8]
                # Limit on iterations within solver
                max_iter = elem[9]
                # One versus rest decision boundary or one versus one (options: 'ovr' and 'ovo')
                decision_function_shape = elem[10]
                # Pseudo random number generator seed
                random_state = elem[11]

                model_init = SVC(C=C,kernel=kernel,degree=degree,gamma=gamma,coef0=coef0,
                                shrinking=shrinking,probability=probability,tol=tol,verbose=verbose,
                                max_iter=max_iter,decision_function_shape=decision_function_shape,
                                random_state=random_state)

                scores, model = cross_validation_custom(G_train, train_pheno,model_init, crossval)

                params_dict = {'C':C,'kernel':kernel,'degree':degree,'gamma':gamma,
                                'coef0':coef0,'shrinking':shrinking,'probability':probability,
                                'tol':tol,'verbose':verbose,'max_iter':max_iter,'decision_function_shape':decision_function_shape,
                                'random_state':random_state}

                if scores[0] > best_training_accuracy:
                    best_training_accuracy = scores[0]
                    # Store corresponding parameters
                    best_params = params_dict
                    best_param_arr = elem

            #--------------------------------------------------------------------#
            elif classifiertype == 'SVR':
                    print('Model Type: SVR')
                    # Parameters
                    # Penalty parameter of error term
                    C = elem[0]
                    # Type of kernel to use (options: 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed')
                    kernel = elem[1]
                    # Degree of the polynomial kernel
                    degree = elem[2]
                    # Kernel coefficient for 'rbf', 'poly', and 'sigmoid'
                    gamma = elem[3]
                    # Independent term in kernel function for 'poly' and 'sigmoid'
                    coef0 = elem[4]
                    # Shrinking heuristic
                    shrinking = elem[5]
                    # Stopping criterion
                    tol_SVR = elem[6]
                    # Verbose output
                    verbose = elem[7]
                    # Limit on iterations within solver
                    max_iter = elem[8]
                    # Epsilon for the margin
                    epsilon = elem[9]

                    model_init = SVR(C=C,kernel=kernel,degree=degree,gamma=gamma,coef0=coef0,
                                    shrinking=shrinking,tol=tol_SVR,verbose=verbose,
                                    max_iter=max_iter,epsilon=epsilon)

                    scores, model = cross_validation_custom(G_train, train_pheno,model_init, crossval)

                    params_dict = {'C':C,'kernel':kernel,'degree':degree,'gamma':gamma,
                                    'coef0':coef0,'shrinking':shrinking,'tol':tol_SVR,'verbose':verbose,
                                    'max_iter':max_iter,'epsilon':epsilon}

                    if scores[0] > best_training_accuracy:
                        best_training_accuracy = scores[0]
                        # Store corresponding parameters
                        best_params = params_dict
                        best_param_arr = elem
            #--------------------------------------------------------------------#
            elif classifiertype == 'kSVM':
                print('Model Type: kernel-SVM')
                # Parameters
                # Penalty parameter of error term
                C = elem[0]
                # Type of kernel to use (options: 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed')
                kernel = elem[1]
                # Degree of the polynomial kernel
                degree = elem[2]
                # Kernel coefficient for 'rbf', 'poly', and 'sigmoid'
                gamma = elem[3]
                # Independent term in kernel function for 'poly' and 'sigmoid'
                coef0 = elem[4]
                # Shrinking heuristic
                shrinking = elem[5]
                # Enable probability estimates
                probability = elem[6]
                # Stopping criterion
                tol = elem[7]
                # Verbose output
                verbose = elem[8]
                # Limit on iterations within solver
                max_iter = elem[9]
                # One versus rest decision boundary or one versus one (options: 'ovr' and 'ovo')
                decision_function_shape = elem[10]
                # Pseudo random number generator seed
                random_state = elem[11]

                model_init = SVC(C=C,kernel=kernel,degree=degree,gamma=gamma,coef0=coef0,
                                shrinking=shrinking,probability=probability,tol=tol,verbose=verbose,
                                max_iter=max_iter,decision_function_shape=decision_function_shape,
                                random_state=random_state)

                scores, model = cross_validation_custom(G_train, train_pheno,model_init, crossval)

                params_dict = {'C':C,'kernel':kernel,'degree':degree,'gamma':gamma,
                                'coef0':coef0,'shrinking':shrinking,'probability':probability,
                                'tol':tol,'verbose':verbose,'max_iter':max_iter,'decision_function_shape':decision_function_shape,
                                'random_state':random_state}

                if scores[0] > best_training_accuracy:
                    best_training_accuracy = scores[0]
                    # Store corresponding parameters
                    best_params = params_dict
                    best_param_arr = elem
            #--------------------------------------------------------------------#
            elif classifiertype == 'RidgeSVM':
                print('Model Type: Ridge-SVM')
                # Parameters
                # Regularization strength
                alpha = elem[0]
                # Intercept of the model
                fit_intercept = elem[1]
                # Normalize (subtract the mean and divide by L2 norm)
                normalize = elem[2]
                # Max iterations for the solver
                max_iter = elem[3]
                # Precision of solution
                tol_RidgeSVM = elem[4]
                # Solver to use (options: 'auto', 'svd', 'cholesky', 'lsqr','sparse_cg', 'sag', and 'saga')
                solver_RidgeSVM = elem[5]
                # Random seed
                random_state = elem[6]


                model_init = Ridge(alpha=alpha,fit_intercept=fit_intercept,max_iter=max_iter,tol=tol_RidgeSVM,solver=solver_RidgeSVM)

                scores, model = cross_validation_custom(G_train, train_pheno,model_init, crossval)

                params_dict = {'alpha':alpha,'fit_intercept':fit_intercept,'max_iter':max_iter,'tolerance':tol_RidgeSVM,'solver':solver_RidgeSVM}

                if scores[0] > best_training_accuracy:
                    best_training_accuracy = scores[0]
                    # Store corresponding parameters
                    best_params = params_dict
                    best_param_arr = elem
            #--------------------------------------------------------------------#
            elif classifiertype == 'RFESVM':
                print('Model Type: RFE-SVM')
                # Parameters
                # Penalty parameter of error term
                C = elem[0]
                # Number of features to select (None is actually half)
                n_features_to_select = elem[1]
                # Number of features to remove at each iteration
                step = elem[2]
                # Verbosity of output
                verbose = elem[3]

                estimator_object = SVC(kernel='linear',C=C)
                model_init = RFE(estimator=estimator_object, n_features_to_select=n_features_to_select,step=step,verbose=verbose)

                scores, model = cross_validation_custom(G_train, train_pheno,model_init, crossval)

                params_dict = {'C':C,'estimator_object':estimator_object,'n_features_to_select':n_features_to_select,'step':step,'verbose':verbose}

                if scores[0] > best_training_accuracy:
                    best_training_accuracy = scores[0]
                    # Store corresponding parameters
                    best_params = params_dict
                    best_param_arr = elem

            #--------------------------------------------------------------------#
            elif classifiertype == 'RidgeCV':
                print('Model Type: RidgeCV')
                # Parameters
                # Array of alpha values to try with built-in cross-validation
                alphas = elem[0]
                # Calculate intercept or not (not when data is centered)
                fit_intercept = elem[1]
                # Ignored when fit_intercept is False
                normalize = elem[2]
                # Scorer callable object
                scoring = elem[3]
                # Cross-Validation scheme (None = LOOCV,integer for number of folds, others...)
                cv = elem[4]
                # Generalized CV (see documentation)
                gcv_mode = elem[5]
                # Cross validation values corresponding to each alpha
                store_cv_values = elem[6]

                # model_init = RidgeCV(alphas=alphas,fit_intercept=fit_intercept,normalize=normalize,scoring=scoring,cv=cv,gcv_mode=gcv_mode,store_cv_values=store_cv_values)

                # Wont need this step but still need to fill 'scores' with appropriate values
                # scores, model = cross_validation_custom(G_train, train_pheno, model_init, crossval)

                # Test based on BEST alphas

                model = RidgeCV(alphas=alphas,fit_intercept=fit_intercept,normalize=normalize,scoring=scoring,cv=cv,gcv_mode=gcv_mode,store_cv_values=store_cv_values)

                params_dict = {'alphas':alphas,'fit_intercept':fit_intercept,'normalize':normalize,'scoring':scoring,'cv':cv,'gcv_mode':gcv_mode,'store_cv_values':store_cv_values}

                model.fit(G_train, train_pheno)

                # Just so it doesnt break given the report filling setup (might be able to remove)
                scores = [model.score(G_train, train_pheno)]*4

                if scores[0] > best_training_accuracy:
                    best_training_accuracy = scores[0]
                    # Store corresponding parameters
                    best_params = params_dict
                    best_param_arr = elem
            #--------------------------------------------------------------------#
            elif classifiertype == 'RandomForest':
                print('Model Type: Random Forest')
                # Parameters
                # Number of trees in the forest
                n_estimators = elem[0]
                # Function to measure quality of the split (options: 'gini' or 'entropy')
                criterion = elem[1]
                # Max depth of the trees
                max_depth = elem[2]
                # Minimum number of samples required to split a node
                min_samples_split = elem[3]
                # Minimum number of samples required to be at leaf node
                min_samples_leaf = elem[4]
                # Minimum weighted fraction of the sum total of weights to be at leaf node
                min_weight_fraction_leaf = elem[5]
                # Number of features when looking at each split (options: int, float, 'auto', 'sqrt', 'log2', and None)
                max_features = elem[6]
                # Grow trees with this many 'best' nodes
                max_leaf_nodes = elem[7]
                # A node will split if it induces a decrease of impurity greater than this value
                min_impurity_decrease = elem[8]
                # Threshold for stopping in tree growth
                min_impurity_split = elem[9]
                # Bootstrap samples when building tree
                bootstrap = elem[10]
                # Out-of-Bag samples to estimate accuracy
                oob_score = elem[11]
                # Number of jobs to run in parallel for fit and predict
                n_jobs = elem[12]
                # Random seed
                random_state = elem[13]
                # Verbosity of the output
                verbose = elem[14]
                # Weights for the classes
                class_weight = elem[15]

                model_init = RandomForestClassifier(n_estimators=n_estimators,criterion=criterion,max_depth=max_depth,
                                                    min_samples_split=min_samples_split,min_samples_leaf=min_samples_leaf,
                                                    min_weight_fraction_leaf=min_weight_fraction_leaf,max_features=max_features,
                                                    max_leaf_nodes=max_leaf_nodes,min_impurity_decrease=min_impurity_decrease,
                                                    min_impurity_split=min_impurity_split,bootstrap=bootstrap,oob_score=oob_score,
                                                    n_jobs=n_jobs,random_state=0,verbose=verbose,class_weight=class_weight)

                scores, model = cross_validation_custom(G_train, train_pheno,model_init, crossval)

                params_dict = {'n_estimators':n_estimators,'criterion':criterion,'max_depth':max_depth,'min_samples_split':min_samples_split,'min_samples_leaf':min_samples_leaf,'min_weight_fraction_leaf':min_weight_fraction_leaf,
                                'max_features':max_features,'max_leaf_nodes':max_leaf_nodes,'min_impurity_decrease':min_impurity_decrease,'min_impurity_split':min_impurity_split,'bootstrap':bootstrap,'oob_score':oob_score,
                                'n_jobs':n_jobs,'random_state':random_state,'verbose':verbose,'class_weight':class_weight}

                if scores[0] > best_training_accuracy:
                    best_training_accuracy = scores[0]
                    # Store corresponding parameters
                    best_params = params_dict
                    best_param_arr = elem
            #--------------------------------------------------------------------#
            else:
                print('Not a classifier covered in this script. See MLGWAS.rmd for the options...')
                sys.exit(0)
            #--------------------------------------------------------------------#

        #--------------------------------------------------------------------#
        if classifiertype == 'LDA':

            # Parameters
            # Solver to use (options: 'svd', 'lsqr', and 'eigen')
            solver_LDA = best_param_arr[0]
            # Shrinkage parameter (options: None, 'auto' and 0 < float < 1)
            shrinkage = best_param_arr[1]
            # Class priors
            priors = best_param_arr[2]
            # Number of components for dimensionality reduction
            n_components = best_param_arr[3]
            # Class covariance matrix (only in 'svd' solver)
            store_covariance = best_param_arr[4]
            # Threshold for rank estimation (only in 'svd' solver)
            tol_LDA = best_param_arr[5]

            model = LinearDiscriminantAnalysis(solver=solver_LDA,shrinkage=shrinkage,tol=tol_LDA)

            params_dict = {'solver':solver_LDA,'shrinkage':shrinkage,'tolerance':tol_LDA}

            model.fit(G_train, train_pheno)
            # test the model on your holdout test data
            ytestpred = model.predict(G_test)
            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
            rp_model = rp.Model(model_name='lda', model_info=[' '], training_metrics=scores,
                           testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                            prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])
        #--------------------------------------------------------------------#
        elif classifiertype == 'QDA':

            # Parameters
            # Class priors
            priors = best_param_arr[0]
            # Regularization
            reg_param = best_param_arr[1]
            # Class covariance matrix
            store_covariance = best_param_arr[2]
            # Threshold for rank estimation
            tol_QDA = best_param_arr[3]

            model = QuadraticDiscriminantAnalysis(reg_param=reg_param,tol=tol_QDA)

            params_dict = {'regularization':reg_param,'tolerance':tol_QDA}

            model.fit(G_train, train_pheno)
            # test the model on your holdout test data
            ytestpred = model.predict(G_test)
            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
            rp_model = rp.Model(model_name='qda', model_info=[' '], training_metrics=scores,
                           testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                            prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])
            # sys.exit(0)
        #--------------------------------------------------------------------#
        elif classifiertype == 'PolyReg':

            # Parameters
            # Intercept of the model
            fit_intercept = best_param_arr[0]
            # Normalize (subtract the mean and divide by L2 norm)
            normalize = best_param_arr[1]
            # Degree of polynomial features
            degree = best_param_arr[2]
            # Interaction of features
            interaction_only = best_param_arr[3]
            # Include bias column (feature with all zeros)
            include_bias = best_param_arr[4]

            model = LinearRegression(fit_intercept=fit_intercept)
            poly = PolynomialFeatures(degree=degree,interaction_only=interaction_only,include_bias=include_bias)

            params_dict = {'fit_intercept':fit_intercept,'degree':degree,'interaction_only':interaction_only,'include_bias':include_bias}

            G_train_ply = poly.fit_transform(G_train)
            G_test_ply = poly.fit_transform(G_test)
            model.fit(G_train_ply, train_pheno)
            # test the model on your holdout test data
            ytestpred = model.predict(G_test_ply)
            # discretization of the continuous predictions
            ytestpred = np.asarray([0 if counter < 1 else 1 for counter in ytestpred])
            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
            # print(ytestpred.round())
            rp_model = rp.Model(model_name='polyreg', model_info=[' '], training_metrics=scores,
                           testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                            prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])

            # sys.exit(0)
        #--------------------------------------------------------------------#
        elif classifiertype == 'Ridge':

            # Parameters
            # Regularization strength
            alpha = best_param_arr[0]
            # Intercept of the model
            fit_intercept = best_param_arr[1]
            # Normalize (subtract the mean and divide by L2 norm)
            normalize = best_param_arr[2]
            # Max iterations for the solver
            max_iter = best_param_arr[3]
            # Precision of solution
            tol_Ridge = best_param_arr[4]
            # Solver to use (options: 'auto', 'svd', 'cholesky', 'lsqr','sparse_cg', 'sag', and 'saga')
            solver_Ridge = best_param_arr[5]
            # Random seed
            random_state = best_param_arr[6]

            model = Ridge(alpha=alpha,fit_intercept=fit_intercept,max_iter=max_iter,tol=tol_Ridge,solver=solver_Ridge)

            params_dict = {'alpha':alpha,'fit_intercept':fit_intercept,'max_iter':max_iter,'tolerance':tol_Ridge,'solver':solver_Ridge}

            model.fit(G_train, train_pheno)
            # test the model on your holdout test data
            ytestpred = model.predict(G_test)

            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)

            rp_model = rp.Model(model_name='ridge', model_info=[' '], training_metrics=scores,
                           testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                            prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])
            # sys.exit(0)
        #--------------------------------------------------------------------#
        elif classifiertype == 'Lasso':

            # Parameters
            # Regularization strength
            alpha = best_param_arr[0]
            # Intercept of the model
            fit_intercept = best_param_arr[1]
            # Normalize (subtract the mean and divide by L2 norm)
            normalize = best_param_arr[2]
            # Max iterations for the solver
            max_iter = best_param_arr[3]
            # Precision of solution
            tol_Lasso = best_param_arr[4]
            # Set coefficients to be positive
            positive = best_param_arr[5]
            # Random seed
            random_state = best_param_arr[6]
            # Convergence method (options: 'cyclic' or 'random')
            selection = best_param_arr[7]

            model = Lasso(alpha=alpha,fit_intercept=fit_intercept,max_iter=max_iter,tol=tol_Lasso,positive=positive,selection=selection)

            params_dict = {'alpha':alpha,'fit_intercept':fit_intercept,'max_iter':max_iter,'tolerance':tol_Lasso,'positive':positive,'selection':selection}

            model.fit(G_train, train_pheno)
            # test the model on your holdout test data
            ytestpred = model.predict(G_test)
            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
            rp_model = rp.Model(model_name='lasso', model_info=[' '], training_metrics=scores,
                           testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                            prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])
            # sys.exit(0)
        #--------------------------------------------------------------------#
        elif classifiertype == 'Elastic':

            # Parameters
            # Regularization strength
            alpha = best_param_arr[0]
            # Mixing parameter
            l1_ratio = best_param_arr[1]
            # Intercept of the model
            fit_intercept = best_param_arr[2]
            # Normalize (subtract the mean and divide by L2 norm)
            normalize = best_param_arr[3]
            # Max iterations for the solver
            max_iter = best_param_arr[4]
            # Precision of solution
            tol_Elastic = best_param_arr[5]
            # Set coefficients to be positive
            positive = best_param_arr[6]
            # Random seed
            random_state = best_param_arr[7]
            # Convergence method (options: 'cyclic' or 'random')
            selection = best_param_arr[8]

            model = ElasticNet(alpha=alpha,l1_ratio=l1_ratio,fit_intercept=fit_intercept,max_iter=max_iter,tol=tol_Elastic,positive=positive,selection=selection)

            params_dict = {'alpha':alpha,'l1_ratio':l1_ratio,'fit_intercept':fit_intercept,'max_iter':max_iter,'tolerance':tol_Elastic,'positive':positive,'selection':selection}

            model.fit(G_train, train_pheno)
            # test the model on your holdout test data
            ytestpred = model.predict(G_test)
            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
            rp_model = rp.Model(model_name='elastic', model_info=[' '], training_metrics=scores,
                           testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                            prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])
            # sys.exit(0)
        #--------------------------------------------------------------------#
        elif classifiertype == 'SVM':

            # Parameters
            # Penalty parameter of error term
            C = best_param_arr[0]
            # Type of kernel to use (options: 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed')
            kernel = best_param_arr[1]
            # Degree of the polynomial kernel
            degree = best_param_arr[2]
            # Kernel coefficient for 'rbf', 'poly', and 'sigmoid'
            gamma = best_param_arr[3]
            # Independent term in kernel function for 'poly' and 'sigmoid'
            coef0 = best_param_arr[4]
            # Shrinking heuristic
            shrinking = best_param_arr[5]
            # Enable probability estimates
            probability = best_param_arr[6]
            # Stopping criterion
            tol = best_param_arr[7]
            # Verbose output
            verbose = best_param_arr[8]
            # Limit on iterations within solver
            max_iter = best_param_arr[9]
            # One versus rest decision boundary or one versus one (options: 'ovr' and 'ovo')
            decision_function_shape = best_param_arr[10]
            # Pseudo random number generator seed
            random_state = best_param_arr[11]

            model = SVC(C=C,kernel=kernel,degree=degree,gamma=gamma,coef0=coef0,
                            shrinking=shrinking,probability=probability,tol=tol,verbose=verbose,
                            max_iter=max_iter,decision_function_shape=decision_function_shape,
                            random_state=random_state)

            params_dict = {'C':C,'kernel':kernel,'degree':degree,'gamma':gamma,
                            'coef0':coef0,'shrinking':shrinking,'probability':probability,
                            'tol':tol,'verbose':verbose,'max_iter':max_iter,'decision_function_shape':decision_function_shape,
                            'random_state':random_state}

            model.fit(G_train, train_pheno)
            # test the model on your holdout test data
            ytestpred = model.predict(G_test)
            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
            rp_model = rp.Model(model_name='svm', model_info=[' '], training_metrics=scores,
                           testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                            prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])
            # sys.exit(0)

        #--------------------------------------------------------------------#
        elif classifiertype == 'SVR':

                # Parameters
                # Penalty parameter of error term
                C = best_param_arr[0]
                # Type of kernel to use (options: 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed')
                kernel = best_param_arr[1]
                # Degree of the polynomial kernel
                degree = best_param_arr[2]
                # Kernel coefficient for 'rbf', 'poly', and 'sigmoid'
                gamma = best_param_arr[3]
                # Independent term in kernel function for 'poly' and 'sigmoid'
                coef0 = best_param_arr[4]
                # Shrinking heuristic
                shrinking = best_param_arr[5]
                # Stopping criterion
                tol_SVR = best_param_arr[6]
                # Verbose output
                verbose = best_param_arr[7]
                # Limit on iterations within solver
                max_iter = best_param_arr[8]
                # Epsilon for the margin
                epsilon = best_param_arr[9]

                model = SVR(C=C,kernel=kernel,degree=degree,gamma=gamma,coef0=coef0,
                                shrinking=shrinking,tol=tol_SVR,verbose=verbose,
                                max_iter=max_iter,epsilon=epsilon)

                params_dict = {'C':C,'kernel':kernel,'degree':degree,'gamma':gamma,
                                'coef0':coef0,'shrinking':shrinking,'tol':tol_SVR,'verbose':verbose,
                                'max_iter':max_iter,'epsilon':epsilon}

                model.fit(G_train, train_pheno)
                # test the model on your holdout test data
                ytestpred = model.predict(G_test)
                ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
                rp_model = rp.Model(model_name='svr', model_info=[' '], training_metrics=scores,
                               testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                                prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                                prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                                prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])
                # sys.exit(0)
        #--------------------------------------------------------------------#
        elif classifiertype == 'kSVM':

            # Parameters
            # Penalty parameter of error term
            C = best_param_arr[0]
            # Type of kernel to use (options: 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed')
            kernel = best_param_arr[1]
            # Degree of the polynomial kernel
            degree = best_param_arr[2]
            # Kernel coefficient for 'rbf', 'poly', and 'sigmoid'
            gamma = best_param_arr[3]
            # Independent term in kernel function for 'poly' and 'sigmoid'
            coef0 = best_param_arr[4]
            # Shrinking heuristic
            shrinking = best_param_arr[5]
            # Enable probability estimates
            probability = best_param_arr[6]
            # Stopping criterion
            tol = best_param_arr[7]
            # Verbose output
            verbose = best_param_arr[8]
            # Limit on iterations within solver
            max_iter = best_param_arr[9]
            # One versus rest decision boundary or one versus one (options: 'ovr' and 'ovo')
            decision_function_shape = best_param_arr[10]
            # Pseudo random number generator seed
            random_state = best_param_arr[11]

            model = SVC(C=C,kernel=kernel,degree=degree,gamma=gamma,coef0=coef0,
                            shrinking=shrinking,probability=probability,tol=tol,verbose=verbose,
                            max_iter=max_iter,decision_function_shape=decision_function_shape,
                            random_state=random_state)

            params_dict = {'C':C,'kernel':kernel,'degree':degree,'gamma':gamma,
                            'coef0':coef0,'shrinking':shrinking,'probability':probability,
                            'tol':tol,'verbose':verbose,'max_iter':max_iter,'decision_function_shape':decision_function_shape,
                            'random_state':random_state}

            model.fit(G_train, train_pheno)
            # test the model on your holdout test data
            ytestpred = model.predict(G_test)
            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
            rp_model = rp.Model(model_name='ksvm', model_info=[' '], training_metrics=scores,
                           testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                            prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])
            # sys.exit(0)
        #--------------------------------------------------------------------#
        elif classifiertype == 'RidgeSVM':

            # Parameters
            # Regularization strength
            alpha = best_param_arr[0]
            # Intercept of the model
            fit_intercept = best_param_arr[1]
            # Normalize (subtract the mean and divide by L2 norm)
            normalize = best_param_arr[2]
            # Max iterations for the solver
            max_iter = best_param_arr[3]
            # Precision of solution
            tol_RidgeSVM = best_param_arr[4]
            # Solver to use (options: 'auto', 'svd', 'cholesky', 'lsqr','sparse_cg', 'sag', and 'saga')
            solver_RidgeSVM = best_param_arr[5]
            # Random seed
            random_state = best_param_arr[6]

            model = Ridge(alpha=alpha,fit_intercept=fit_intercept,max_iter=max_iter,tol=tol_RidgeSVM,solver=solver_RidgeSVM)

            params_dict = {'alpha':alpha,'fit_intercept':fit_intercept,'max_iter':max_iter,'tolerance':tol_RidgeSVM,'solver':solver_RidgeSVM}

            model.fit(G_train, train_pheno)
            # test the model on your holdout test data
            ytestpred = model.predict(G_test)
            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
            rp_model = rp.Model(model_name='ridgesvm', model_info=[' '], training_metrics=scores,
                           testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                            prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])
            # sys.exit(0)
        #--------------------------------------------------------------------#
        elif classifiertype == 'RFESVM':

            # Parameters
            # Penalty parameter of error term
            C = best_param_arr[0]
            # Number of features to select (None is actually half)
            n_features_to_select = best_param_arr[1]
            # Number of features to remove at each iteration
            step = best_param_arr[2]
            # Verbosity of output
            verbose = best_param_arr[3]

            estimator_object = SVC(kernel='linear',C=C)
            model = RFE(estimator=estimator_object, n_features_to_select=n_features_to_select,step=step,verbose=verbose)

            params_dict = {'C':C,'estimator_object':estimator_object,'n_features_to_select':n_features_to_select,'step':step,'verbose':verbose}

            model.fit(G_train, train_pheno)
            # test the model on your holdout test data
            ytestpred = model.predict(G_test)
            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
            rp_model = rp.Model(model_name='rfesvm', model_info=[' '], training_metrics=scores,
                           testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                            prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])
            # sys.exit(0)

        #--------------------------------------------------------------------#
        elif classifiertype == 'RidgeCV':

            # Parameters
            # Array of alpha values to try with built-in cross-validation
            alphas = best_param_arr[0]
            # Calculate intercept or not (not when data is centered)
            fit_intercept = best_param_arr[1]
            # Ignored when fit_intercept is False
            normalize = best_param_arr[2]
            # Scorer callable object
            scoring = best_param_arr[3]
            # Cross-Validation scheme (None = LOOCV,integer for number of folds, others...)
            cv = best_param_arr[4]
            # Generalized CV (see documentation)
            gcv_mode = best_param_arr[5]
            # Cross validation values corresponding to each alpha
            store_cv_values = best_param_arr[6]

            # model_init = RidgeCV(alphas=alphas,fit_intercept=fit_intercept,normalize=normalize,scoring=scoring,cv=cv,gcv_mode=gcv_mode,store_cv_values=store_cv_values)

            # Wont need this step but still need to fill 'scores' with appropriate values
            # scores, model = cross_validation_custom(G_train, train_pheno, model_init, crossval)

            # Test based on BEST alphas

            model = RidgeCV(alphas=alphas,fit_intercept=fit_intercept,normalize=normalize,scoring=scoring,cv=cv,gcv_mode=gcv_mode,store_cv_values=store_cv_values)

            params_dict = {'alphas':alphas,'fit_intercept':fit_intercept,'normalize':normalize,'scoring':scoring,'cv':cv,'gcv_mode':gcv_mode,'store_cv_values':store_cv_values}

            model.fit(G_train, train_pheno)

            # Just so it doesnt break given the report filling setup (might be able to remove)
            scores = [model.score(G_train, train_pheno)]*4

            # test the model on your holdout test data
            ytestpred = model.predict(G_test)
            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
            rp_model = rp.Model(model_name='ridgecv', model_info=[' '], training_metrics=scores,
                           testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                            prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])
            # sys.exit(0)
        #--------------------------------------------------------------------#
        elif classifiertype == 'RandomForest':

            # Parameters
            # Number of trees in the forest
            n_estimators = best_param_arr[0]
            # Function to measure quality of the split (options: 'gini' or 'entropy')
            criterion = best_param_arr[1]
            # Max depth of the trees
            max_depth = best_param_arr[2]
            # Minimum number of samples required to split a node
            min_samples_split = best_param_arr[3]
            # Minimum number of samples required to be at leaf node
            min_samples_leaf = best_param_arr[4]
            # Minimum weighted fraction of the sum total of weights to be at leaf node
            min_weight_fraction_leaf = best_param_arr[5]
            # Number of features when looking at each split (options: int, float, 'auto', 'sqrt', 'log2', and None)
            max_features = best_param_arr[6]
            # Grow trees with this many 'best' nodes
            max_leaf_nodes = best_param_arr[7]
            # A node will split if it induces a decrease of impurity greater than this value
            min_impurity_decrease = best_param_arr[8]
            # Threshold for stopping in tree growth
            min_impurity_split = best_param_arr[9]
            # Bootstrap samples when building tree
            bootstrap = best_param_arr[10]
            # Out-of-Bag samples to estimate accuracy
            oob_score = best_param_arr[11]
            # Number of jobs to run in parallel for fit and predict
            n_jobs = best_param_arr[12]
            # Random seed
            random_state = best_param_arr[13]
            # Verbosity of the output
            verbose = best_param_arr[14]
            # Weights for the classes
            class_weight = best_param_arr[15]

            model = RandomForestClassifier(n_estimators=n_estimators,criterion=criterion,max_depth=max_depth,
                                                min_samples_split=min_samples_split,min_samples_leaf=min_samples_leaf,
                                                min_weight_fraction_leaf=min_weight_fraction_leaf,max_features=max_features,
                                                max_leaf_nodes=max_leaf_nodes,min_impurity_decrease=min_impurity_decrease,
                                                min_impurity_split=min_impurity_split,bootstrap=bootstrap,oob_score=oob_score,
                                                n_jobs=n_jobs,random_state=0,verbose=verbose,class_weight=class_weight)

            params_dict = {'n_estimators':n_estimators,'criterion':criterion,'max_depth':max_depth,'min_samples_split':min_samples_split,'min_samples_leaf':min_samples_leaf,'min_weight_fraction_leaf':min_weight_fraction_leaf,
                            'max_features':max_features,'max_leaf_nodes':max_leaf_nodes,'min_impurity_decrease':min_impurity_decrease,'min_impurity_split':min_impurity_split,'bootstrap':bootstrap,'oob_score':oob_score,
                            'n_jobs':n_jobs,'random_state':random_state,'verbose':verbose,'class_weight':class_weight}

            model.fit(G_train, train_pheno)
            # test the model on your holdout test data
            ytestpred = model.predict(G_test)
            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
            rp_model = rp.Model(model_name='randomforest', model_info=[' '], training_metrics=scores,
                           testing_metrics=[accuracy_score(test_pheno, ytestpred),
                                            prfs(test_pheno, ytestpred.round(), average='micro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='macro')[-2],
                                            prfs(test_pheno, ytestpred.round(), average='weighted')[-2]])
            # sys.exit(0)
        #--------------------------------------------------------------------#
        else:
            print('Not a classifier covered in this script. See MLGWAS.rmd for the options...')
            sys.exit(0)
        #--------------------------------------------------------------------#

        metrics = ['accuracy', 'f1_micro', 'f1_macro', 'f1_weighted']
        # headers for output files
        dataset_name = str(prefix)
        file_prefix = dataset_name+'_'+classifiertype
        # time_string = datetime.datetime.now()
        rn = random.randint(1,200000)

        outdir = file_prefix+'_'+str(value)+'_Results_'+str(rn)
        os.mkdir(outdir)

        doc = rp.generate_pdf_with_name(file_prefix,rn)
        num_indivs = len(G_test) + len(G_train)
        n_snps = G_test.shape[1]

        # build the confusion matrix
        stats_var = prfs(test_pheno, ytestpred, average='macro')
        cnfmat = confusion_matrix(test_pheno, ytestpred)
        rp.fill_document(doc, [rp_model], metrics, cnfmat, dataset_name, num_indivs, n_snps, G_train.shape, G_test.shape,5)

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
        print(classification_report(test_pheno, ytestpred))

        # Move the files to their own directory
        absolute_prefix = file_prefix+'_'+str(rn)

        # f = open(absolute_prefix+'_parameters.txt','w')
        # f.write('Parameter List: \n')
        # for key in params_dict:
        #     f.write('{'+str(key)+','+str(params_dict[key])+'}\n')
        # f.close()
        # shutil.move(absolute_prefix+'_parameters.txt', outdir+'/'+absolute_prefix+'_parameters.txt')
        shutil.move(absolute_prefix+'.tex', outdir+'/'+absolute_prefix+'.tex')
        shutil.move(absolute_prefix+'.pdf', outdir+'/'+absolute_prefix+'.pdf')

        best_accuracy = (float(cnfmat[0][0]+cnfmat[1][1]))/(cnfmat[0][0]+cnfmat[0][1]+cnfmat[1][0]+cnfmat[1][1])
        best_cnfmat = cnfmat
        # Also store the corresponding case(postive) and control(negative) accuracies
        # TN/(TN+FP)
        control_accuracy = (float(cnfmat[0][0]))/(cnfmat[0][0]+cnfmat[0][1])
        # TP/(TP+FN)
        case_accuracy = (float(cnfmat[1][1]))/(cnfmat[1][0]+cnfmat[1][1])
        # Store corresponding f1 score
        f1score = round(stats_var[1],3)

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
    ################################ ML Stuff #################################
