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

import os, datetime, random, sys

def classify(trainingdata, train_pheno, testingdata, best_param_arr, classifiertype):

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

        model.fit(trainingdata, train_pheno)
        # test the model on your holdout test data
        ytestpred = model.predict(testingdata)
        ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
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

        model.fit(trainingdata, train_pheno)
        # test the model on your holdout test data
        ytestpred = model.predict(testingdata)
        ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
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

        trainingdata_ply = poly.fit_transform(trainingdata)
        testingdata_ply = poly.fit_transform(testingdata)
        model.fit(trainingdata_ply, train_pheno)
        # test the model on your holdout test data
        ytestpred = model.predict(testingdata_ply)
        # discretization of the continuous predictions
        ytestpred = np.asarray([0 if counter < 1 else 1 for counter in ytestpred])
        ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
        # print(ytestpred.round())
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

        model.fit(trainingdata, train_pheno)
        # test the model on your holdout test data
        ytestpred = model.predict(testingdata)

        ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
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

        model.fit(trainingdata, train_pheno)
        # test the model on your holdout test data
        ytestpred = model.predict(testingdata)
        ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
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

        model.fit(trainingdata, train_pheno)
        # test the model on your holdout test data
        ytestpred = model.predict(testingdata)
        ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
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

        model.fit(trainingdata, train_pheno)
        # test the model on your holdout test data
        ytestpred = model.predict(testingdata)
        ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
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

            model.fit(trainingdata, train_pheno)
            # test the model on your holdout test data
            ytestpred = model.predict(testingdata)
            ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
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

        model.fit(trainingdata, train_pheno)
        # test the model on your holdout test data
        ytestpred = model.predict(testingdata)
        ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
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

        model.fit(trainingdata, train_pheno)
        # test the model on your holdout test data
        ytestpred = model.predict(testingdata)
        ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
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

        model.fit(trainingdata, train_pheno)
        # test the model on your holdout test data
        ytestpred = model.predict(testingdata)
        ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
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
        # scores, model = cross_validation_custom(trainingdata, train_pheno, model_init, crossval)

        # Test based on BEST alphas

        model = RidgeCV(alphas=alphas,fit_intercept=fit_intercept,normalize=normalize,scoring=scoring,cv=cv,gcv_mode=gcv_mode,store_cv_values=store_cv_values)



        params_dict = {'alphas':alphas,'fit_intercept':fit_intercept,'normalize':normalize,'scoring':scoring,'cv':cv,'gcv_mode':gcv_mode,'store_cv_values':store_cv_values}

        model.fit(trainingdata, train_pheno)

        # Just so it doesnt break given the report filling setup (might be able to remove)
        scores = [model.score(trainingdata, train_pheno)]*4

        # test the model on your holdout test data
        ytestpred = model.predict(testingdata)
        ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
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

        model.fit(trainingdata, train_pheno)
        # test the model on your holdout test data
        ytestpred = model.predict(testingdata)
        ytestpred = np.where(ytestpred.round() >= 1, 1, 0)
    #--------------------------------------------------------------------#
    else:
        print('Not a classifier covered in this script. See MLGWAS.rmd for the options...')
        sys.exit(0)
    #--------------------------------------------------------------------#

    return ytestpred

if __name__ == '__main__':

    print("HI")
