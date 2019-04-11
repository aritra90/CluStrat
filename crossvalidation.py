from sklearn.model_selection import train_test_split, KFold
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

def crossvalidation(trainingdata, train_pheno, classifiertype, crossval, elem):

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

        scores, model = cross_validation_custom(trainingdata, train_pheno,model_init, crossval, classifiertype)

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

        scores, model = cross_validation_custom(trainingdata, train_pheno,model_init, crossval, classifiertype)

        model = QuadraticDiscriminantAnalysis(reg_param=reg_param,tol=tol_QDA)

        params_dict = {'regularization':reg_param,'tolerance':tol_QDA}
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

        scores, model = cross_validation_custom(trainingdata, train_pheno,model_init, crossval, classifiertype, degree)

        params_dict = {'fit_intercept':fit_intercept,'degree':degree,'interaction_only':interaction_only,'include_bias':include_bias}
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

        scores, model = cross_validation_custom(trainingdata, train_pheno,model_init, crossval, classifiertype)

        params_dict = {'alpha':alpha,'fit_intercept':fit_intercept,'max_iter':max_iter,'tolerance':tol_Ridge,'solver':solver_Ridge}
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

        scores, model = cross_validation_custom(trainingdata, train_pheno,model_init, crossval, classifiertype)

        params_dict = {'alpha':alpha,'fit_intercept':fit_intercept,'max_iter':max_iter,'tolerance':tol_Lasso,'positive':positive,'selection':selection}
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

        scores, model = cross_validation_custom(trainingdata, train_pheno,model_init, crossval, classifiertype)

        params_dict = {'alpha':alpha,'l1_ratio':l1_ratio,'fit_intercept':fit_intercept,'max_iter':max_iter,'tolerance':tol_Elastic,'positive':positive,'selection':selection}
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

        scores, model = cross_validation_custom(trainingdata, train_pheno,model_init, crossval, classifiertype)

        params_dict = {'C':C,'kernel':kernel,'degree':degree,'gamma':gamma,
                        'coef0':coef0,'shrinking':shrinking,'probability':probability,
                        'tol':tol,'verbose':verbose,'max_iter':max_iter,'decision_function_shape':decision_function_shape,
                        'random_state':random_state}
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

            scores, model = cross_validation_custom(trainingdata, train_pheno,model_init, crossval, classifiertype)

            params_dict = {'C':C,'kernel':kernel,'degree':degree,'gamma':gamma,
                            'coef0':coef0,'shrinking':shrinking,'tol':tol_SVR,'verbose':verbose,
                            'max_iter':max_iter,'epsilon':epsilon}
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

        scores, model = cross_validation_custom(trainingdata, train_pheno,model_init, crossval, classifiertype)

        params_dict = {'C':C,'kernel':kernel,'degree':degree,'gamma':gamma,
                        'coef0':coef0,'shrinking':shrinking,'probability':probability,
                        'tol':tol,'verbose':verbose,'max_iter':max_iter,'decision_function_shape':decision_function_shape,
                        'random_state':random_state}
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

        scores, model = cross_validation_custom(trainingdata, train_pheno,model_init, crossval, classifiertype)

        params_dict = {'alpha':alpha,'fit_intercept':fit_intercept,'max_iter':max_iter,'tolerance':tol_RidgeSVM,'solver':solver_RidgeSVM}
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

        scores, model = cross_validation_custom(trainingdata, train_pheno,model_init, crossval, classifiertype)

        params_dict = {'C':C,'estimator_object':estimator_object,'n_features_to_select':n_features_to_select,'step':step,'verbose':verbose}
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
        # scores, model = cross_validation_custom(trainingdata, train_pheno, model_init, crossval)

        # Test based on BEST alphas

        model = RidgeCV(alphas=alphas,fit_intercept=fit_intercept,normalize=normalize,scoring=scoring,cv=cv,gcv_mode=gcv_mode,store_cv_values=store_cv_values)

        params_dict = {'alphas':alphas,'fit_intercept':fit_intercept,'normalize':normalize,'scoring':scoring,'cv':cv,'gcv_mode':gcv_mode,'store_cv_values':store_cv_values}

        model.fit(trainingdata, train_pheno)

        # Just so it doesnt break given the report filling setup (might be able to remove)
        scores = [model.score(trainingdata, train_pheno)]*4
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

        scores, model = cross_validation_custom(trainingdata, train_pheno,model_init, crossval, classifiertype)

        params_dict = {'n_estimators':n_estimators,'criterion':criterion,'max_depth':max_depth,'min_samples_split':min_samples_split,'min_samples_leaf':min_samples_leaf,'min_weight_fraction_leaf':min_weight_fraction_leaf,
                        'max_features':max_features,'max_leaf_nodes':max_leaf_nodes,'min_impurity_decrease':min_impurity_decrease,'min_impurity_split':min_impurity_split,'bootstrap':bootstrap,'oob_score':oob_score,
                        'n_jobs':n_jobs,'random_state':random_state,'verbose':verbose,'class_weight':class_weight}
    #--------------------------------------------------------------------#
    else:
        print('Not a classifier covered in this script. See MLGWAS.rmd for the options...')
        sys.exit(0)
    #--------------------------------------------------------------------#
    return scores, params_dict

def cross_validation_custom(X, y, model_init, crossval, classifiertype, degree=2, verbose=False):
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
            print('Degree: ', degree)
            poly = PolynomialFeatures(degree=degree,interaction_only=True)
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

if __name__ == '__main__':

    print("HI")
