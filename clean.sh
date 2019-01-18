#!/bin/bash

rm -fr PRK_LDA_*
rm -fr PRK_QDA_*
rm -fr PRK_PolyReg_*
rm -fr PRK_Ridge_*
rm -fr PRK_Lasso_*
rm -fr PRK_Elastic_*
rm -fr PRK_SVM_*
rm -fr PRK_kSVM_*
rm -fr PRK_RidgeSVM_*
rm -fr PRK_RFESVM_*
rm -fr PRK_RandomForest_*
rm -f PKHCGSRRS_extracted_*

rm -fr SCZ_LDA_*
rm -fr SCZ_QDA_*
rm -fr SCZ_PolyReg_*
rm -fr SCZ_Ridge_*
rm -fr SCZ_Lasso_*
rm -fr SCZ_Elastic_*
rm -fr SCZ_SVM_*
rm -fr SCZ_kSVM_*
rm -fr SCZ_RidgeSVM_*
rm -fr SCZ_RFESVM_*
rm -fr SCZ_RandomForest_*
rm -f SZHCGSRRS_extracted_*

rm -f *.out 
rm -f *bestparameters.txt
rm -f *bestaccuracy.txt
rm -f *f1score.txt
rm -f *positiveaccuracy.txt
rm -f *negativeaccuracy.txt
#END
