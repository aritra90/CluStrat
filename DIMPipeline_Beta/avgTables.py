import argparse
import numpy as np
import pandas as pd
import os

def parse_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-exp", "--experiments", dest='numExps', action='store', help="Put the number of experiments ran...",metavar="EXP")
	args = parser.parse_args()
	return args.numExps

if __name__ == '__main__':
	numExps = int(parse_arguments())


	############# Case Accuracy Table #############
	caseFileArr = [""]*numExps
	casedfArr = [0]*numExps
	for i in range(0,len(caseFileArr)):
		caseFileArr[i] = "Experiments/Exp_"+str(i+1)+"_PRK_caseaccuracy.csv"
		casedfArr[i] = pd.read_csv(caseFileArr[i], sep=r'\t', engine='python')
		# print(casedfArr[i])
		# print("")

	caseAvgdf = pd.concat(casedfArr[:]).groupby(level=0).mean()
	# print(caseAvgdf)
	caseStdDevdf = pd.concat(casedfArr[:]).groupby(level=0).std()
	# print(caseStdDevdf)

	caseAvgdf.to_csv("PRK_caseaccuracy.csv", sep='\t',index=False)
	caseStdDevdf.to_csv("PRK_caseaccuracy_stddev.csv", sep='\t')


	############# Control Accuracy Table #############
	controlFileArr = [""]*numExps
	controldfArr = [0]*numExps
	for i in range(0,len(controlFileArr)):
		controlFileArr[i] = "Experiments/Exp_"+str(i+1)+"_PRK_controlaccuracy.csv"
		controldfArr[i] = pd.read_csv(controlFileArr[i], sep=r'\t', engine='python')
		# print(controldfArr[i])
		# print("")

	controlAvgdf = pd.concat(controldfArr[:]).groupby(level=0).mean()
	# print(controlAvgdf)
	controlStdDevdf = pd.concat(controldfArr[:]).groupby(level=0).std()
	# print(controlStdDevdf)

	controlAvgdf.to_csv("PRK_controlaccuracy.csv", sep='\t',index=False)
	controlStdDevdf.to_csv("PRK_controlaccuracy_stddev.csv", sep='\t')


	############# Training Accuracy Table #############
	trainingFileArr = [""]*numExps
	trainingdfArr = [0]*numExps
	for i in range(0,len(trainingFileArr)):
		trainingFileArr[i] = "Experiments/Exp_"+str(i+1)+"_PRK_trainingaccuracy.csv"
		trainingdfArr[i] = pd.read_csv(trainingFileArr[i], sep=r'\t', engine='python')
		# print(trainingdfArr[i])
		# print("")

	trainingAvgdf = pd.concat(trainingdfArr[:]).groupby(level=0).mean()
	# print(trainingAvgdf)
	trainingStdDevdf = pd.concat(trainingdfArr[:]).groupby(level=0).std()
	# print(trainingStdDevdf)

	trainingAvgdf.to_csv("PRK_trainingaccuracy.csv", sep='\t',index=False)
	trainingStdDevdf.to_csv("PRK_trainingaccuracy_stddev.csv", sep='\t')


	############# Testing Accuracy Table #############
	testingFileArr = [""]*numExps
	testingdfArr = [0]*numExps
	for i in range(0,len(testingFileArr)):
		testingFileArr[i] = "Experiments/Exp_"+str(i+1)+"_PRK_testingaccuracy.csv"
		testingdfArr[i] = pd.read_csv(testingFileArr[i], sep=r'\t', engine='python')
		# print(testingdfArr[i])
		# print("")

	testingAvgdf = pd.concat(testingdfArr[:]).groupby(level=0).mean()
	# print(testingAvgdf)
	testingStdDevdf = pd.concat(testingdfArr[:]).groupby(level=0).std()
	# print(testingStdDevdf)

	testingAvgdf.to_csv("PRK_testingaccuracy.csv", sep='\t',index=False)
	testingStdDevdf.to_csv("PRK_testingaccuracy_stddev.csv", sep='\t')


	############# F1 score Table #############
	f1scoreFileArr = [""]*numExps
	f1scoredfArr = [0]*numExps
	for i in range(0,len(f1scoreFileArr)):
		f1scoreFileArr[i] = "Experiments/Exp_"+str(i+1)+"_PRK_f1score.csv"
		f1scoredfArr[i] = pd.read_csv(f1scoreFileArr[i], sep=r'\t', engine='python')
		# print(f1scoredfArr[i])
		# print("")

	f1scoreAvgdf = pd.concat(f1scoredfArr[:]).groupby(level=0).mean()
	# print(f1scoreAvgdf)
	f1scoreStdDevdf = pd.concat(f1scoredfArr[:]).groupby(level=0).std()
	# print(f1scoreStdDevdf)

	f1scoreAvgdf.to_csv("PRK_f1score.csv", sep='\t',index=False)
	f1scoreStdDevdf.to_csv("PRK_f1score_stddev.csv", sep='\t')
