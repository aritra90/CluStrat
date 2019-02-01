###################################################################################
#####	Author: Aritra Bose
#####	Date: 26/01/2019
##### Spaghetti Plots for sensitivity and specificity for different classifiers
###################################################################################

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import math
import pandas as pd
import os

def parse_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--file", dest='datafile', action='store', help="Put the path to the PCs",metavar="PCF")
	args = parser.parse_args()
	return args.datafile

if __name__ == '__main__':
	datahandle = parse_arguments()
	base = os.path.basename(datahandle)
	fname = os.path.splitext(base)[0]
	#change the delimiter as you will
	df = pd.read_csv(datahandle, sep='\t')
	xvals = [100,400,800,1200,1600,2000]

	plt.style.use('bmh')
	fig, ax = plt.subplots(1)
	headers = df.dtypes.index

	color_arr = ['b','g','r','c','m','y','k','darkorange','goldenrod','gray','purple']

	# multiple line plot; Put the entry (1,1) of the table as values.
	#A sample table looks like:
	#################################################################
	# Values	LDA	QDA	PolyReg(k=3) PolyReg(k=8) SVM(poly=3)	SVM-linear	SVM-rbf	RFESVM	etc...
	# 100
	# 400
	# 800
	# ...
	# ...
	ctr = 0
	for column in df.drop('Values', axis=1):
		l1 = ax.plot(df['Values'], df[column], '-', marker='*', ms=8,linewidth=3.5, alpha=0.75, label = column,color=color_arr[ctr])
		ctr = ctr + 1

	savename = fname+"_spagplot.png"
	ax.set_yticks(ax.get_yticks()[1::2])
	ax.legend(prop={'size': 10},loc='center left', bbox_to_anchor=(1, 0.5))
	plt.xticks(xvals,fontsize=8)
	plt.yticks(fontsize=8)\
	#plt.legend( prop={'size': 10})
	plt.xlabel('# of SNPs', fontsize=10)
	plt.ylabel('Accuracy', fontsize=10)
	fig.suptitle(str(fname), fontsize=20)
	fig.savefig(savename, bbox_inches='tight', dpi=800)
