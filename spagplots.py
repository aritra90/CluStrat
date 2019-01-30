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
	#threads = [1,2,3,4,5,6,7,8,9,10]
	#plt.style.use('fivethirtyeight')
	plt.style.use('bmh')
	fig, ax = plt.subplots(1)
	headers = df.dtypes.index
	# multiple line plot; Put the entry (1,1) of the table as values. 
	#A sample table looks like: 
	#################################################################
	# Values	LDA	QDA	PolyReg(k=3) PolyReg(k=8) SVM(poly=3)	SVM-linear	SVM-rbf	RFESVM	etc...
	# 100
	# 400
	# 800
	# ...
	# ...
	for column in df.drop('Values', axis=1):
		l1 = ax.plot(df['Values'], df[column], '-', marker='*', ms=8,linewidth=3.5, alpha=0.75)
	#for column in df.drop('PCs', axis=1):
	#	l2 = ax.plot(df['PCs'], df[column], '.-', marker='^', ms=8,linewidth=3.5, alpha=0.75)	
		
	#print(ax.get_yticks())
	#ax.set_ylim([1,3])
	savename = fname+"_spagplot.png"
	ax.set_yticks(ax.get_yticks()[1::2])
	plt.xticks(threads,fontsize=8)
	plt.yticks(fontsize=8)\
	plt.legend( prop={'size': 10})
	plt.xlabel('# of threads', fontsize=10)
	plt.ylabel('Speed-up', fontsize=10)
	fig.savefig(savename, bbox_inches='tight', dpi=800)
	