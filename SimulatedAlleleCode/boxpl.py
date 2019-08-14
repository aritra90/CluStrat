import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import seaborn as sns
import pandas as pd 
import os

#Run boxplot with structured input like PSD.txt attached. 

def parse_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--file1", dest='ht', action='store', help="Put the path to the stats",metavar="PCF")
	args = parser.parse_args()
	return args.ht
	
rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
	
if __name__ == '__main__':
	ht = parse_arguments()
	ht_pd = pd.read_csv(ht, sep='\t')
	print(ht_pd)
	base = os.path.basename(ht)
	fname = os.path.splitext(base)[0]
	sns.set(style="ticks", palette="deep")
	sns.set_style("darkgrid")
	#tips = sns.load_dataset("tips")
	fig = sns.boxplot(x="Proportion", y="Spurious",
            hue="Model",
            data=ht_pd)
	sns.despine(offset=10, trim=True)
	plt.legend(loc="upper right", fontsize = "x-small")
	fig.get_figure().savefig(fname+"_spurious.png", bbox_inches='tight', dpi=600)

