import numpy as np
import pandas as pd
import random as rng
import sys, itertools, glob, warnings, scipy, os, math, time, timeit, subprocess, re
import matplotlib as pyplot
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import optimize
warnings.filterwarnings("ignore")
from matplotlib.patches import Patch
sns.set(style='darkgrid')
sns.set_palette("pastel")
pyplot.rcParams['figure.dpi']=600
pyplot.rcParams.update({'font.size': 10})
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title
from scipy.interpolate import interp1d


maxclust = sys.argv[1]
mhdist = sys.argv[2]
model = sys.argv[3]




# causal barplots
data = []
methods = []
causal_fns = ['causal_clustrat', 'causal_plink2', 'causal_armitage', 'causal_eigenstrat']
for elem in causal_fns:
    method = elem.split('_')[1]
    stuff = pd.read_csv(elem, header = None).values.flatten()
    data.append(stuff)
    methods.append(method)
data = pd.DataFrame(data).T
data.columns = causal_fns
ax = sns.boxplot(data=data)
# ax = sns.stripplot(data=data, color="orange", jitter=0.2, size=2.5)
ax.set(xticklabels=[])
ax.set(ylabel="number of causal assocs.")
legend_elements = [Patch(facecolor='#d0bbff', edgecolor='dimgrey',label='Color Patch'),
                    Patch(facecolor='#ffb482', edgecolor='dimgrey',label='Color Patch'),
                    Patch(facecolor='#8de5a1', edgecolor='dimgrey',label='Color Patch'),
                    Patch(facecolor='#ff9f9b', edgecolor='dimgrey',label='Color Patch'),]
plt.legend(handles = legend_elements, labels = methods)
plt.savefig("compare_others_"+str(model)+"_mhdist"+str(mhdist)+"_clust"+str(maxclust)+"_causal.png", bbox_inches = 'tight')



# sns.color_palette("pastel", 4)
# sns.color_palette("pastel").as_hex()[1]
# sns.color_palette("pastel").as_hex()[2]
# sns.color_palette("pastel").as_hex()[3]
# sns.color_palette("pastel").as_hex()[4]

plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

# spurious barplots
data = []
methods = []
spurions_fns = ['spurious_clustrat', 'spurious_plink2', 'spurious_armitage', 'spurious_eigenstrat']
for elem in spurions_fns:
    method = elem.split('_')[1]
    stuff = pd.read_csv(elem, header = None).values.flatten()
    data.append(stuff)
    methods.append(method)
data = pd.DataFrame(data).T
data.columns = spurions_fns
ax = sns.boxplot(data=data)
# ax = sns.stripplot(data=data, color="orange", jitter=0.2, size=2.5)
ax.set(xticklabels=[])
legend_elements = [Patch(facecolor='#d0bbff', edgecolor='dimgrey',label='Color Patch'),
                    Patch(facecolor='#ffb482', edgecolor='dimgrey',label='Color Patch'),
                    Patch(facecolor='#8de5a1', edgecolor='dimgrey',label='Color Patch'),
                    Patch(facecolor='#ff9f9b', edgecolor='dimgrey',label='Color Patch'),]
plt.legend(handles = legend_elements, labels = methods)
plt.savefig("compare_others_"+str(model)+"_mhdist"+str(mhdist)+"_clust"+str(maxclust)+"_spurious.png", bbox_inches = 'tight')








# end of code
