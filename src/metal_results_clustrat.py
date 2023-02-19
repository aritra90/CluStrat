import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys, csv, os, math, subprocess, itertools, time, random, warnings, array, shutil, glob
warnings.filterwarnings(action="ignore")

def get_hits(df, pval_thresh, method_type, model):
    if method_type == 'metal':
        # hits = df.loc[df['P-value'] < pval_thresh]
        # hits = hits.sort_values('P-value')

        # checking spurious vs. causal
        # causal_df = pd.read_csv("/scratch/bell/mcburch/MOSAIC_SIM/causal_variants.txt", header = None)
        causal_df = pd.read_csv("/depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/"+str(model)+"_simdata/causal_variants.txt", header = None)
        thesnps = np.reshape(causal_df.values, causal_df.values.shape[0])
        true_hits = df.loc[df['MarkerName'].isin(thesnps)]
        hits = df.loc[df['P-value'] < pval_thresh]
        return true_hits, hits, hits.loc[hits['MarkerName'].isin(thesnps)]

        # return hits
    else:
        # hits = df.loc[df['P'] < pval_thresh]
        # hits = hits.sort_values('P')

        # checking spurious vs. causal
        # causal_df = pd.read_csv("/scratch/bell/mcburch/MOSAIC_SIM/causal_variants.txt", header = None)
        causal_df = pd.read_csv("/depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/"+str(model)+"_simdata/causal_variants.txt", header = None)
        thesnps = np.reshape(causal_df.values, causal_df.values.shape[0])
        true_hits = df.loc[df['ID'].isin(thesnps)]
        hits = df.loc[df['P'] < pval_thresh]
        return true_hits, hits, hits.loc[hits['ID'].isin(thesnps)]

        # return hits

if __name__ == '__main__':
    # pval_thresh = 1e-6
    metal_filenm = str(sys.argv[1])
    method_type = str(sys.argv[2])
    pval_thresh = float(sys.argv[3])
    model = str(sys.argv[4])
    df = pd.read_csv(metal_filenm, delim_whitespace=True)
    actual, hits, causal_hits = get_hits(df, pval_thresh, method_type, model)
    if method_type == 'metal':
        # print("clustrat-metal")
        hits = hits.sort_values('P-value')
        print(hits)
        # hits.MarkerName.to_csv("sig_snps_" + str(method_type), sep='\t', index = False, header = None)

        # causal
        print("causal:", causal_hits.shape[0])
        # spurious
        print("spurious:", hits.shape[0] - causal_hits.shape[0])

        with open('causal_clustrat', 'a+') as f:
            f.write(str(causal_hits.shape[0])+'\n')

        with open('spurious_clustrat', 'a+') as f:
            f.write(str(hits.shape[0] - causal_hits.shape[0])+'\n')
    else:
        # print("plink2")
        hits = hits.sort_values('P')
        print(hits)
        # hits.ID.to_csv("sig_snps_" + str(method_type), sep='\t', index = False, header = None)

        # causal
        print("causal:", causal_hits.shape[0])
        # spurious
        print("spurious:", hits.shape[0] - causal_hits.shape[0])

        with open('causal_plink2', 'a+') as f:
            f.write(str(causal_hits.shape[0])+'\n')

        with open('spurious_plink2', 'a+') as f:
            f.write(str(hits.shape[0] - causal_hits.shape[0])+'\n')
