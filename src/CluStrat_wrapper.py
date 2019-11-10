from plinkio import plinkfile
import parseplink as pp
import numpy as np

import os, datetime, random, sys, shutil, itertools, math, argparse, time, subprocess
from scipy.sparse.linalg import svds
import normalize, getMH, CluStrat
allmodels = ['BN', 'PSD', 'TGP']
def msg(name=None):
    return '''CluStrat_wrapper.py
         >> python3 CluStrat_wrapper.py --sim 1 --prop 10,20,70 --trait 1 --model BN --ver 0 --plot 0 --pval 0.0001 --numclust 10 --size 1000,1000 
         >> python3 CluStrat_wrapper.py --dir /testing/test_data --pval 0.0001 --numclust 10
         *above the sample data is in the directory 'testing/' and the prefix to the PLINK format data is 'test_data'.
        '''

def parse_arguments():
    parser = argparse.ArgumentParser(usage=msg())

    parser.add_argument("-s", "--sim", dest='simulate', action='store', help="Indicate whether to use simulated data (--sim 1)."+
                        " If not using simulated data, indicate real dataset using '-d/--dir' flag.",
                        metavar="SIM")

    parser.add_argument("-d", "--dir", dest='realdata_directory', action='store', help="Put the path to the real dataset.",
                        metavar="DIR")

    parser.add_argument("-p", "--pval", dest='pvalue', action='store', help="Enter the desired p-value threshold",
                        metavar="PV")

    parser.add_argument("-pr", "--prop", dest='prop', action='store', help="Enter comma separated proportion of variation",
                        metavar="PROP")

    parser.add_argument("-tf", "--trait", dest='trait_flag', action='store', help="Enter the trait flag -- 1 for binary, 0 for continuous",
                        metavar="TR")

    parser.add_argument("-m", "--model", dest='model', action='store', help="Enter the desired simulation model",
                        metavar="MODEL")

    parser.add_argument("-sz", "--size", dest='size', action='store', help="Enter the desired simulation matrix dimensions (individuals by SNPs)",
                        metavar="SIZE")

    parser.add_argument("-v", "--ver", dest='verbose', action='store', help="Set for verbose output with timing profile",
                        metavar="VER")

    parser.add_argument("-pf", "--plot", dest='plot_flag', action='store', help="flag for plotting",
                        metavar="PFLAG")

    parser.add_argument("-nc", "--numclust", dest='numclust', action='store', help="Enter comma separated number of clusters to train with",
                        metavar="PV")
    args = parser.parse_args()

    return args

def read_handlers(datahandle):
    # print('Real dataset: ',datahandle)
    plink_file = plinkfile.open(datahandle)
    geno, pheno = pp.read_n_parse(plink_file)

    return geno, pheno

if __name__ == '__main__':
    begin_time = time.time()

    ########################### Parsing Input Args ###########################
    print('##################### Parsing in arguments...\n')
    args = parse_arguments()

    if args.pvalue:
        try:
            pvalue = float(args.pvalue)
        except ValueError:
            print("Usage: Pvalue flag should be a float (0.001 or 1e-3)")
    else:
        print("Usage: Pvalue flag should be a float (0.001 or 1e-3)")
        sys.exit(1)

    if args.numclust:
        try:
            dele = [int(args.numclust)]
        except ValueError:
            print("Usage: Number of clusters flag must be an integer")
            sys.exit(1)
    else:
        print("Usage: Number of clusters flag must be an integer")
        sys.exit(1)

    ######################### Loading/Simulating Data #########################
    if args.simulate == "1":
        print('##################### Simulating data...\n')
        # Simulating data from all scenarios using Python 2.7
        if args.prop:
            v = args.prop.split(',')
            try:
                v_set = [int(i) for i in v]

                if len(v) != 3 or np.sum(v_set) != 100:
                    print("Usage: -pr 10,20,70 means 10% genetic, 20% environmental and 70% noise variance.")
                    sys.exit(1)

            except ValueError:
                print("Usage: -pr 10,20,70 means 10% genetic, 20% environmental and 70% noise variance.")
                sys.exit(1)
        else:
            print("Usage: -pr 10,20,70 means 10% genetic, 20% environmental and 70% noise variance.")
            sys.exit(1)

        if args.size:
            dims = args.size.split(',')
            try:
                sim_shape = [int(i) for i in dims]

                if len(dims) != 2:
                    print("Usage: -sz 1000,10000 means 1k individuals and 10k SNPs")
                    sys.exit(1)

            except ValueError:
                print("Usage: -sz 1000,10000 means 1k individuals and 10k SNPs")
                sys.exit(1)
        else:
            print("Usage: -sz 1000,10000 means 1k individuals and 10k SNPs")
            sys.exit(1)

        if args.trait_flag != "1" and args.trait_flag != "0":
            print("Usage: -tf 0 (for continuous) or 1 (for binary)")
            sys.exit(1)
        elif args.model not in allmodels:
            print("Usage: Model flag should either BN (Balding-Nichols) | PSD (Pritchard-Stephens-Donnelly) | TGP (1000 Genomes Project)")
            sys.exit(1)
        elif args.verbose != "0" and args.verbose != "1":
            print("Usage: Verbose should either be 0 or 1")
            sys.exit(1)
        elif args.plot_flag != "0" and args.plot_flag != "1":
            print("Usage: Plot flag can either be 0 or 1")
            sys.exit(1)

        str_cmd = "python data_simulate.py --model " + str(args.model) + " --prop " + str(args.prop) + " --pheno " + str(args.trait_flag) + " --size " + str(args.size)
        # COMMAND = "python data_simulate.py --model BN --prop 20,10,70 --pheno continuous"
        subprocess.call(str_cmd, shell=True)

        if args.verbose == "1":
           print('##################### Loading data...\n')
        # Grab a file from the simulated data directory
        file_handle = "simfile_"+str(args.model)+"_"+str(v[0])+"_"+str(v[1])+"_"+str(v[2])+"_"+str(args.trait_flag)
        print('Given dataset: '+file_handle)
        load_time = time.time()
        X, pheno = read_handlers(file_handle)
        print("Loaded genotype matrix of dimension ", X.shape)
        print('Loading time (secs): ', (time.time()-load_time))
        print(' ')
        # pass

    elif args.realdata_directory:
        print('##################### Loading data...\n')
        print('Given dataset: '+args.realdata_directory)
        load_time = time.time()
        file_handle = str(args.realdata_directory)
        X, pheno = read_handlers(file_handle)
        print("Loaded genotype matrix of dimension ", X.shape)
        print('Loading time (secs): ', (time.time()-load_time))
        print(' ')
        # pass

    else:
        print("Indicate whether to simulate data or to use real data to run CluStrat...")
        sys.exit(0)


    # Set the number of individuals, SNPs and pvalue
    m = X.shape[0]
    n = X.shape[1]

    ############################# Normalize Data #############################
    print('##################### Normalizing data...\n')
    norm_time = time.time()
    normX,_ = normalize.norm(X,0)
    print('Normalizing time (secs): ', (time.time()-norm_time))
    print(' ')
    Y = pheno

    # Plotting PCA (top 2 PCs) of the data
    if args.plot_flag == 1:
        plot_flag = 1
    else:
        plot_flag = 0

    if (plot_flag == 1):
        svd_time = time.time()
        temp = np.matmul(normX,normX.T)
        U, Sig, _ = svds(temp, k =10)
        print('SVD time (mins): ', (time.time()-svd_time)/60.0)
        print(' ')
        print(np.sqrt(Sig))
        # PCA PLOTS
        # grab indices for cases and controls
        idxcase = np.where(Y == 1)[0]
        idxcontr = np.where(Y == 0)[0]

        plt.style.use('seaborn-darkgrid')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        scale = 200.0 * np.random.rand(750)
        print('# Cases: '+str(len(idxcase)))
        print('# Controls: '+str(len(idxcontr)))
        ax.scatter(U[idxcase, 0], U[idxcase, 1], marker='o', color='red', s=10, label='case')
        ax.scatter(U[idxcontr, 0], U[idxcontr, 1], marker='*', color='blue', s=10, label='control')
        ax.legend()
        ax.set_xlabel('PC 1')
        ax.set_ylabel('PC 2')
        plt.savefig('CluStrat_top2PCs.png', bbox_inches='tight', dpi=600)


    ######################### Compute Distance Matrix #########################
    print('##################### Getting distance matrix...\n')
    dist_time = time.time()
    # 1000 choose 2 pairs (calculating distance between each pair)
    #D = squareform(pdist(normX))
    D = getMH.MH(normX)
    print('Calculating distance matrix time (mins): ', (time.time()-dist_time)/60.0)
    print(' ')
    # print(D)
    # dele = [3, 5]
    # variable to control different numbers of clusters
    # dele = [10,12] # = [args.numclust]?
    d = 2

    ######################### Run CluStrat Algorithm #########################
    print('##################### Running CluStrat...')

    # Grabbing SNP rsIDs and chromosome from the data
    SNPids = []
    chromids = []
    with open(file_handle+'.bim', 'r') as f:
        for line in f:
            elems = line.split('\t')
            # print(elems)
            chromids.append(elems[0])
            SNPids.append(elems[1])

    clu_time = time.time()
    # SP, CS, clustcount, pvals, sigidx = CluStrat.cluster(X, D, d, Y, pvalue, dele, [SNPids, chromids])
    SP, CS, clustcount = CluStrat.cluster(X, D, d, Y, pvalue, dele, 0, [SNPids, chromids])
    print('CluStrat time (mins): ', (time.time()-clu_time)/60.0)

    CS = list(CS)
    SP = list(SP)

    SPmax = max(SP)
    SP3 = SPmax
    #CS3 = CS[maxidx]
    print('Total time (mins): ', (time.time() - begin_time)/60.0)
