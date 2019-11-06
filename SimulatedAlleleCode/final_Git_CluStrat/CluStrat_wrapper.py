from plinkio import plinkfile
import parseplink as pp
import numpy as np

import os, datetime, random, sys, shutil, itertools, math, argparse, time, subprocess

import normalize, getMH, CluStrat

def msg(name=None):
    return '''CluStrat_wrapper.py
         >> python3 CluStrat_wrapper.py --sim 1
         >> python3 CluStrat_wrapper.py --dir /testing/test_data
         *above the sample data is in the directory 'testing/' and the prefix to the PLINK format data is 'test_data'.
        '''

def parse_arguments():
    parser = argparse.ArgumentParser(usage=msg())

    parser.add_argument("-s", "--sim", dest='simulate', action='store', help="Indicate whether to use simulated data (--sim 1)."+
                        " If not using simulated data, indicate real dataset using '-d/--dir' flag.",
                        metavar="SIM")

    parser.add_argument("-d", "--dir", dest='realdata_directory', action='store', help="Put the path to the real dataset.",
                        metavar="DIR")

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

    ######################### Loading/Simulating Data #########################
    if args.simulate == "1":
        print('##################### Simulating data...\n')
        # Simulating data from all scenarios using Python 2.7
        COMMAND = "python data_simulate.py --model BN --prop 2 --pheno continuous"
        subprocess.call(COMMAND, shell=True)
        print('##################### Loading data...\n')
        # Grab a file from the simulated data directory
        file_handle = "sim_plinkfiles/BN/proportion2/continuous/"+str(random.choice(os.listdir("sim_plinkfiles/BN/proportion2/continuous/"))).split('.')[0]
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
    pvalue = 1e-7

    ############################# Normalize Data #############################
    print('##################### Normalizing data...\n')
    norm_time = time.time()
    normX,_ = normalize.norm(X,0)
    print('Normalizing time (secs): ', (time.time()-norm_time))
    print(' ')
    Y = pheno


    # Plotting PCA (top 2 PCs) of the data
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
    dele = [10,12]
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
