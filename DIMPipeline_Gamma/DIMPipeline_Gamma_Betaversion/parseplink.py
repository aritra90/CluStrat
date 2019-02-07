
from plinkio import plinkfile
import numpy as np

def read_n_parse(plink_file):
	if not plink_file.one_locus_per_row( ):
     		print( "This script requires that snps are rows and samples columns." )
     		exit( 1 )

	sample_list = plink_file.get_samples( )
	locus_list = plink_file.get_loci( )

	stdout = 'The number of samples is ' + repr(len(sample_list)) + ' and '+ repr(len(locus_list)) + ' markers.'
	print(stdout)
	print('')

	l = [] 
	x = []
	for genotype in plink_file:
        	l.append( genotype )
	for sample in sample_list:
		x.append(sample.phenotype)

	pheno = np.asarray(x)
	pheno = pheno.astype(float)	
	G = np.asarray(l)
	G = G.reshape((len(locus_list),len(sample_list)))
	G = G.astype(float)
	G[G == 3] = np.nan
	G[G == 2] = 3.0
	G[G == 1] = 2.0
	G[G == 0] = 1.0
	G = np.transpose(G)
	means = np.nanmean(G,axis=0)
	G = G - means
	G[np.isnan(G)] = 0.0
	return G, pheno
