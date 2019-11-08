<img src="src/purduecs.jpg" align="right" />

# CluStrat: 

This software performs agglomerative hierarchical clustering with the 
Mahalanobis distance based Genetic Relationship Matrix (GRM) representing the 
population-level covariance (LD) for the genetic markers. 

## Description

Genome-wide association studies (GWAS) have been extensively used to 
estimate the signed effects of trait-associated alleles and compute 
polygenic risk scores. Recently, it has been made evident that more 
rigorous and sophisticated methods for standard population structure correction 
are needed . Here, we provide a correction technique for complex population structure 
while leveraging the linkage disequilibrium (LD) induced distances between individuals. We implement CluStrat, which performs agglomerative hierarchical clustering using the Mahalanobis distance based Genetic Relationship Matrix (GRM) representing the population-level covariance (LD) for the SNPs. Here, we provide a comprehensive guide to stratification and subsequent disorder trait prediction or estimation utilizing the underlying LD structure of the genotypes.

## Getting Started

### Dependencies

* This software was developed using Python 3 (Python 2 for simulating datasets). 
* There are many imports for the software and you just need to make sure your pip package has the necessary packages installed. For example:
```
python3 -m pip install plinkio
```

### Installing

* Files can be installed from this repository.

## Executing programs
* The main file to run this software is CluStrat_wrapper.py. This code allows you to indicate whether you want to run CluStrat on simulated data or on real data (how to run below). The simulated data is fixed to one scenario by default but can be altered.  
```
python3 CluStrat_wrapper.py --sim 1
```
```
python3 CluStrat_wrapper.py --dir example/test_data 
```
```
python3 CluStrat_wrapper.py --help
```
<!--
* Another important file that can be run is StratCompare.py. This script can be run to compare Armitage Trend CHISQ, EigenStrat, Gemma and Emmax methods with CluStrat on simulated data. The paths to the various software packages need to be edited accordingly to where they are on your machines.
```
python StratCompare.py 1
```
* The last file to run is geneAnnot.r. This Rscript can be run taking the ouput file of CluStrat.py and the number of desired top annotations as input. The code uses biomaRt to extract gene annotations for the significant SNPs found by CluStrat. 
```
Rscript geneAnnot.r CluStrat_signficantSNPs_dele0.txt 40
```
-->

## Output 
* The output of CluStrat_wrapper.py are the chromosome number, SNP rsIDs and p-values from ridge regression. The format is the following: 
```
chrom SNPs p-values
3 rs2875479 1.2204460492567813e-16
6 rs3456713 4.220446565656313e-15
5 rs1987654 3.9854337898655673e-13
...
```

## Notes

* The simulation code (data_simulate.py) is ran using Python 2. When running with Python 3, a segmentation fault occurs when trying to save the simulated data in PLINK format using the libplinkio library. We are currently working on making this fix as Python 2 will be deprecated soon. 

  Reference: https://github.com/mfranberg/libplinkio

* Another note when running CluStrat is to adjust the clustering depth based on the dendogram to get an appropriate number of desired clusters. During the execution, the dendogram plot is saved so you can halt the execution to view the plot, adjust the depth accordingly and re-run the code.

## Authors and Correspondence 

Aritra Bose (a dot bose @ ibm dot com)

Myson Burch (mcburch @ purdue dot edu)

<!--

```
* The output of StratCompare.py is ... The format is the following:
```
```
* The output from running geneAnnot.r are the annotations from biomaRt. The format is the following:
```
    refsnp_id ensembl_gene_stable_id associated_gene
1   rs6699993
2  rs12049279
3   rs6689517
4   rs6427623
5  rs12096958
6   rs2794867                        RPL13AP11,CNTN2
7   rs6692892        ENSG00000143353
8  rs10449246
9   rs6696837        ENSG00000198626
10  rs6696837                LRG_402
...
```

## Version History
* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release
## License
This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details
## Acknowledgments
Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
-->

## Acknowledgments
