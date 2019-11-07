<img src="icon.png" align="right" />

# CluStrat: 

This software performs agglomerative hierarchical clustering with the Mahalanobis distance based Genetic Relationship Matrix (GRM) representing the population-level covariance (LD) for the genetic markers. 

## Description

Genome-wide association studies (GWAS) have been extensively used to estimate the signed effects of trait-associated alleles and compute polygenic risk scores. Recently, it has been made evident that more rigorous and sophisticated methods for standard population structure correction are needed . Here, we provide a correction technique for complex population structure while leveraging the linkage disequilibrium (LD) induced distances between individuals. We implement CluStrat, which performs agglomerative hierarchical clustering using the Mahalanobis distance based Genetic Relationship Matrix (GRM) representing the population-level covariance (LD) for the SNPs. Here, we provide a comprehensive guide to stratification and subsequent disorder trait prediction or estimation utilizing the underlying LD structure of the genotypes.

## Getting Started

### Dependencies

* This software was developed using Python 3 (Python 2 for simulating datasets). 
* There are many imports for the software and you just need to make sure your pip package has the necessary packages installed. For example:
```
python3 -m pip install plinkio
```

### Installing

* Files can be installed from this repository.

### Executing program
```
python3 CluStrat_wrapper.py --sim 1
```
```
python3 CluStrat_wrapper.py --dir test/test_data 
```

## Help

Any advise for common problems or issues.
```
python3 CluStrat_wrapper.py --help
```

## Notes

## Authors 

Aritra Bose (email@site.com)

Myson Burch (email@site.com)

<!---
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
