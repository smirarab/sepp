[![SEPP_github_tests](https://github.com/smirarab/sepp/actions/workflows/sepp_tests.yml/badge.svg)](https://github.com/smirarab/sepp/actions/workflows/sepp_tests.yml) [![Coverage Status](https://coveralls.io/repos/github/smirarab/sepp/badge.svg?branch=master)](https://coveralls.io/github/smirarab/sepp?branch=master) [![Anaconda-Server Badge](https://anaconda.org/bioconda/sepp/badges/version.svg)](https://anaconda.org/bioconda/sepp)
------------------------------------
Summary
------------------------------------
This repository includes code for SEPP, TIPP, UPP, HIPPI.  The three methods use ensembles of Hidden Markov Models (HMMs) in different ways, each focusing on a different problem.

Each of these related tools has its own README file.

[README.SEPP.md](README.SEPP.md)
* **SEPP** stands for "SATe-enabled Phylogenetic Placement", and addresses the problem of phylogenetic placement of short reads into reference alignments and trees.

[README.UPP.md](README.UPP.md)
* **UPP** stands for "Ultra-large alignments using Phylogeny-aware Profiles", and addresses the problem of alignment of very large datasets, potentially containing fragmentary data. UPP can align datasets with up to 1,000,000 sequences.

[README.HIPPI.md](README.HIPPI.md)
* **HIPPI** stands for "Highly Accurate Protein Family Classification with Ensembles of HMMs", and addresses the problem of classifying query sequences to protein families.

[README.TIPP.md](https://github.com/TeraTrees/TIPP/)
* **TIPP** stands for "Taxonomic Identification and Phylogenetic Profiling", and addresses the problem of taxonomic identification and abundance profiling of metagenomic data. We have moved TIPP to be a separate package from SEPP. TIPP package can be accessed [here](https://github.com/TeraTrees/TIPP/).

**NOTE:** All these programs heavily rely on [HMMER](http://hmmer.org/). Please cite HMMER when citing these tools as well and mention the version of the HMMER used. 

---------------------------------------------
Bugs and Errors
---------------------------------------------
SEPP, TIPP, UPP, HIPPI are under active research development at UIUC by the Warnow Lab and former student Siavash Mirarab (now at UCSD). Please report any errors on the GitHub issues page. 

