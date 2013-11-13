------------------------------------
Summary
------------------------------------
This repository includes code for SEPP and TIPP.

Each of these related tools have their own README files [README.SEPP.md](README.SEPP.md) and [README.TIPP.md](README.TIPP.md).

* **SEPP** stands for "SATe-enabled phylogenetic placement", and addresses the problem of phylogenetic placement of short reads into reference alignment and trees. 
* **UPP** stands for "Ultra-large alignment using SEPP", and addresses the problem of aligning ultra-large alignments, both full length and fragmentary sequences.  UPP randomly samples a backbone set and estimates an alignment and tree on the backbone set.  The remaining sequences are aligned to the backbone set.
* **TIPP** stands for "Taxon Identification using Phylogenetic Placement". TIPP classifies unknown fragmentary sequences into a given taxonomy. TIPP uses phylogenetic placement on a reference alignment/tree for classification. It requires that input fragments come from the same gene as the full length sequences in the reference alignment. 

---------------------------------------------
Bugs and Errors
---------------------------------------------
SEPP, UPP, and TIPP are under active research development at UTCS by the Warnow Lab (and especially with her PhD students Siavash Mirarab and Nam Nguyen). Please report any errors to Siavash Mirarab (smirarab@gmail.com) and Nam Nguyen (namphuon@cs.utexas.edu).

