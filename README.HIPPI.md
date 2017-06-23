------------------------------------
Summary
------------------------------------

HIPPI stands for `Highly Accurate Protein Family Classification with Ensembles of HMMs`, and so is a method for the following problem:

Protein family identification:
- Input: A query sequence `q` and a set of alignments and trees on a protein family `F`
- Output: A score of `q` against all families `F`

HIPPI is a modification of UPP for scoring protein sequences against a protein family database.  HIPPI operates in two steps.  In the first step, HIPPI builds an ensemble of HMMs on an input protein family.  In the next step, HIPPI scores the query sequences against all the ensemble of HMMs for the protein family. By pipelining these step for all protein families, one can obtain the score of the query sequences against all families.  HIPPI is in active development, and scripts will soon be made available to pipeline these steps into a single command.

Developers: Nam Nguyen, Michael Nute, Siavash Mirarab, and Tandy Warnow.

###Data: 
Nguyen, Nam-phuong (2016): `HIPPI Dataset`. University of Illinois at Urbana-Champaign. https://doi.org/10.13012/B2IDB-6795126_V1

###Publication:
Nam Nguyen, Michael Nute, Siavash Mirarab, Keerthana Kumar, and Tandy Warnow. `HIPPI: Highly Accurate Protein Family Classification with Ensembles of profile Hidden Markov Models`. Accepted to RECOMB CG 2016 ().


### Note and Acknowledgment: 
- HIPPI bundles the following program into its distribution:
  1. hmmer: http://hmmer.janelia.org/    
- HIPPI uses the [Dendropy](http://pythonhosted.org/DendroPy/) package. 
- HIPPI uses some code from [SATe](http://phylo.bio.ku.edu/software/sate/sate.html).

-------------------------------------
Installation
-------------------------------------
This section details steps for installing and running HIPPI. We have run HIPPI on Linux and MAC. If you experience difficulty installing or running the software, please contact Nam Nguyen or Siavash Mirarab.

Requirements:
-------------------
Before installing the software you need to make sure the following programs are installed on your machine.  

1. Python: Version > 2.7. 
2. SEPP: Version > 3.0. 

Installation Steps:
-------------------
HIPPI is a part of the SEPP distribution package.  By installing SEPP, HIPPI is automatically installed. (see [SEPP readme] (https://github.com/smirarab/sepp/blob/master/README.SEPP.md)).  

Common Problems:
-------------------
1.  HIPPI requires SEPP to be installed.  If HIPPI is not running, first check to see if SEPP was installed correctly.

---------------------------------------------
Running HIPPI
---------------------------------------------
To run HIPPI, invoke the `run_ensembles.py` script from the `bin` sub-directory of the location in which you installed the Python packages. To see options for running the script, use the command:

`python <bin>/run_ensembles.py -h`

The general command for running HIPPI is:

`run_ensembles.py -a input_alignment -t input_tree -f input_query_sequences -A decomp_size - m amino -D 0.60`

where decomp_size is 10% of the original input alignment.  This will run HIPPI(10%,40%) as described in the main paper.  

The main output of HIPPI output_scores.csv.  This lists the score of the query sequences against all the ensemble of HMMs.  

We are currently building a pipeline script to streamline this process.

---------------------------------------------
Bugs and Errors
---------------------------------------------
HIPPI is under active research development at UIUC by the Warnow Lab (and especially with her former PhD students Siavash Mirarab and Nam Nguyen). Please report any errors or requests to Michael Nute (nute2@illinois.edu), Siavash Mirarab (smirarab@gmail.com) and Nam Nguyen (ndn006@eng.ucsd.edu).

