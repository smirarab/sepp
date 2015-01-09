------------------------------------
Summary
------------------------------------

UPP stands for `Ultra-large alignments using Phylogeny-aware Profiles`, and so is a method for the following problem:

Alignment:
- Input: A set of unaligned sequences `S`
- Output: An alignment `A` on `S`

UPP is a modification of SEPP for performing alignments of ultra-large and fragmentary datasets.  UPP operates in four steps.  In the first step, UPP partitions set `S` into a backbone set and a query set and computes an alignment and tree on the backbone set.  In the next step, UPP decomposes the backbone alignment into ensemble of profile Hidden Markov models (HMM).  The third step, UPP searches for the best alignment of the query sequence to each HMM.  The final step inserts the query sequence into the backbone alignment using the best scoring HMM.  Our study shows that UPP results in accurate alignments and ML trees estimated on the alignments, is robust to datasets containing both fragmentary and full-length sequences, and is fast enough to produce an alignment on 1,000,000 sequences in two days.

Developers: Tandy Warnow, Nam Nguyen, and Siavash Mirarab

###Publication:
Nam Nguyen, Siavash Mirarab, Keerthana Kumar, and Tandy Warnow. `Ultra-large alignments using ensembles of Hidden Markov Models`. Research in Computational Molecular Biology (2015): accepted.


### Note and Acknowledgment: 
- UPP bundles the following two programs into its distribution:
  1. hmmer: http://hmmer.janelia.org/
- UPP uses the [Dendropy](http://pythonhosted.org/DendroPy/) package. 
- UPP uses some code from [SATe](http://phylo.bio.ku.edu/software/sate/sate.html).

-------------------------------------
Installation
-------------------------------------
This section details steps for installing and running UPP. We have run UPP on Linux and MAC. If you experience difficulty installing or running the software, please contact one of us (Tandy Warnow, Nam Nguyen, or Siavash Mirarab).

Requirements:
-------------------
Before installing the software you need to make sure the following programs are installed on your machine.

1. Python: Version > 2.6. 

Installation Steps:
-------------------
UPP is a part of the SEPP distribution package.  First install SEPP.  Once done, do the following. 

1. Configure: run `sudo python setup.py upp`. 

The last step creates a ~/.sepp/upp.config config file. Since this is specific to a user, each user that runs sepp needs to execute the last step. 

Common Problems:
-------------------
1.  UPP requires SEPP to be installed.  If UPP is not running, first check to see if UPP was installed correctly.

---------------------------------------------
Running UPP
---------------------------------------------
To run UPP, invoke the `run_upp.py` script from the `bin` sub-directory of the location in which you installed the Python packages. To see options for running the script, use the command:

`python <bin>/run_upp.py -h`

The general command for running TIPP is:

`python <bin>/run_upp.py -s <unaligned_sequences> `

UPP can also be run using a configuration file. Sample configuration files and input files can be found under test/unittest/data/upp/. Change to that directory to run SEPP on the sample files. To run using command options, run

`python <bin>/run_upp.py -s input.fas`

and to run using a configuration file, run

`python <bin>/run_upp.py -c sample.config`

The main outputs of UPP are two alignment files, <prefix>_alignment.fasta and <prefix>_alignment_masked.fasta.  The  <prefix>_alignment.fasta file is the alignment of the query sequences.  The <prefix>_alignment_masked.fasta is the masked
alignment file; sites in the query sequence that are marked as non-homologous to sites in the backbone are removed.  

By setting SEPP_DEBUG environmental variable to `True`, you can instruct SEPP to output more information that can be helpful for debugging.  

---------------------------------------------
Bugs and Errors
---------------------------------------------
UPP is under active research development at UTCS by the Warnow Lab (and especially with her PhD students Siavash Mirarab and Nam Nguyen). Please report any errors to Siavash Mirarab (smirarab@gmail.com) and Nam Nguyen (namphuon@cs.utexas.edu).

