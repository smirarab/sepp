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
- UPP bundles the following program into its distribution:
  1. hmmer: http://hmmer.janelia.org/    
- UPP uses the [Dendropy](http://pythonhosted.org/DendroPy/) package. 
- UPP uses [PASTA] (https://github.com/smirarab/PASTA/) to estimate the backbone alignment and tree 
- UPP uses some code from [SATe](http://phylo.bio.ku.edu/software/sate/sate.html).

-------------------------------------
Installation
-------------------------------------
This section details steps for installing and running UPP. We have run UPP on Linux and MAC. If you experience difficulty installing or running the software, please contact one of us (Tandy Warnow, Nam Nguyen, or Siavash Mirarab).

Requirements:
-------------------
Before installing the software you need to make sure the following programs are installed on your machine.  

1. Python: Version > 2.6. 
2. SEPP: Version > 1.0. 
3. PASTA: Version > 1.0. 

Installation Steps:
-------------------
UPP is a part of the SEPP distribution package.  First install SEPP (see [SEPP readme] (https://github.com/smirarab/sepp/blob/master/README.SEPP.md)).  Next install [PASTA] (https://github.com/smirarab/PASTA/).  Once done, do the following. 

1. Configure: run `sudo python setup.py upp`. 

The last step creates a ~/.sepp/upp.config config file. Since this is specific to a user, each user that runs sepp needs to execute the last step. 

Common Problems:
-------------------
1.  UPP requires SEPP to be installed.  If UPP is not running, first check to see if SEPP was installed correctly.
2.  UPP requires PASTA to be installed and the run_pasta.py executable to be on the path.  

---------------------------------------------
Running UPP
---------------------------------------------
To run UPP, invoke the `run_upp.py` script from the `bin` sub-directory of the location in which you installed the Python packages. To see options for running the script, use the command:

`python <bin>/run_upp.py -h`

The general command for running TIPP is:

`python <bin>/run_upp.py -s <unaligned_sequences>`

By default, this will automatically select 1,000 sequences to be in the backbone set, generate a PASTA alignment and tree, and then align the remaining sequences to the backbone alignment.  UPP can also be run using a configuration file. Sample configuration files and input files can be found under test/unittest/data/upp/. Change to that directory to run SEPP on the sample files. To run using command options, run

The main outputs of UPP are two alignment files, <prefix>_alignment.fasta and <prefix>_alignment_masked.fasta.  The  <prefix>_alignment.fasta file is the alignment of the unaligned sequences.  The <prefix>_alignment_masked.fasta is the masked alignment file; non-homologous sites in the query set are removed.  

The secondary outputs are the backbone alignment and tree (always named as pasta.fasta and pasta.fasttree).

To run a small test example with 1,000 sequences, run the following command:

`python <bin>/run_upp.py -s initial.fas -B 100`

This will generate a backbone set of 100 sequences and align the remaining 900 sequences to the backbone alignment.

To run using a configuration file, run

`python <bin>/run_upp.py -c sample.config`

By setting SEPP_DEBUG environmental variable to `True`, you can instruct SEPP to output more information that can be helpful for debugging.  

---------------------------------------------
Advance options
---------------------------------------------
To run UPP with a pre-computed backbone alignment and tree, run

`python <bin>/run_upp.py -s input.fas -a <alignment_file> -t <tree_file>`

To filter fragments from the backbone selection process, run

`python <bin>/run_upp.py -s input.fas -M <median_full_length>`

UPP will only include sequences within 25% of the median length in the backbone set.

---------------------------------------------
Bugs and Errors
---------------------------------------------
UPP is under active research development at UIUC by the Warnow Lab (and especially with her PhD students Siavash Mirarab and postdoc Nam Nguyen). Please report any errors or requests to Siavash Mirarab (smirarab@gmail.com) and Nam Nguyen (namphuon@illinois.edu).

