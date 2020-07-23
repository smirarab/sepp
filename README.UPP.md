------------------------------------
Summary
------------------------------------

UPP stands for `Ultra-large alignments using Phylogeny-aware Profiles`, and so is a method for the following problem:

Alignment:
- Input: A set of unaligned sequences `S`
- Output: An alignment `A` on `S`

UPP is a modification of SEPP for performing alignments of ultra-large and fragmentary datasets.  UPP operates in four steps.  In the first step, UPP partitions set `S` into a backbone set and a query set and computes an alignment and tree on the backbone set using [PASTA](https://github.com/smirarab/pasta) (Mirarab et al., RECOMB 2014 and Journal of Computational Biology 2014), which is a direct improvement to SATe (Liu et al., Science 2009 and Systematic Biology 2012).  In the next step, UPP decomposes the backbone alignment into an ensemble of profile Hidden Markov Models (HMMs).  The third step in UPP searches for the best alignment of the query sequence to each HMM.  The final step inserts the query sequence into the backbone alignment using the best scoring HMM.  Our study shows that UPP results in accurate alignments, and that ML trees estimated on the alignments are also highly accurate. UPP has good accuracy on datasets that contain fragmentary sequences. 

UPP(Default): The default version selects 1000 sequences at random for the backbone alignment. If the dataset has at most 1000 sequences, this means that UPP(Default) is identical to PASTA. 

UPP(Fast): We have designed a fast version of UPP that uses a backbone with at most 100 sequences. The default version uses a backbone of 1000 sequences. The fast version can produce an alignment on 10,000 sequences in less than an hour using 12 processors, and on 1,000,000 sequences in less than 12 days; the default version requires seven hours on 10,000 sequences and would take an estimated 120 days on 1,000,000 sequences. The difference in accuracy between UPP(Fast) and UPP(Default) depends on the rate of evolution -- datasets with low to moderate evolutionary diameters can be analyzed well with UPP(Fast); otherwise, we recommend the use of UPP(Default). However, on large datasets, UPP(Default) will take nearly ten times as much running time.

FRAGMENTARY DATASETS: UPP can be used in default mode, which will select the backbone sequences randomly and without trying to restrict the backbone to full length sequences. However, if the dataset contains fragments, then UPP should be used in a mode that restricts the backbone to just the "full-length" sequences. To do this, you will need to provide UPP with an estimate of the full length of sequences for your locus. See Advanced Usage information about how to do this.

SUPPLYING YOUR OWN SEED ALIGNMENT AND TREE: If you have a pre-computed seed alignment (with phylogenetic tree), you can provide this to UPP. See Advanced Usage information below about how to do this.

PARALLEL IMPLEMENTATION: UPP is embarrassingly parallel. See Advanced Usage information about how to do this.



Developers: Nam Nguyen, Siavash Mirarab, and Tandy Warnow with valuable contributions from Metin Balaban.

###Publication:
Nam Nguyen, Siavash Mirarab, Keerthana Kumar, and Tandy Warnow. `Ultra-large alignments using Phylogeny Aware Profiles`. Accepted to RECOMB 2015 (Research in Computational Molecular Biology 2015) and Genome Biology.


### Note and Acknowledgment: 
- UPP bundles the following program into its distribution:
  1. hmmer: http://hmmer.janelia.org/    
- UPP uses the [Dendropy](http://pythonhosted.org/DendroPy/) package. 
- UPP uses [PASTA] (https://github.com/smirarab/PASTA/) to estimate the backbone alignment and tree 
- UPP uses some code from [SATe](http://phylo.bio.ku.edu/software/sate/sate.html).

-------------------------------------
Installation
-------------------------------------
This section details steps for installing and running UPP. We have run UPP on Linux and MAC. If you experience difficulty installing or running the software, please contact Nam Nguyen or Siavash Mirarab.

Requirements:
-------------------
Before installing the software you need to make sure the following programs are installed on your machine.  

1. Python: Version > 2.7 (including python 3). 
2. SEPP: Version > 3.0. 
3. PASTA: Version > 1.0. 

Installation Steps:
-------------------
UPP is a part of the SEPP distribution package.  

1. Install SEPP (see [SEPP readme](https://github.com/smirarab/sepp/blob/master/README.SEPP.md)).  
2. Install [PASTA](https://github.com/smirarab/PASTA/) and make sure the run_pasta.py executable is on the PATH variable.
4. Configure: run `python setup.py upp` or `python setup.py upp -c` (you should use `-c` if you used `-c` when you installed SEPP). 

The last step creates an `upp.config` config file. It is important that you use `-c` here if you used `-c` when installing SEPP and otherwise, not use `-c`. 

Common Problems:
-------------------
1.  UPP requires SEPP to be installed.  If UPP is not running, first check to see if SEPP was installed correctly.
2.  UPP requires PASTA to be installed and the run_pasta.py executable to be on the path.  

---------------------------------------------
Running UPP
---------------------------------------------
To run UPP, invoke the `run_upp.py` script from the `bin` sub-directory of the location in which you installed the Python packages. To see options for running the script, use the command:

`python <bin>/run_upp.py -h`

The general command for running UPP is:

`python <bin>/run_upp.py -s <unaligned_sequences>`

This will run UPP(Default) as described in the main paper.  This will automatically select up to 1,000 sequences to be in the backbone set, generate a PASTA alignment and tree, and then align the remaining sequences to the backbone alignment.  UPP can also be run using a configuration file. 

The main outputs of UPP are two alignment files, `<prefix>_alignment.fasta` and `<prefix>_alignment_masked.fasta`.  The  `<prefix>_alignment.fasta` file is the alignment of the unaligned sequences.  The `<prefix>_alignment_masked.fasta` is the masked alignment file; non-homologous sites in the query set are removed.  

The secondary outputs are the backbone alignment and tree (always named as pasta.fasta and pasta.fasttree) and the list of insertion columns (named `<prefix>_insertion_columns.txt`).

Sample configuration files and input files can be found under `test/unittest/data/upp/`. Change to that directory to run UPP on the sample files.  To run UPP(Fast) on a small test example with 1,000 sequences, run the following command from the `test/unittest/data/upp/` directory:

`python <bin>/run_upp.py -s initial.fas -B 100`

This will generate a backbone set of 100 sequences and align the remaining 900 sequences to the backbone alignment.

To run using a configuration file, run

`python <bin>/run_upp.py -c sample.config`

By setting SEPP_DEBUG environmental variable to `True`, you can instruct SEPP to output more information that can be helpful for debugging.  

---------------------------------------------
Advanced Usage Options
---------------------------------------------
To run UPP(Fast) as described in the main paper, run

`python <bin>/run_upp.py -s input.fas -B 100`

UPP currently assumes that the input sequences are nucleotide sequences.  To
select the input datatype, run

`python <bin>/run_upp.py -s input.fas -m [dna|rna|amino]`

To run UPP with a pre-computed backbone alignment and tree, run

`python <bin>/run_upp.py -s input.fas -a <alignment_file> -t <tree_file>`

If you have a pre-computed alignment, but not a tree, FastTree can be
run to generate the tree (FastTree is packaged with PASTA).  The general command
for running FastTree is:

FastTree backbone_alignment.fasta > backbone_tree

To filter fragments from the backbone selection process, run

`python <bin>/run_upp.py -s input.fas -M <median_full_length>`

UPP will only include sequences in the backbone set that are within 25% of the median length provided.

To run the parallelized version of UPP, run

`python <bin>/run_upp.py -s input.fas -x <cpus>`

If nucleotide sequences are known for input amino acid sequences (backbone and query), backtranslation of extended 
alignment is performed via the following command:

`python <bin>/run_upp.py -s input.fas -a <alignment_file> -t <tree_file> -b nucleotide_sequences.fas`

where nucleotide_sequences.fas contains unaligned DNA sequences of both backbone and query sequences. 
The command will create two DNA alignment files: `<prefix>_backtranslated_alignment.fasta` and 
`<prefix>_backtranslated_alignment_masked.fasta` .

---------------------------------------------
Bugs and Errors
---------------------------------------------
UPP is under active research development at UIUC by the Warnow Lab (and especially with her former PhD students Siavash Mirarab and Nam Nguyen). Please report any errors or requests to Siavash Mirarab (smirarab@gmail.com) and Nam Nguyen (ndn006@eng.ucsd.edu).

