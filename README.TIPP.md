------------------------------------
Summary
------------------------------------

TIPP stands for `Taxonomic identification and phylogenetic profiling`, and so is a method for the following problems:

Taxonomic identification:
- Input: A query sequence `q`
- Output: The taxonomic lineage of `q`.

Abundance profiling:
- Input: A set of query sequences `Q`
- Output: An abundance profile estimated on `Q`


TIPP is a modification of SEPP for classifying query sequences using phylogenetic placement.  TIPP inserts the query sequences into a taxonomic tree and uses the insertion location to identify the reads.  The novel idea behind TIPP is that rather than using the single best alignment and placement for taxonomic identification, we use a collection of alignments and placements.  Our study shows that TIPP provides improved classification accuracy on novel sequences and on sequences with evolutionary divergent datasets.  TIPP can also be used for abundance estimation by computing an abundance profile on the reads binned to the 30 gene reference dataset.

Developers: Tandy Warnow, Nam Nguyen, and Siavash Mirarab

###Publication:
Nguyen, N., Mirarab, S., Liu, B., Pop, M., & Warnow, T.   TIPP: Taxonomic Identification and Phylogenetic Profiling.  Under revision.


### Note and Acknowledgment: 
- TIPP bundles the following two programs into its distribution:
  1. pplacer: http://matsen.fhcrc.org/pplacer/
  2. hmmer: http://hmmer.janelia.org/
  3. EPA: http://sco.h-its.org/exelixis/software.html
- TIPP uses the [Dendropy](http://pythonhosted.org/DendroPy/) package. 
- TIPP uses some code from [SATe](http://phylo.bio.ku.edu/software/sate/sate.html).

-------------------------------------
Installation
-------------------------------------
This section details steps for installing and running TIPP. We have run TIPP on Linux and MAC. If you experience difficulty installing or running the software, please contact one of us (Tandy Warnow, Nam Nguyen, or Siavash Mirarab).

Requirements:
-------------------
Before installing the software you need to make sure the following programs are installed on your machine.

1. Python: Version > 2.6. 
2. Java: Version > 1.5
3. Blast: Version > 2.2.2

Installation Steps:
-------------------
TIPP is a part of the SEPP distribution package.  First install SEPP.  Once done, do the following. 

1. Download the reference dataset available at www.cs.utexas.edu/~phylo/software/sepp/tipp.zip
2. Unzip it to a directory
3. Set the environment variable REFERENCE to point to the location of the reference directory
4. Set the environment variable BLAST to point to the location of blastn
5. Configure: run `sudo python setup.py tipp`. 


The last step creates a ~/.sepp/tipp.config config file. Since this is specific to a user, each user that runs sepp needs to execute the last step. 

Common Problems:
-------------------
1. The last step by default requires root access to the system. If you do not have root access, invoke the setup script as follows: `python setup.py install --prefix=/some/path/on/your/system`, where `/some/path/on/your/system` is the path to a directory on your system to which you do have read and write access. If you use the `--prefix` option, you must ensure that the `lib/python2.x/site-packages` subdirectory (where `x` denotes the minor version number of your Python install) of the directory you specify following `--prefix=` is on Python's search path. To add a directory to Python's search path, modify your PYTHONPATH environment variable.

2. TIPP relies on pplacer and HMMER for alignment and placement steps. These tools are packaged with SEPP. If for some reason the packaged version of HMMER and pplacer do not run in your environment, you need to download and build those programs for your system (see below for links), and point TIPP to them. To point TIPP to your installation of hmmer and placer modify ~/.sepp/tipp.config. 
   pplacer: http://matsen.fhcrc.org/pplacer/
   hmmer://hmmer.janelia.org/

3.  TIPP relies on blastn for the binning of metagenomic reads.  This needs to be installed separately.  To point BLAST to your installation of blastn, modify ~/.sepp/tipp.config. 
   blast: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

3.  TIPP performs abundance profiling uses a set of 30 marker genes.  This needs to be downloaded separately.  Download the reference dataset and unzip it to a directory.  Point the REFERENCE environment variable to this directory before installing TIPP.  You can manually point TIPP to the reference directory by modifying the ~/.sepp/tipp.config file. 
   reference datasets: www.cs.utexas.edu/~phylo/software/sepp/tipp.zip
   
---------------------------------------------
Running TIPP
---------------------------------------------
To run TIPP, invoke the `run_tipp.py` script from the `bin` sub-directory of the location in which you installed the Python packages. To see options for running the script, use the command:

`python <bin>/run_tipp.py -h`

The general command for running SEPP is:

`python <bin>/run_tipp.py -t <tree_file> -a <alignment_file> -f <fragment_file> -r <raxml_info_file> -A <alignment_set_size> -P <placement_set_size>`

TIPP can also be run using a configuration file. Sample configuration files and input files can be found under test/unittest/data/mock/. Change to that directory to run SEPP on the sample files. To run using command options, run

`python <bin>/run_sepp.py -t test.tree -a test.fasta -f test.fas -r test.RAxML_info -A 250 -P 250`

and to run using a configuration file, run

`python <bin>/run_sepp.py -c sample.config`

The main output of SEPP is a .json file, created according to pplacer format. Please refer to pplacer website (currently http://matsen.github.com/pplacer/generated_rst/pplacer.html#json-format-specification) for more information on the format of the josn file. Also note that pplacer package provides a program called guppy that can read .json files and perform downstream steps such as visualization.

In addition to the .json file, SEPP outputs alignments of fragments to reference sets. There could be multiple alignment files created, each corresponding to a different placement subset. 

By setting SEPP_DEBUG environmental variable to `True`, you can instruct SEPP to output more information that can be helpful for debugging.  

---------------------------------------------
Bugs and Errors
---------------------------------------------
SEPP is under active research development at UTCS by the Warnow Lab (and especially with her PhD students Siavash Mirarab and Nam Nguyen). Please report any errors to Siavash Mirarab (smirarab@gmail.com) and Nam Nguyen (namphuon@cs.utexas.edu).

