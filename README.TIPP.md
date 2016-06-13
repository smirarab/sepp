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


TIPP is a modification of SEPP for classifying query sequences using phylogenetic placement.  TIPP inserts the query sequences into a taxonomic tree and uses the insertion location to identify the reads.  The novel idea behind TIPP is that rather than using the single best alignment and placement for taxonomic identification, we use a collection of alignments and placements and consider statistical support for each alignment and placement.  Our study shows that TIPP provides improved classification accuracy on novel sequences and on sequences with evolutionarily divergent datasets.  TIPP can also be used for abundance estimation by computing an abundance profile on the reads binned to the 30 gene reference dataset.

Developers: Nam Nguyen, Siavash Mirarab, and Tandy Warnow.

###Publication:
Nguyen, Nam , Siavash Mirarab, Bo Liu, Mihai Pop, and Tandy Warnow. `TIPP: Taxonomic identification and phylogenetic profiling`. Bioinformatics (2014). [doi:10.1093/bioinformatics/btu721](http://bioinformatics.oxfordjournals.org/content/30/24/3548.full.pdf).


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

1. Python: Version > 2.7. 
2. Java: Version > 1.5
3. Blast: Version > 2.2.2

Installation Steps:
-------------------
TIPP is a part of the SEPP distribution package.  First install SEPP.  Once done, do the following. 

1. Download the reference dataset available at www.cs.utexas.edu/~phylo/software/sepp/tipp.zip
2. Unzip it to a directory
3. Set the environmental variables that will be used to create the configuration file.  The environmental variable can be set using the following command (shell-dependent)
  `export VARIABLE_NAME=/path/to/file`  (bash shell)
  `setenv VARIABLE_NAME /path/to/file` (tcsh shell)  
  2a. Set the environment variable REFERENCE to point to the location of the reference directory (i.e., the location of the tipp folder generated from unzipping the tipp.zip file)
  2b. Set the environment variable BLAST to point to the location of blastn
4. Configure: run `sudo -E python setup.py tipp`. 
  


The last step creates a ~/.sepp/tipp.config config file. Since this is specific to a user, each user that runs TIPP needs to execute the last step. 

Common Problems:
-------------------
1.  TIPP requires SEPP to be installed.  If TIPP is not running, first check to see if TIPP was installed correctly.

2.  TIPP relies on blastn for the binning of metagenomic reads.  This needs to be installed separately.  To point BLAST to your installation of blastn, modify ~/.sepp/tipp.config. 
   blast: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

3.  TIPP performs abundance profiling uses a set of 30 marker genes.  This needs to be downloaded separately.  Download the reference dataset and unzip it to a directory.  Point the REFERENCE environment variable to this directory before installing TIPP.  You can manually point TIPP to the reference directory by modifying the ~/.sepp/tipp.config file. 
   
---------------------------------------------
Running TIPP
---------------------------------------------
To run TIPP, invoke the `run_tipp.py` script from the `bin` sub-directory of the location in which you installed the Python packages. To see options for running the script, use the command:

`python <bin>/run_tipp.py -h`

TIPP is a general pipeline for classifying reads belonging to a specific marker gene.  We provide precomputed marker gene datasets for a collection of genes found in the tipp.zip archive.  

The general command for running TIPP for a specific marker is:

`python <bin>/run_tipp.py -R <reference_marker> -f <fragment_file>`

The main output of TIPP is a _classification.txt file that annotation for each read.  In addition, TIPP outputs a .json file with the placements, created according to pplacer format. Please refer to pplacer website (currently http://matsen.github.com/pplacer/generated_rst/pplacer.html#json-format-specification) for more information on the format of the josn file. Also note that pplacer package provides a program called guppy that can read .json files and perform downstream steps such as visualization.

In addition to the .json file, TIPP outputs alignments of fragments to reference sets. There could be multiple alignment files created, each corresponding to a different placement subset. 

By setting SEPP_DEBUG environmental variable to `True`, you can instruct SEPP to output more information that can be helpful for debugging.  

---------------------------------------------
Bugs and Errors
---------------------------------------------
TIPP is under active research development at UIUC by the Warnow Lab (and especially with her PhD students Siavash Mirarab and Nam Nguyen). Please report any errors to Siavash Mirarab (smirarab@gmail.com) and Nam Nguyen (namphuon@illinois.edu).

