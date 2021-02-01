Taxonomic Identification and Phylogenetic Profiling (TIPP)
==========================================================
TIPP is a method for the following problems:

Taxonomic identification:
+ Input: A query sequence *q*
+ Output: The taxonomic lineage of *q*

Abundance profiling:
+ Input: A set *Q* of query sequences
+ Output: An abundance profile estimated on *Q*

TIPP is a modification of SEPP for classifying query sequences (i.e. reads) using phylogenetic placement. TIPP inserts each read into a taxonomic tree and uses the insertion location to identify the taxonomic lineage of the read. The novel idea behind TIPP is that rather than using the single best alignment and placement for taxonomic identification, we use a collection of alignments and placements and consider statistical support for each alignment and placement. Our study shows that TIPP provides improved classification accuracy on novel sequences and on sequences with evolutionarily divergent datasets. TIPP can also be used for abundance estimation by computing an abundance profile on the reads binned to marker genes in a reference dataset. TIPP2 provides an new reference dataset with 40 marker genes, assembled from the NCBI RefSeq database (learn more [here](https://github.com/shahnidhi/TIPP_reference_package)). In addition, TIPP2 updates how query sequences (i.e. reads) are mapped to marker genes. This repository corresponds to TIPP2, and henceforth we use the terms TIPP and TIPP2 interchangeably.

Developers of TIPP: Nam Nguyen, Siavash Mirarab, Nidhi Shah, Erin Molloy, and Tandy Warnow.

### Publications:
Nguyen, Nam, Siavash Mirarab, Bo Liu, Mihai Pop, and Tandy Warnow, "TIPP: Taxonomic identification and phylogenetic profiling," *Bioinformatics*, 2014. [doi:10.1093/bioinformatics/btu721](http://bioinformatics.oxfordjournals.org/content/30/24/3548.full.pdf).

Shah, Nidhi, Erin K. Molloy, Mihai Pop, and Tandy Warnow, "TIPP2: metagenomic taxonomic profiling using phylogenetic markers," *Bioinformatics*, 2020. [doi:10.1093/bioinformatics/btab023](https://doi.org/10.1093/bioinformatics/btab023)

### Note and Acknowledgment: 
- TIPP bundles the following two programs into its distribution:
	- pplacer: http://matsen.fhcrc.org/pplacer/
	- hmmer: http://hmmer.janelia.org/
	- EPA: http://sco.h-its.org/exelixis/software.html
- TIPP uses the [Dendropy](http://pythonhosted.org/DendroPy/) package. 
- TIPP uses some code from [SATe](http://phylo.bio.ku.edu/software/sate/sate.html).

-------------------------------------

Installing TIPP
===============
This section details steps for installing and running TIPP. We have run TIPP on Linux and MAC. If you experience difficulty installing or running the software, please contact one of us (Tandy Warnow or Siavash Mirarab).

Requirements:
-------------
Before installing the software you need to make sure the following programs are installed on your machine.

- Python: Version > 2.7 
- Java: Version > 1.5
- Blast: Version > 2.2.2

Installation Steps:
-------------------
TIPP is a part of the SEPP distribution package. First install SEPP. Once done, do the following. 

1. Download and decompress the reference dataset available at [https://obj.umiacs.umd.edu/tipp/tipp2-refpkg.tar.gz](https://obj.umiacs.umd.edu/tipp/tipp2-refpkg.tar.gz).
2. Set the environment variables (`REFERENCE` and `BLAST`) that will be used to create the configuration file. The `REFERENCE` environment variable should point to the location of the reference dataset (i.e. it should point to the `tipp2-refpkg` directory with its full path). The `BLAST` environment variable should to point to the location of the binary: `blastn`. Environment variables can be set using the following (shell-dependent) commands:
	- `export VARIABLE_NAME=/path/to/file` (bash shell)
	- `setenv VARIABLE_NAME /path/to/file` (tcsh shell)
3. Create the TIPP configuration file by running the command: `python setup.py tipp` or `python setup.py tipp -c`. This  creates a `tipp.config` config file. It is important that you use `-c` here if you used `-c` when installing SEPP and otherwise, not use `-c`. 


Common Problems:
----------------
1. TIPP requires SEPP to be installed. If TIPP is not running, first check to see if SEPP was installed correctly by typing `run_sepp.py -h` to see it's user options.
2. If you are installing sepp with the `--user` flag, then you may need to set your `PYTHONPATH` environmental variable. To run TIPP, you may also need to include your python bin directory in your `PATH` environmental variable. If you are having issues installing or running sepp, try adding the following lines to your bash profile file `~/.bash_profile`:
```
export PYTHONPATH="$HOME/.local/lib/python3.7/site-packages:$PYTHONPATH"
export PATH="$HOME/.local/bin:$PATH"
```
where `3.7` is replaced with the version of python that you are running.

3. Several cluster systems that we have worked on have different commands to run different versions of python e.g. `python` runs python version 2.X and `python3` runs python version 3.Y. In this scenario, you need to be careful to be consistent with the python command (e.g. `python` or `python3`) when installing/configuring sepp and tipp.

4. TIPP relies on `blastn` for the binning of metagenomic reads, so BLAST needs to be downloaded and installed separately (learn more [here](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)). Then, point the `BLAST` environment variable to your installation of `blastn`. Alternatively, you can manually point TIPP to the `blastn` installation by modifying the `tipp.config` file. 

5. TIPP performs abundance profiling uses a set of 40 marker genes. This reference dataset needs to be downloaded separately from [here](https://obj.umiacs.umd.edu/tipp/tipp2-refpkg.tar.gz). Then, point the `REFERENCE` environment variable to the decompressed directory before installing TIPP. Alternatively, you can manually point TIPP to the reference dataset by modifying the `tipp.config` file. 

---------------------------------------------

Running TIPP
============
TIPP is a general pipeline for classifying reads belonging to a specific marker gene.  We provide precomputed marker gene datasets for a collection of genes found in the tipp.zip archive.  

The general command for running TIPP for a specific marker is:

`run_tipp.py -R <reference_marker> -f <fragment_file>`

To see options for running the script, use the command:

`run_tipp.py -h`

The main output of TIPP is a `_classification.txt` file that annotation for each read. 
TIPP also outputs other information about the fragments, including the alignments (note that there could be multiple alignment files created, each corresponding to a different placement subset) as well as the phylogenetic placements.
Placements are given in the `.json` file, created according to *pplacer* format. Please refer to [pplacer website](http://matsen.github.com/pplacer/generated_rst/pplacer.html#json-format-specification) for more information on the form at of the `.json` file. Also note that *pplacer* package provides a program called *guppy* that can read `.json` files and  perform downstream steps such as visualization.

By setting `SEPP_DEBUG` environment variable to `True`, you can instruct SEPP to output more information that can be helpful for debugging.  

This [tutorial](tutorial/tipp-tutorial.md) contains examples of running TIPP for read classification as well as abundance profiling using the script `run_abundance.py`.

---------------------------------------------

Bugs and Errors
===============
TIPP is under active research development at UIUC by the Warnow Lab. Please report any errors to Tandy Warnow (warnow@illinois.edu) and Siavash Mirarab (smirarab@ucsd.edu).

