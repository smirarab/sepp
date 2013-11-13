------------------------------------
Summary
------------------------------------

UPP stands for `Ultra-large alignment using SEPP`, and so is a method for the following problem:

- Input: tree `T` and alignment `A` for a set full-length sequences, and set `Q` of full-length sequences to align.

- Output: alignment of each query sequence `Q` to the alignment `A`.

UPP operates by using a divide-and-conquer strategy adopted from SATe [Liu et. al., Science, 2009](http://www.sciencemag.org/content/324/5934/1561.abstract) to improve the alignment produced by running HMMER (code by Sean Eddy). Our study shows that UPP provides improved alignment accuracy compared to other profile based methods and can efficiently align very large datasets.

Developers: Tandy Warnow, Nam Nguyen, and Siavash Mirarab

###Publication:
N. Nguyen, S. Mirarab, and T. Warnow, UPP: Ultra-large alignment using SEPP, in preparation.

### Note and Acknowledgment: 
- UPP bundles the following two programs into its distribution:
  1. hmmer: http://hmmer.janelia.org/
- UPP uses the [Dendropy](http://pythonhosted.org/DendroPy/) package. 
- UPP uses some code from [SATe](http://phylo.bio.ku.edu/software/sate/sate.html).
- UPP requires SATe for generating alignments and backbones if none are provided from the user.

-------------------------------------
Installation
-------------------------------------
This section details steps for installing and running UPP. We have run UPP on Linux and MAC. If you experience difficulty installing or running the software, please contact one of us (Tandy Warnow, Nam Nguyen, or Siavash Mirarab).

Note that by installing UPP, you also install SEPP.

Requirements:
-------------------
Before installing the software you need to make sure the following programs are installed on your machine.

1. Python: Version > 2.6. 
2. Java: Version > 1.5

Optional (Automatic backbone alignment and tree generation)
-------------------
If you want to have UPP automatically estimate backbone alignments and trees for alignment, you need to have SATe version 2.2.7 or greater installed on your machine before the installation of UPP.

Set the SATE environment variable to point to the SATe executable, i.e., the run_sate.py python script.  This will cause the default configuration file to include the SATE location and will allow the automatic use of SATE for backbone alignment and tree generation.

If you want later have the capability of automatically generating backbone alignments and trees, then the default configuration file (found in ~/.sepp/main.config by default) should be changed to include the following line:
[sate]
path=/path/to/run_sate.py


Installation Steps:
-------------------
UPP is distributed as Python source code. Once you have the above required software installed, do the following. 

1. Obtain the latest UPP distribution from git repository (using `git clone` or by simply downloading the Zip file). If you downloaded the zip file, uncompress the distribution file.
2. Go to the distribution directory
3. Install: run `sudo python setup.py install`. 
4. Configure: run `sudo python setup.py config`. 

The last step creates a ~/.sepp/ directory and put the default config file under ~/.sepp/main.config. Since this is specific to a user, each user that runs sepp needs to execute the last step. 

Common Problems:
-------------------
1. The last step by default requires root access to the system. If you do not have root access, invoke the setup script as follows: `python setup.py install --prefix=/some/path/on/your/system`, where `/some/path/on/your/system` is the path to a directory on your system to which you do have read and write access. If you use the `--prefix` option, you must ensure that the `lib/python2.x/site-packages` subdirectory (where `x` denotes the minor version number of your Python install) of the directory you specify following `--prefix=` is on Python's search path. To add a directory to Python's search path, modify your PYTHONPATH environment variable.

2. UPP relies on HMMER for alignment. This tool is packaged with UPP. If for some reason the packaged version of HMMER and pplacer do not run in your environment, you need to download and build those programs for your system (see below for links), and point UPP to them. To point sepp to your installation of hmmer modify ~/.sepp/main.config. 
   hmmer://hmmer.janelia.org/

3. UPP requires SATe installed so that backbone alignments and trees can be generated when none are provided.  SATe is available from [here](http://phylo.bio.ku.edu/software/sate/sate.html).

---------------------------------------------
Running UPP
---------------------------------------------
To run UPP, invoke the exhaustive_upp.py` script from the `bin` sub-directory of the location in which you installed the Python packages. To see options for running the script, use the command:

`python <bin>/exhausive_upp.py -h`

The general command for running UPP is:

`python <bin>/exhaustive_upp.py -t <tree_file> -a <alignment_file> -s <query_sequence_file>`

where the tree_file is the backbone tree, the alignment_file is the backbone alignment, and the query_sequence_file is the remaining sequences to be aligned.

UPP can also be run using a configuration file. Sample configuration files and input files can be found under test/unittest/data/upp/. Change to that directory to run UPP on the sample files. To run using command options, run

`python <bin>/exhaustive_upp.py -t test.tree -a test.fasta -s query.fas`

and to run using a configuration file, run

`python <bin>/exhaustive_upp.py -c sample.config`

The main output of UPP are two alignment .fasta files, *_alignment.fasta and *_alignment_masked.fasta.  The first file is the alignment of the query sequences to the backbone.  The second file is the same alignment with insertion columns masked out.  We suggest using the second file for phylogenetic inference as characters found in insertion columns should not be treated as homologous. 

Finally, UPP can automatically estimate the backbone alignment and tree.  To do so, run UPP with just the `-s` command; UPP will select a backbone alignment size and generate a backbone alignment and tree on a random sampling of the sequences.  For example, 

`python <bin>/exhaustive_upp.py -s initial.fas`

Note that this command may take a while because a SATe tree will be estimated on a backbone size of 100 sequences.  To speed this up, use a backbone of size 20 (`-B 20`).

By setting SEPP_DEBUG environmental variable to `True`, you can instruct UPP to output more information that can be helpful for debugging.  

---------------------------------------------
Bugs and Errors
---------------------------------------------
UPP is under active research development at UTCS by the Warnow Lab (and especially with her PhD students Siavash Mirarab and Nam Nguyen). Please report any errors to Siavash Mirarab (smirarab@gmail.com) and Nam Nguyen (namphuon@cs.utexas.edu).

