------------------------------------
Summary
------------------------------------

SEPP stands for `SATe-enabled phylogenetic placement`, and so is a method for the following problem:

- Input: tree `T` and alignment `A` for a set of full-length gene sequences, and set `X` of (potentially fragmentary) query sequences for the same gene

- Output: placement of each fragment in `X` into the tree T, and alignment of each fragment in `X` to the alignment `A`.

SEPP operates by using a divide-and-conquer strategy adopted from SATe-II ([Liu et al., Systematic Biology 2012](http://sysbio.oxfordjournals.org/content/61/1/90.full.pdf+html?sid=dd32838d-89dc-4bda-8008-6f948146341f) and [Liu et. al., Science, 2009](http://www.sciencemag.org/content/324/5934/1561.abstract)) to construct an Ensemble of Hidden Markov Models (HMMs) to represent the input multiple sequence alignment `A`.
It then computes the fit of each query sequence in `X` to each HMM in the ensemble, and uses the highest scoring HMM to add the sequence to the input tree `T`. This technique improves the accuracy of the placements of the query sequences compared to using a single HMM to represent the input alignment. SEPP uses tools in HMMER to construct HMMs, compute the fit of sequences to HMMs, and add sequences to the alignment `A` (code by Sean Eddy). UPP uses pplacer (code by Erick Matsen) to add query sequences to the input tree `T`, after they are added to the alignment `A`.   SEPP is also used in other software, including TIPP (taxonomic identical using phylogenetic placement) and UPP (ultra-large alignments using phylogeny-aware profiles).

Developers: Siavash Mirarab, Tandy Warnow, and Nam Nguyen, with valuable contributions from Uyen Mai, Daniel McDonald and Stefan Janssen.

### Publication:
S. Mirarab, N. Nguyen, and T. Warnow, SEPP: SATe-enabled phylogenetic placement, Proceedings of the Pacific Symposium of Biocomputing 2012, pages 247-58 [http://www.ncbi.nlm.nih.gov/pubmed/22174280#](http://www.ncbi.nlm.nih.gov/pubmed/22174280#).


### Documentations and related pages

- SEPP on green genes: to run SEPP on green genes, it would be easier to use: [wiki](https://github.com/smirarab/sepp/wiki/SEPP-on-Greengenes)
- SEPP [tutorial](tutorial/sepp-tutorial.md).
- 

### Note and Acknowledgment: 
- SEPP bundles the following two programs into its distribution:
  1. [pplacer](http://matsen.fhcrc.org/pplacer/)
  2. [hmmer](http://hmmer.janelia.org/)
  3. [EPA](http://sco.h-its.org/exelixis/software.html)
- SEPP uses the [Dendropy](http://pythonhosted.org/DendroPy/) package. 
- SEPP uses some code from [SATe](http://phylo.bio.ku.edu/software/sate/sate.html).


-------------------------------------
Installation
-------------------------------------
This section details steps for installing and running SEPP. We have run SEPP on Linux and MAC. If you experience difficulty installing or running the software, please contact one of us (Tandy Warnow, Nam Nguyen, or Siavash Mirarab).

Requirements:
-------------------
Before installing the software you need to make sure the following programs are installed on your machine.

1. Python: Version 2.7 or later (including python 3). 
2. Java: Version > 1.5

Installation Steps:
-------------------
SEPP is distributed as Python source code. Once you have the above required software installed, do the following. 

**Note:** these installation steps have recently changed

1. Obtain the latest SEPP distribution from git repository (using `git clone` or by simply downloading the Zip file). If you downloaded the zip file, uncompress the distribution file.
2. Go to the distribution directory
3. Configure: run `python setup.py config` (or `python setup.py config -c` to avoid using the home directory). 
4. Install: run `python setup.py install`. 

The third step creates a `~/.sepp/` directory, puts the default config file under `~/.sepp/main.config`, and puts all the binary executables under it as well. 

* Since this (default location) is specific to a user, each user that runs sepp needs to execute the last step. 
* To have a self-contained install instead of using the home directory, you should run step 3 with the `-c` option.
With `-c`, the main config file and the binary files are kept under the place where you have downloaded SEPP. Therefore that
directory should not be removed later if `-c` is used. 
* To use any other location for the config file and the binary files, change the file `home.path` after step 3 and before step 4. 
You need to also update the paths inside the `.sepp/main.config` file after step 3 and before step 4. 

Common Problems:
-------------------
1. The last step by default requires root access to the system. If you do not have root access, invoke the setup script as follows: `python setup.py install --prefix=/some/path/on/your/system`, where `/some/path/on/your/system` is the path to a directory on your system to which you do have read and write access. If you use the `--prefix` option, you must ensure that the `lib/python2.x/site-packages` subdirectory (where `x` denotes the minor version number of your Python install) of the directory you specify following `--prefix=` is on Python's search path. To add a directory to Python's search path, modify your PYTHONPATH environment variable.

2. SEPP relies on pplacer and HMMER for alignment and placement steps. These tools are packaged with SEPP. If for some reason the packaged version of HMMER and pplacer do not run in your environment, you need to download and build those programs for your system (see below for links), and point SEPP to them. To point sepp to your installation of hmmer and placer modify `~/.sepp/main.config`. 
   * pplacer: <http://matsen.fhcrc.org/pplacer/>
   * hmmer: <http://hmmer.janelia.org/>


---------------------------------------------
Running SEPP
---------------------------------------------
To run SEPP, invoke the `run_sepp.py` script from the `bin` sub-directory of the location in which you installed the Python packages. To see options for running the script, use the command:

`python <bin>/run_sepp.py -h`

The general command for running SEPP is:

`python <bin>/run_sepp.py -t <tree_file> -a <alignment_file> -f <fragment_file> -r <raxml_info_file> -A <alignment_set_size> -P <placement_set_size> `

SEPP can also be run using a configuration file. Sample configuration files and input files can be found under `test/unittest/data/simulated/`. Change to that directory to run SEPP on the sample files. To run using command options, run

`python <bin>/run_sepp.py -t test.tree -a test.fasta -f test.fas -r test.RAxML_info -A 250 -P 250`

and to run using a configuration file, run

`python <bin>/run_sepp.py -c sample.config`

**Output:** 

The main output of SEPP is a .json file, created according to pplacer format. Please refer to pplacer website (currently http://matsen.github.com/pplacer/generated_rst/pplacer.html#json-format-specification) for more information on the format of the josn file. Also note that pplacer package provides a program called guppy that can read .json files and perform downstream steps such as visualization.

In addition to the .json file, SEPP outputs alignments of fragments to reference sets. There could be multiple alignment files created, each corresponding to a different placement subset. 

Finally, SEPP internally renames internal node names to safe names because guppy cannot handle some complicated names. SEPP outputs a python script that can 
be used to rename your internal nodes back to the original value. To update the json file using this script, you can run:

```
cat [the name of .json/.tre/.xml file with mapped names]| python output__rename-json.py > [name of the relabelled file]
```

By setting SEPP_DEBUG environmental variable to `True`, you can instruct SEPP to output more information that can be helpful for debugging.  

---------------------------------------------
Bugs and Errors
---------------------------------------------
SEPP is under active research development at UIUC by the Warnow Lab (and especially with her former PhD students Siavash Mirarab and Nam Nguyen). Please report any errors to Siavash Mirarab (smirarab@gmail.com) and Nam Nguyen (ndn006@eng.ucsd.edu).


