Introduction
===

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

---
Installation
===

You have two options for installing TIPP. 

 - **Windows:** If you have a Windows machine, currently using the Virtual Machine (VM) image we provide is your only option. 
 - **Linux:** and **MAC:**  If you have Linux (or other \*nix systems) or MAC, you can still use VM, but downloading the code from github and installing it is what we strongly recommend. 
 
 ### 1. From Source Code
Current version of TIPP has been developed and tested entirely on Linux and MAC. 
Windows won't work currently (future versions may or may not support Windows). 

You need to have:

- Python (version 2.7 or later)
- Blast (version 2.2.2 or later)
- Java (version > 1.5 or later)

**Installation of SEPP**:
TIPP is a part of the SEPP distribution package.  First download and install SEPP:

1. Open a terminal and create a directory where you want to keep SEPP. e.g. `mkdir ~/sepp-code`. Go to this directory. e.g. `cd ~/sepp-code`.

2. Clone the SEPP code repository from our [github repository](https://github.com/smirarab/sepp). For example you can use `git clone https://github.com/smirarab/sepp.git`.
If you don't have git, you can directly download a [zip file from the repository](https://github.com/smirarab/sepp/archive/master.zip)
and decompress it into your desired directory. 

3. `cd sepp` (or `cd sepp-master` if you used the zip file instead of cloning the git repository)

4. Then run:

```
 sudo python setup.py install
```
 
If you don't have root access, remove the `sudo` part and instead  use  `--user` option. Alternativley, you can `--prefix` to install in a different location, but that different location needs to be part of your `PYTHONPATH` environmental variable. 

5. Run the following command:


```
python setup.py config
```

**Installation of TIPP**:

Once done, do the following. 

1. Download the reference dataset available at www.cs.utexas.edu/~phylo/software/sepp/tipp.zip
2. Unzip it to a directory
3. Set the environment variable REFERENCE to point to the location of the reference directory.  This can be performed using:

```
 REFERENCE=/PATH/TO/REFERENCE
 export REFERENCE
```

4. Set the environment variable BLAST to point to the directory containing the location of blastn.  This can be performed using:

```
 BLAST=/PATH/TO/BLASTN/DIRECTORY
 export BLAST
```

5. Run the following command:

```
python setup.py tipp
```

The last step creates a ~/.sepp/tipp.config config file. Since this is specific to a user, each user that runs tipp needs to execute the last step. 

**Common Problems:**
1.  TIPP requires SEPP to be installed.  If TIPP is not running, first check to see if TIPP was installed correctly.

2.  TIPP relies on blastn for the binning of metagenomic reads.  This needs to be installed separately.  To point BLAST to your installation of blastn, modify ~/.sepp/tipp.config. 
   blast: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

3.  TIPP performs abundance profiling uses a set of 30 marker genes.  This needs to be downloaded separately.  Download the reference dataset and unzip it to a directory.  Point the REFERENCE environment variable to this directory before installing TIPP.  You can manually point TIPP to the reference directory by modifying the ~/.sepp/tipp.config file. 
   reference datasets: www.cs.utexas.edu/~phylo/software/sepp/tipp.zip


### 2. From Virtual Machine (VM)

VM Image (mostly for Windows users) is available for [download](http://www.cs.utexas.edu/~phylo/software/PASTA_TIPP_UPP.ova) (2.5 GB). Once the image is downloaded, you need to run it using a VM environment ([VirtualBox](https://www.virtualbox.org/) is a good option). After you install VirtualBox, you just need to use File/import to import the PASTA_TIPP_UPP.ova image that you have downloaded (If your machine has less than 3GB you might want to reduce the memory to something like 512MB). Once VM is imported, you can start it from the Virtualbox. If you are asked to login, the username and passwords are (username: phylolab, password: phylolab). TIPP and UPP are already installed on the VM machine, so you can simply proceed by opening a terminal and running it.

Email `sepp-users@googlegroups.com` for installation issues. 


---------
Contact
===
Email: `ensemble-of-hmms@googlegroups.com` for all issues. 
