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

5. Run the following command from the SEPP directory:


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

5. Run the following command from the SEPP directory:

```
python setup.py tipp
```

The last step creates a ~/.sepp/tipp.config config file. Since this is specific to a user, each user that runs tipp needs to execute the last step. 

**Common Problems:**
1.  TIPP requires SEPP to be installed.  If TIPP is not running, first check to see if TIPP was installed correctly.

2.  TIPP relies on blastn for the binning of metagenomic reads.  This needs to be installed separately.  To point BLAST to your installation of blastn, modify ~/.sepp/tipp.config. 
   blast: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
   You can also download the linux version of blastn from: http://www.cs.utexas.edu/~phylo/software/blastn

3.  TIPP performs abundance profiling uses a set of 30 marker genes.  This needs to be downloaded separately.  Download the reference dataset and unzip it to a directory.  Point the REFERENCE environment variable to this directory before installing TIPP.  You can manually point TIPP to the reference directory by modifying the ~/.sepp/tipp.config file. 
   reference datasets: www.cs.utexas.edu/~phylo/software/sepp/tipp.zip


### 2. From Virtual Machine (VM)

VM Image (mostly for Windows users) is available for [download](http://www.cs.utexas.edu/~phylo/software/PASTA_TIPP_UPP.ova) (2.5 GB). Once the image is downloaded, you need to run it using a VM environment ([VirtualBox](https://www.virtualbox.org/) is a good option). After you install VirtualBox, you just need to use File/import to import the PASTA_TIPP_UPP.ova image that you have downloaded (If your machine has less than 3GB you might want to reduce the memory to something like 512MB). Once VM is imported, you can start it from the Virtualbox. If you are asked to login, the username and passwords are (username: phylolab, password: phylolab). TIPP and UPP are already installed on the VM machine, so you can simply proceed by opening a terminal and running it.

Note that we constantly update our software.  Before running the tutorial, it's best to grab
the most updated version of the software onto the VM machine.  This can be done by opening a terminal in the VM and typing the following commands:


```
cd ~/tools/sepp
git pull
```

If this command fails due to an error that the repository is corrupted, this can be fixed by typing the following series of commands from the SEPP directory:

```
rm -fr .git
git init
git remote add origin https://github.com/smirarab/sepp.git
git fetch
git reset --hard origin/master
```

Finally, if the BLAST environmental variable or the REFERENCE environmental variable
cannot be read by TIPP during the configuration, you can manually edit the ~/.sepp/tipp.config
file to point to the right location.  To do this, change:

```
[blast]
path=None

[reference]
path=None
```

```
[blast]
path=~/bin/

[reference]
path=~/testdata/tipp/
```

---------
Using TIPP
===

If your installation is successful, you should be able to run TIPP by running the following command from any location. Open up a terminal window and type: 

```
run_tipp.py -h
``` 

Running TIPP with the `-h` option produces the list of options available in TIPP. 

The general command for running TIPP for a specific pre-computed marker is:

```
run_tipp.py -R reference_marker -f fragment_file
```

Step 1: Running a test job
---

TIPP currently can only be run from the command line.  We have provided some test data files under the `test/` directory.  A good start is classifying reads from the pyrg gene, a smaller marker gene with only 65 sequences.

```
run_tipp.py -R pyrg -f test/unittest/data/mock/pyrg/pyrg.even.fas -d temp/ -p temp/tmp -o outads -P 30
```

This will run TIPP on the fragmentary sequences that have been binned to the pyrg gene.  
The main output of TIPP is a _classification.txt file that contains the classification of each read.  The classification consists of the name of the read, the NCBI taxonomic id of the classification,the rank of the classification, the name of the classification, and the support of the classification.

```
EAS25_26_1_15_381_1761_0_1,2157,Archaea,superkingdom,1.0000
EAS25_26_1_15_381_1761_0_1,1,root,root,1.0000
EAS25_26_1_15_381_1761_0_1,183925,Methanobacteria,class,1.0000
EAS25_26_1_15_381_1761_0_1,2172,Methanobrevibacter,genus,1.0000
EAS25_26_1_15_381_1761_0_1,28890,Euryarchaeota,phylum,1.0000
EAS25_26_1_15_381_1761_0_1,2158,Methanobacteriales,order,1.0000
EAS25_26_1_15_381_1761_0_1,2173,Methanobrevibacter smithii,species,1.0000
```

For example, the EAS25_26_1_15_381_1761_0_1 is classified as the species Methanobrevibacter smithii with 100% support.  In addition, TIPP outputs a .json file with the placements, created according to pplacer format. Please refer to pplacer website (currently http://matsen.github.com/pplacer/generated_rst/pplacer.html#json-format-specification) for more information on the format of the josn file. Also note that pplacer package provides a program called guppy that can read .json files and perform downstream steps such as visualization.

In addition to the .json file, TIPP outputs alignments of fragments to reference sets. There could be multiple alignment files created, each corresponding to a different placement subset. 

Step 2: Analyzing a metagenomic dataset
---

Step 1: Testing that SEPP is correctly installed:
---

SEPP can also be run using a configuration file. Sample configuration files and input files can be found under test/unittest/data/mock/. Change to that directory to run SEPP on the sample files. To run using command options, run

`run_sepp.py -t test.tree -a test.fasta -f test.fas -r test.RAxML_info -A 250 -P 250`

and to run using a configuration file, run

`python run_sepp.py -c sample.config`

The main output of SEPP is a .json file, created according to pplacer format. Please refer to pplacer website (currently http://matsen.github.com/pplacer/generated_rst/pplacer.html#json-format-specification) for more information on the format of the josn file. Also note that pplacer package provides a program called guppy that can read .json files and perform downstream steps such as visualization.

In addition to the .json file, SEPP outputs alignments of fragments to reference sets. There could be multiple alignment files created, each corresponding to a different placement subset. 

By setting SEPP_DEBUG environmental variable to `True`, you can instruct SEPP to output more information that can be helpful for debugging.  


Step 2: Testing that TIPP is correctly installed:
Go to the $SEPP/test/unittest/data/mock/pyrg directory, where $SEPP is the SEPP home directory
run 

`run_abundance.py -f pyrg.even.fas -c ~/.sepp/tipp.config -d out`

The system should create an output directory in the same directory labeled out, and within the directory, files called abundance.*.csv will be created.

### Running TIPP for profiling:

To run default TIPP on a set of metagenomic reads, run the command: 

`run_abundance.py -f fragment_file -c ~/.sepp/tipp.config -d output_directory`

The output will be tab delimited files that estimate the abundance at a given taxonomic level. For example, go to the $SEPP/test/unittest/data/mock/mixed directory, where $SEPP is the SEPP home directory and run 

`run_abundance.py -f facs_simhc.short.fas -c ~/.sepp/tipp.config -d out`

Below is an example of the abundance.species.csv output

`taxa    abundance
Agrobacterium tumefaciens       0.0370
Alcanivorax borkumensis 0.0370
Alcanivorax sp. DG881   0.0062
Anabaena variabilis     0.0432
Archaeoglobus fulgidus  0.0556
Bdellovibrio bacteriovorus      0.0247
Brucella abortus        0.0062
Burkholderia cenocepacia        0.0062
Campylobacter jejuni    0.0802
Candidatus Blochmannia floridanus       0.0432
Candidatus Phytoplasma australiense     0.0062
...
Pseudomonas fluorescens 0.0741
Pseudomonas putida      0.0062
Roseobacter denitrificans       0.0062
Streptomyces coelicolor 0.0494
Streptomyces ghanaensis 0.0062
Streptomyces griseoflavus       0.0062
Streptomyces scabiei    0.0123
Sulfolobus tokodaii     0.0494
Xanthomonas oryzae      0.0062
Yersinia pestis 0.0062
unclassified    0.0185`

This profile estimates that Pseudomonas fluorescens make up 7.4% of the species abundance, and 1.9% of the species could not be classified.


---------
Contact
===
Email: `ensemble-of-hmms@googlegroups.com` for all issues. 
