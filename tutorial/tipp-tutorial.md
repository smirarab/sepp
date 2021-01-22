TIPP Tutorial
=============

TIPP, which stands for *Taxonomic identification and phylogenetic profiling*, is a method for the following problems:

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

You have two options for installing TIPP. 

 - **Windows:** If you have a Windows machine, currently using the Virtual Machine (VM) image we provide is your only option.
 - **Linux:** and **MAC:** If you have Linux (or other \*nix systems) or MAC, you can still use VM, but downloading the code from github and installing it is what we strongly recommend. 
 
Installing From Source Code
---------------------------
The current version of TIPP has been developed and tested entirely on Linux and MAC. 
Windows does not work currently, and future versions may or may not support Windows. 

Before installing the software you need to make sure the following programs are installed on your machine.

- Python (version 2.7 or later, including python 3)
- Blast (version 2.10.1 or later)
- Java (version > 1.5 or later)

#### Step 1: Install SEPP
TIPP is a part of the SEPP distribution package. First download and install SEPP:

1. Open a terminal and create a directory where you want to keep SEPP. e.g. `mkdir ~/sepp-code`. Go to this directory. e.g. `cd ~/sepp-code`.
2. Clone the SEPP code repository from our [github repository](https://github.com/smirarab/sepp). For example you can use `git clone https://github.com/smirarab/sepp.git`.
If you don't have git, you can directly download a [zip file from the repository](https://github.com/smirarab/sepp/archive/master.zip)
and decompress it into your desired directory. 
3. Go to the sepp repository by using the command `cd sepp` (or `cd sepp-master` if you used the zip file instead of cloning the git repository). From the SEPP directory, configure SEPP by running the command: `python setup.py config` or alternatively, `python setup.py config -c`. The first command puts all the config files under your home directory while the second version uses the current path. 
4. Install SEPP by running the command: `sudo python setup.py install`. If you don't have root access, remove the `sudo` part and instead  use  `--user` option. Alternativley, you can `--prefix` to install in a different location, but that different location needs to be part of your `PYTHONPATH` environmental variable.


#### Step 2: Install TIPP

1. Download and decompress the reference dataset available at [https://obj.umiacs.umd.edu/tipp/tipp2-refpkg.tar.gz](https://obj.umiacs.umd.edu/tipp/tipp2-refpkg.tar.gz). For example, you could use the following commands:
  ```
  wget https://obj.umiacs.umd.edu/tipp/tipp2-refpkg.tar.gz
  tar xvzf tipp2-refpkg.tar.gz
  ```
2. Set the environment variable `REFERENCE` to point to the location of the reference directory. This can be performed using the command: 
  ```
  export REFERENCE=/PATH_TO_REFERENCE/tipp2-refpkg
  ```
3. Set the environment variable `BLAST` to point to the `blastn` binary.  This can be performed using the command:
  ```
  export BLAST=/PATH_TO_BLAST/blastn
  ```
4. Lastly, configure TIPP by running the following command from the SEPP directory:
  ```
  python setup.py tipp 
  ```
  or
  ```
  python setup.py tipp -c
  ```

**NOTE:** It's important that you use either use `-c` for both SEPP and TIPP or don't use it for both. 


**Common Problems:**

1. TIPP requires SEPP to be installed. If TIPP is not running, first check to see if TIPP was installed correctly.

2. TIPP relies on `blastn` for the binning of metagenomic reads, so BLAST needs to be downloaded and installed separately (learn more [here](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)). Then, point the `BLAST` environment variable to your installation of `blastn`. Alternatively, you can manually point TIPP to the `blastn` installation by modifying the `tipp.config` file. 
   
3. TIPP performs abundance profiling uses a set of 40 marker genes. This reference dataset needs to be downloaded separately from [here](https://obj.umiacs.umd.edu/tipp/tipp2-refpkg.tar.gz). Then, point the `REFERENCE` environment variable to the decompressed directory before installing TIPP. Alternatively, you can manually point TIPP to the reference dataset by modifying the `tipp.config` file.


Again, if the `BLAST` environment variable or the `REFERENCE` environmental variable cannot be read by TIPP during the configuration, you can manually edit the `tipp.config`
file to point to the right location. If you did NOT use the `-c` flag, then this file is in `~/.sepp/tipp.config`; otherwise, this file is in `/PATH_TO_SEPP_INSTALLATION/.sepp/tipp.config`. To edit the configuration file, change:
```
[blast]
path=None

[reference]
path=None
```
to
```
[blast]
path=/PATH_TO_BLAST/blastn

[reference]
path=/PATH_TO_REFERENCE/tipp2-refpkg
```

Installing from Virtual Machine (VM)
------------------------------------

VM Image (mostly for Windows users) is available for [download](https://drive.google.com/open?id=0B0lcoFFOYQf8U2NZV2Z2RmRaRjQ) (3 GB).
Once the image is downloaded, you need to run it using a VM environment ([VirtualBox](https://www.virtualbox.org/) is a good option).
After you install VirtualBox, you just need to use File/import to import the downloaded phylo.ova image that you have downloaded (If your machine has less than 3GB you might want to reduce the memory to something like 512MB). Once VM is imported, you can start it from the Virtualbox. If you are asked to login, the username and passwords are (username: phylolab, password: phylolab). TIPP is installed on the VM machine, so you can simply proceed by opening a terminal and running it.

Note that we constantly update our software.  Before running the tutorial, it's best to grab the most updated version of the software onto the VM machine. This can be done by opening a terminal in the VM and typing the following commands:
```
cd ~/tools/sepp
git pull
```

---------

Running TIPP
============

If your installation is successful, you should be able to run TIPP by running the following command from any location. Open up a terminal window and type: 

```
run_tipp.py -h
``` 

Running TIPP with the `-h` option produces the list of options available in TIPP. 

The general command for running TIPP for a specific pre-computed marker is:

```
run_tipp.py -R <reference_marker> -f <fragment_file>
```

Task 1: Performing read classification.
---------------------------------------

TIPP currently can only be run from the command line. We have provided some test data files under the `test` directory. A good start is classifying reads from the pyrg gene, a smaller marker gene with only 65 sequences.

```
run_tipp.py -R markers-v1/pyrg \
            -f test/unittest/data/mock/pyrg/pyrg.even.fas \
            -o output \
            -P 30
```

NOTE: The TIPP2 reference package includes three reference data sets:
- markers-v1 = reference dataset benchmarked in the 2014 paper
- markers-v2 = additional markers distributed with the 2014 paper
- markers-v3 = reference dataset created and benchmarked in the 2020 TIPP paper

This command runs TIPP on the fragmentary sequences that have been binned to the pyrg gene from the reference dataset benchmarked in the 2014 paper. In other words, the TIPP classification algorithm will use the pre-computed alignment and tree that has been estimated on the known bacterial pyrg genes.

The main output of TIPP is a `output_classification.txt` file that contains the classification of each read.  The classification consists of the name of the read, the NCBI taxonomic id of the classification,the rank of the classification, the name of the classification, and the support of the classification.

```
EAS25_26_1_15_381_1761_0_1,2157,Archaea,superkingdom,1.0000
EAS25_26_1_15_381_1761_0_1,1,root,root,1.0000
EAS25_26_1_15_381_1761_0_1,183925,Methanobacteria,class,1.0000
EAS25_26_1_15_381_1761_0_1,2172,Methanobrevibacter,genus,1.0000
EAS25_26_1_15_381_1761_0_1,28890,Euryarchaeota,phylum,1.0000
EAS25_26_1_15_381_1761_0_1,2158,Methanobacteriales,order,1.0000
EAS25_26_1_15_381_1761_0_1,2173,Methanobrevibacter smithii,species,1.0000
```

For example, the `EAS25_26_1_15_381_1761_0_1` is classified as the species Methanobrevibacter smithii with 100% support. By default, TIPP requires 95% placement support to classify a sequence. If we lower the classification threshold to 0 (using the `-pt` flag), we can see all possible classifications for a sequence.

```
run_tipp.py -R markers-v1/pyrg \
            -f test/unittest/data/mock/pyrg/pyrg.even.fas \
            -o output_lower_threshold \
            -P 30 \
            -pt 0.0
```

In addition, TIPP outputs a `.json` file with the placements, created according to pplacer format. Please refer to [pplacer website](http://matsen.github.com/pplacer/generated_rst/pplacer.html#json-format-specification) for more information on the format of the josn file. Also note that pplacer package provides a program called guppy that can read `.json` files and perform downstream steps such as visualization.
 
In addition to the `.json` file, TIPP outputs alignments of fragments to reference sets. There could be multiple alignment files created, each corresponding to a different placement subset. 

Task 2: Converting the result into an abundance profile.
--------------------------------------------------------
The classification file lists all possible classifications for a fragment, even if it has very low support.  In some situations, we only want the most highly supported classifications.  We can use a helper tool to convert the classification file into a profile.  

```
mkdir profile
run_tipp_tool.py -g markers-v1/pyrg \
                 -a profile \
                 -o profile \
                 -p markers-v1/pyrg \
                 -i output_classification.txt \
                 -t 0.95
```

This command will create taxonomic profiles (one for each taxonomic ranking) from the classification results.  Fragments will only be classified if they have at least 95% support for the classification.  Let's start by looking at the file labelled `pyrg.classification` in `profile/markers-v1`.

```
fragment        species genus   family  order   class   phylum
EAS25_26_1_100_940_776_0_1      2173    2172    2159    2158    183925  28890
EAS25_26_1_11_733_1260_0_2      2173    2172    2159    2158    183925  28890
EAS25_26_1_15_381_1761_0_1      2173    2172    2159    2158    183925  28890
EAS25_26_1_15_381_1761_0_2      NA      NA      2206    94695   224756  28890
```

This file lists the classification (shown as NCBI taxonomic ids) of each fragment at each of the taxonomic rankings.  If a fragment does not meet the support threshold (95% in this case), it will be left as unclassified (NA).  

Let's look at `abundance.species.csv`.  The file shows the abundance profiles for the species level.  The file shows that 80% of the reads belong to the species Methanobrevibacter smithii and 19% of the fragments were unclassified at the species level.
```
taxa    abundance
Methanobrevibacter smithii      0.7969
Methanococcus maripaludis       0.0156
unclassified    0.1875
```

Task 3: Running TIPP for abundance profiling.
---------------------------------------------

The previous example shows how to analyze a dataset when the fragments come from a specific gene.  When analyzing shotgun metagenomic reads, however, the reads originate from all across the genome.  Thus, we need to take a different approach for analyzing such a dataset.

TIPP comes with a collection of 30 single copy housekeeping genes which are used for abundance profiling of metagenomic reads.  The pipeline to analyze a metagenomic dataset is to first bin the reads to each of the marker genes, and then run TIPP individual on each of the individual marker genes.  We have simplified this process with a helper script `run_abundance.py`.

The general command for `run_abundance.py` is:

```
run_abundance.py -f <fragment file> \
                 -c <tipp config file> \
                 -d <output directory>
```

The input fragment files must be in FASTA or FASTQ formats with the following extenstions:
- `.fastq` or `.fq` for FASTQ files
- `.fasta`, `.fas`, `.fa`, `.fna` for FASTA files
The output will be tab delimited files that estimate the abundance at a given taxonomic level. 

By default, this command will use the reference dataset benchmarked in the 2020 paper. This is equivalent to specifying the option `-g markers-v3`; you can use the markers benchmarked in the 2014 paper by adding the `-G markers-v1` option.

For example, go to the `test/unittest/data/mock/mixed` directory 
```
cd test/unittest/data/mock/mixed
```
and run 
```
run_abundance.py -G markers-v1 \
                 -f facs_simhc.short.fas \
                 -c ~/.sepp/tipp.config \
                 -d out
```
**NOTE:** If you used the `-c` option when installing TIPP and SEPP, then instead of `~/.sepp/tipp.config`, you would use `/PATH_TO_SEPP_INSTALLATION/.sepp/tipp.config`.

Running this command creates an abundance profile by 
1. assigning the fragments to the marker genes using BLAST (while orienting and trimming reads as discussed in the TIPP2 paper),
2. running TIPP on all reads (using the marker that they were assigned to as a reference), and
3. aggregating the classifications to create an abundance profile at each taxonomic level.

The `markers` directory contains the individual results for each marker.

Below is an example of the `abundance.species.csv` output from the run.

```
taxa	abundance
Agrobacterium tumefaciens	0.0400
Alcanivorax borkumensis	0.0467
Alcanivorax sp. DG881	0.0067
Anabaena variabilis	0.0400
Archaeoglobus fulgidus	0.0533
Bdellovibrio bacteriovorus	0.0333
Campylobacter jejuni	0.0733
Candidatus Blochmannia floridanus	0.0600
Clostridium acetobutylicum	0.0600
Escherichia coli	0.0333
Francisella tularensis	0.0267
Haemophilus influenzae	0.0067
Lactococcus lactis	0.0200
Methanoculleus marisnigri	0.0267
Nitrosomonas europaea	0.0400
Pasteurella multocida	0.0533
Pseudomonas aeruginosa	0.0400
Pseudomonas entomophila	0.0267
Pseudomonas fluorescens	0.0533
Streptomyces coelicolor	0.0333
Sulfolobus tokodaii	0.0400
unclassified	0.1867
```

This profile estimates that `Pseudomonas fluorescens` make up 5.3% of the species abundance and that 18.7% of the fragments could not be classified even though they may be classified at other levels.

Task 4: Looking a reference dataset.
------------------------------------

Let's take a look at the files within a reference dataset, starting with the taxonomy file. Go to the `16S_bacteria.refpkg` directory (can be found in the `16S-bacteria-v1` directory of the `tipp2-refpkg.tar.gz` archive).

```
head -n2 all_taxon.taxonomy
"tax_id","parent_id","rank","tax_name","root","below_root","below_below_root","superkingdom","below_superkingdom","below_below_superkingdom","below_below_below_superkingdom","superphylum","phylum","below_phylum","below_below_phylum","subphylum","class","below_class","below_below_class","below_below_below_class","subclass","order","below_order","below_below_order","suborder","below_suborder","family","below_family","below_below_family","below_below_below_family","subfamily","tribe","genus","below_genus","subgenus","species_group","species_subgroup","species"
"1","1","root","root","1","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""
```

This file represents a taxonomy as a comma delimited file. The first line is the header that describes each column. The important files are the unique id of each clade, the clade name, rank, and the parent of the clade.

Next is the backbone alignment and tree files (`sate.fasta` and `sate.taxonomy`). These are the full-length sequences from the known organisms. The sequences in this file are mapped to the taxonomy through the `species.mapping` file shown below.

```
seqname,tax_id
S000438419,10
S000539682,100
```

Each sequence name is mapped to the unique id in the taxonomy file.

Finally, in order to find the best placement, we need the model parameters of the taxonomic tree. This can be generated by RAxML using the `-f e` or the `-g` option.  

Thus, specialized marker datasets can be generated for any organisms, not just bacteria, by providing these files.

Task 5: Running a 16S amplicon analysis.
----------------------------------------

Finally, we have included a 16S reference marker gene that can be used to analyze 16S amplicon data. Below is an example of running TIPP on 16S amplicon data.
```
run_tipp.py -R 16S-bacteria-v1/16S_bacteria \
            -f test/unittest/data/mock/16S_bacteria/human_gut_16S.fas \
            -o 16s \
            -A 1000 \
            -P 1000
```

As in the previous example, you can convert the classification results into a more easily digestible format using the `run_tipp_tool.py` script.
```
run_tipp_tool.py -g 16S-bacteria-v1/16S_bacteria \
                 -a profile_16s \
                 -o profile_16s \
                 -p 16S-bacteria-v1_16S_bacteria \
                 -i 16s_classification.txt \
                 -t 0.95
```

---------

Bugs and Errors
===============
TIPP is under active research development at UIUC by the Warnow Lab. Please report any errors to Tandy Warnow (warnow@illinois.edu) and Siavash Mirarab (smirarab@ucsd.edu).
