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
export REFERENCE=/PATH/TO/REFERENCE
```

4. Set the environment variable BLAST to point to blastn.  This can be performed using:

```
export BLAST=/PATH/TO/DIRECTORY/blastn
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
path=~/bin/blastn

[reference]
path=~/testdata/tipp/
```

### 2. From Virtual Machine (VM)

VM Image (mostly for Windows users) is available for [download](http://www.cs.utexas.edu/~phylo/software/WARNOW_LAB.ova) (2.5 GB). Once the image is downloaded, you need to run it using a VM environment ([VirtualBox](https://www.virtualbox.org/) is a good option). After you install VirtualBox, you just need to use File/import to import the WARNOW_LAB.ova image that you have downloaded (If your machine has less than 3GB you might want to reduce the memory to something like 512MB). Once VM is imported, you can start it from the Virtualbox. If you are asked to login, the username and passwords are (username: phylolab, password: phylolab). TIPP is installed on the VM machine, so you can simply proceed by opening a terminal and running it.

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
run_tipp.py -R pyrg -f test/unittest/data/mock/pyrg/pyrg.even.fas  -o output -P 30
```

This will run TIPP on the fragmentary sequences that have been binned to the pyrg gene.  This will use the pre-computed alignment and tree that has been estimated on the known bacterial pyrg genes.

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

For example, the EAS25_26_1_15_381_1761_0_1 is classified as the species Methanobrevibacter smithii with 100% support.  By default, TIPP requires
95% placement support to classify a sequence.  If we lower the classification
threshold, we can see all possible classifications for a sequence.

```
run_tipp.py -R pyrg -f test/unittest/data/mock/pyrg/pyrg.even.fas -o lower_threshold -P 30 -pt 0.0
```

In addition, TIPP outputs a .json file with the placements, created according to pplacer format. Please refer to pplacer website (currently http://matsen.github.com/pplacer/generated_rst/pplacer.html#json-format-specification) for more information on the format of the josn file. Also note that pplacer package provides a program called guppy that can read .json files and perform downstream steps such as visualization.
 
In addition to the .json file, TIPP outputs alignments of fragments to reference sets. There could be multiple alignment files created, each corresponding to a different placement subset. 

Step 2: Converting the result into an abundance profile
---
The classification file lists all possible classifications for a fragment, even if it has very low support.  In some situations, we only want the most highly supported classifications.  We can use a helper tool to convert the classification file into a profile.  

```
mkdir profile
run_tipp_tool.py -g pyrg -a profile -o profile -p pyrg -i output_classification.txt -t 0.95
```

This command will create taxonomic profiles (one for each taxonomic ranking) from the classification results.  Fragments will only be classified if they have at least 95% support for the classification.  Let's start by looking at the file labelled pyrg.classification

```
fragment        species genus   family  order   class   phylum
EAS25_26_1_100_940_776_0_1      2173    2172    2159    2158    183925  28890
EAS25_26_1_11_733_1260_0_2      2173    2172    2159    2158    183925  28890
EAS25_26_1_15_381_1761_0_1      2173    2172    2159    2158    183925  28890
EAS25_26_1_15_381_1761_0_2      NA      NA      2206    94695   224756  28890
```

This file lists the classification (shown as NCBI taxonomic ids) of each fragment at each of the taxonomic rankings.  If a fragment does not meet the support threshold (95% in this case), it will be left as unclassified (NA).  

Let's look at abundance.species.csv.  The file shows the abundance profiles for the species level.  The file shows that 80% of the reads belong to the species Methanobrevibacter smithii and 19% of the fragments were unclassified at the species level.
```
taxa    abundance
Methanobrevibacter smithii      0.7969
Methanococcus maripaludis       0.0156
unclassified    0.1875
```

Step 3: Running TIPP for profiling:
---

The previous example shows how to analyze a dataset when the fragments come from a specific gene.  When analyzing shotgun metagenomic reads, however, the reads originate from all across the genome.  Thus, we need to take a different approach for analyzing such a dataset.

TIPP comes with a collection of 30 single copy housekeeping genes which are used for abundance profiling of metagenomic reads.  The pipeline to analyze a metagenomic dataset is to first bin the reads to each of the marker genes, and then run TIPP individual on each of the individual marker genes.  We have simplified this process with a helper script run_abundance.py.

The general command for run_abudnance.py is:

```
run_abundance.py -f fragment_file -c ~/.sepp/tipp.config -d output_directory
```

The output will be tab delimited files that estimate the abundance at a given taxonomic level. For example, go to the test/unittest/data/mock/mixed directory and run 

```
run_abundance.py -f facs_simhc.short.fas -c ~/.sepp/tipp.config -d out
```

Running this command bin the fragments to the marker genes using BLAST.  Once the sequences have been binned to the marker genes, the direction of the sequences will be estimated by aligning the forward and reverse complemented sequences against the HMM.  The orientation with the best HMM score is selected.  Finally, TIPP is run on the sequence to classify it.  The markers directory contains the individual results for each marker.

Below is an example of the abundance.species.csv output from the run.

```
taxa    abundance
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
unclassified    0.0185
```

This profile estimates that Pseudomonas fluorescens make up 7.4% of the species abundance, and 1.9% of the species could not be classified.

Step 4: Looking a reference dataset:
---

Let's take a look at the files within a reference dataset, starting with the taxonomy file.  Go to the 16S_bacteria.refpkg directory (can be found in the tipp.zip archive)

```
head all_taxon.taxonomy

"tax_id","parent_id","rank","tax_name","root","below_root","below_below_root","superkingdom","below_superkingdom","below_below_superkingdom","below_below_below_superkingdom","superphylum","phylum","below_phylum","below_below_phylum","subphylum","class","below_class","below_below_class","below_below_below_class","subclass","order","below_order","below_below_order","suborder","below_suborder","family","below_family","below_below_family","below_below_below_family","subfamily","tribe","genus","below_genus","subgenus","species_group","species_subgroup","species"
"1","1","root","root","1","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""
...
```

This file represents a taxonomy as a comma delimited file.  The first line is the header that describes each column.  The important files are the unique id of each clade, the clade name, rank, and the parent of the clade.

Next is the backbone alignment and tree files (sate.fasta and sate.taxonomy).  These are the full-length sequences from the known organisms.  The sequences in this file are mapped to the taxonomy through the species.mapping file shown below.

```
seqname,tax_id
S000438419,10
S000539682,100
```

Each sequence name is mapped to the unique id in the taxonomy file.

Finally, in order to find the best placement, we need the model parameters of the taxonomic tree.  This can be generated by RAxML using the `-f e` option.  

Thus, specialized marker datasets can be generated for any organisms, not just bacteria, by providing these files.

Step 5: 16S amplicon analysis:
---

Finally, we have included a 16S reference marker gene that can be used to analyze 16S amplicon data.  Below is an example of running TIPP on 16S amplicon data.

```
run_tipp.py -R 16S_bacteria -f test/unittest/data/mock/16s_bacteria/human_gut_16S.fas -o 16s -A 1000 -P 1000
```

As in the previous example, you can convert the classification results into a more easily digestible format using the run_tipp_tool.py script:

```
run_tipp_tool.py -g 16_bacteria -a profile -o -p 16_bacteria -i 16s_classification.txt -t 0.95
```

---------
Contact
===
Post all questions, comments, requests to: https://groups.google.com/forum/#!forum/ensemble-of-hmms


