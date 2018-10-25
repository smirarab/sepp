Introduction to SEPP 
=================================

SEPP [1] stands for SATé-Enabled Phylogenetic
Placement and addresses the problem of phylogenetic placement for
meta-genomic short reads. More precisely, SEPP addresses the
following problem.

* **Input**:  
    i) *backbone* tree `T` and backbone alignment `A` for a set of
    full-length gene sequences
    ii) the set `X` of fragmentary sequences from
    the same gene as the backbone
* **Output**: the placement of each fragment in `X` onto the tree T and the alignment of
    all fragment in `X` to the alignment `A`.

Phylogenetic placement adds unknown (short) fragments into a phylogenetic
tree and hence helps identifying species included in a metagenomic
dataset. Phylogenetic placement involves two steps: alignment of short
fragments to full length sequence alignment `A` (*backbone* alignment) and then placement of aligned short
reads on the fixed tree `T` (*backbone* tree).

SEPP operates by using a divide-and-conquer strategy adopted
from SATé [2] and SATé-II [3] to improve the alignment of
fragments to the backbone alignment (produced by running
HMMER [4]). It then places each fragment into the
user-provided tree using pplacer [5]. 

Our studies shows
that SEPP provides improved accuracy for quickly evolving
genes as compared to other methods. For more information see the [paper](http://www.ncbi.nlm.nih.gov/pubmed/22174280). 



Installing SEPP 
============================

First you are going to setup SEPP on your machines.
SEPP currently runs only on Linux and Mac. Those running a
Windows need to either use cygwin (<http://www.cygwin.com/>), or a
Ubuntu virtual machine image we have created.

**No setup-up Greengrenes version:** on [this page](../sepp-package/) you can find a version of SEPP that doesn't require installation and comes prepackaged ready to add taxa to the greenness dataset. 

Installing on a Linux/MAC Machine
-----------------------------

To download and install SEPP , follow these steps:

1.  Obtain the software code. There are two options
	-  From the zip file:
		1. Download software from
    <https://github.com/smirarab/sepp/zipball/master>
    	2. Copy the archive to your favorite location (we will use
    `~/sepp` in this tutorial)
    	3. Unpack the zip file using your favorite software.
	-  Clone the git repository <https://github.com/smirarab/sepp>
	(we will assume SEPP is clone to `~/sepp` in this tutorial)

2.  Open a terminal and change into the unpacked directory
    (`cd ~/sepp`).

3.  Run the following command:
    
    `
        python setup.py config
    `   
    
    Note that this step needs to be run by each individual user. This
    step creates a `.sepp/main.config` file in
    your home directory and copies the bundled tool under the same
    directory. This file contains defaults settings of SEPP, and can be
    modified by users if needed.

    If you wish not to use the home directory, you can use

    `python setup.py config -c`

    instead. This will ensure that the main config file and binary files
    are written to your current directory instead of the home.

4.  In the new directory, there should be a
    `setup.py` file. Run the following command:

    `sudo python setup.py install`

    If you don’t have root access, instead use the following command (or
    use `–prefix` option; see the `README` file for details)

    `python setup.py install --user`

Refer to the [README](../README.SEPP.md) file for more information regarding the
installation and solutions to common installation problems.

Using Virtual Machine
---------------------

An Ubuntu VM image with SEPP installed on it is available
for download at <https://drive.google.com/file/d/0B0lcoFFOYQf8U2NZV2Z2RmRaRjQ/view?usp=sharing>.

If you were unable to install SEPP on your machine, you can
use this VM image. You first need to copy the VM image to your machine.
Then, open the VM image in your favorite VM software (e.g. VirtualBox)
and start the VM. Once your VM starts, you can run SEPP from
a terminal. The username and password for the virtual machine are
“osboxes” and “osboxes.org”.

Running SEPP 
=========================

SEPP is currently available only as a command line tool, and
so the tutorial is based on this command line usage.

In this section we will run SEPP on a small sample dataset
provided with the software. The sample dataset consists of a
SATé backbone alignment and tree on the “pyrg” marker gene
with only 65 sequences (previously studied in Metaphyler [6]). The
fragments are from a WGS sample of a mock community created by the NIH
Human Microbiome Project (<http://www.hmpdacc.org/HMMC/>). The fragment
file we provide includes only 106 fragments that we found to possibly
belong to “pyrg” marker (based on hits to its HMMER profile).

Run SEPP with `-h` option to see the help 
----------------------------------------------------

-   Make sure you have a terminal open.

-   Make a directory where you will run SEPP and `cd` into it; e.g.
    
    ```    
     mkdir seppRuns
     cd seppRuns
    ```

-   run SEPP with the `-h` option to see a help:

    ```
    run_sepp.py -h
    ````

Sample Datasets: default parameters
------------------------------------

-   Copy test datasets into your directory. 
    Test datasets are part of the distribution and can be found under
    `/test/unittest/data` (use the directory
    where you unpacked SEPP instead of
    `~/sepp`). For example, on Unix run:
    
    ```    
    cp -R ~/sepp/test/unittest/data/mock .
    ```
    
    Note that on the VM, the files to copy are found in
    `~/testdata/sepp/`.

-   Execute the following command to run SEPP on a sample
    biological dataset.

    ```
    run_sepp.py -t mock/pyrg/sate.tre -r mock/pyrg/sate.tre.RAxML_info -a mock/pyrg/sate.fasta -f mock/pyrg/pyrg.even.fas
    ```

    (SEPP should finish running quickly; 10 seconds on my machine)

This runs SEPP on the given example dataset. All the options
provided to SEPP in this example run are mandatory. As you
can see, there are few inputs required to run SEPP . The
following describes the minimum input of SEPP:

* **Backbone tree**:  this is the tree on which SEPP places short fragments.
    This tree should be a binary maximum likelihood (ML) tree in newick
    format; we therefore recommend you estimate the ML tree using
    RAxML [7] or Phyml [8] on the backbone
    alignment. The input tree is given to
    SEPP using `-t` option.
* **Backbone Alignment**: this is a multiple sequence alignment of full length sequences for
    some gene. These sequences need to be for the same gene as the
    fragmentary sequences, as phylogenetic placement only makes sense
    in this case. The backbone alignment needs to be highly accurate,
    since it determines the backbone tree, and the use of RAxML to
    estimate an ML tree on the backbone alignment may help with the
    accuracy. If you obtain a backbone tree in some other way, you
    should ensure that the backbone tree and alignment are on the same
    exact set of taxa. The backbone alignment is provided to
    SEPP using `-a` option, and should be in the Fasta format.
    You can use refseq, available at
    <http://www.ncbi.nlm.nih.gov/RefSeq/> to convert between different
    alignment formats.
* **Stats or info file**
:   this is the info file generated by RAxML (or phyml) when it computes
    the ML tree on the backbone alignment (i.e., the backbone tree).
    This file is required by pplacer (run internally by
    SEPP ) in order to avoid re-estimation of ML parameters.
    To be able to use SEPP you need to make sure you keep
    your info file when you are generating the backbone tree. If you do
    not have the info file (or if you used some other software programs,
    such as PASTA [9], to produce the backbone tree), you can use
    RAxML’s `-f e` option to quickly estimate
    the model parameters (including branch lengths) on your backbone
    tree topology (see [this section](#sec:backbone)). The RAxML info file
    should be provided to SEPP using `-r` option.
* **Fragments file**:   this is a Fasta file containing the actual short fragments that are
    going to be placed. Fragments file should be given to
    SEPP using `-f` option.

### 10% rule

Recall that SEPP operates by dividing the set of taxa in the
backbone alignment and tree into alignment subsets and placement
subsets. In our test run, we did not explicitly set the maximum subset
sizes for alignment and placement subsets. The choice of these two
parameters affects accuracy on the one hand and computational resources
on the other hand, as explained below.

-   Decreasing alignment sizes should increase (and have increased in
    our experience) the accuracy of SEPP . On the other
    hand, smaller subsets increase the running time. This is because
    SEPP needs to score each fragment against all subsets
    independently, and therefore increasing the number of subsets adds
    to the running time. Note that extremely small subsets (i.e. less
    than 10 taxa) have not been tested and are not
    recommended.

-   Increasing placement sizes should result in better accuracy in
    general (although there could be exceptions). If your placement tree
    is very large (thousands or tens of thousands of leaves), the memory
    requirement of pplacer , and hence of SEPP ,
    increases dramatically. Reducing the placement size reduces the
    memory footprint, and hence enables placement on larger trees given
    a fixed amount of memory available. This would be one of the main
    motivations to reduce placement subset size. Reducing the placement
    subset *can* result in reduced running time as well,
    especially if your placement tree has thousands of taxa. For smaller
    trees, the effect of the placement size on the running time is not
    easily predicted, and is practically of less interest.

By default, when alignment and placement subset sizes are not explicitly
specified by user, SEPP uses what we call the “10% rule”
to automatically set those parameters. 10% rule specifies that alignment
and placement subset sizes should be both set to 10% of the number of
full length sequences (i.e. number of leaves in the backbone tree). The
10% rule is just a heuristic setting we have found empirically to give a
reasonable tradeoff *in general* between accuracy and
computational requirements on the datasets we have tried. Users are
encouraged to change subsets sizes based on their available
computational resources, and the desired accuracy, according to the
guidelines outlined above. For example, in our prepackaged SEPP
for the Greengenes dataset, we are using placement sizes that
 are closer to 2% and alignment subset sizes close to 0.5% to allow scaling to reference 
 datasets with many hundreds of thousands of tips. 

### Specifying subset sizes using `-P`, `-A`, `-M`, and `-F` options

Imagine we cannot wait 2 minutes to get results on our test dataset. We
are going to increase the alignment subset so that SEPP runs
faster. The test dataset included 65 full length sequences, and hence
the 10% rule amounts to alignment and placement subsets of maximum size
7. Since our toy example of 65 sequences is very small, it makes sense
to increase the alignment size to a larger number (e.g. 10), and
placement size to the entire dataset (i.e. 65). The maximum alignment
and placement subset sizes are controlled with `-A` and `-P` options,
respectively.

-   Execute the following command to run SEPP with `-A=10`
    and `-P=65`:

    ```
    run_sepp.py -t mock/pyrg/sate.tre -r mock/pyrg/sate.tre.RAxML_info -a mock/pyrg/sate.fasta -f mock/pyrg/pyrg.even.fas -A 10 -P 65 -o run2.A10P65
    ```

    (SEPP should finish in less than a minute.)

In the above run note the `-o` option. This option controls the prefix of
the output files generated by SEPP . We have not looked at
SEPP output yet, and we will do so in a moment. For now,
just be aware that SEPP generates a bunch of output files,
and prefixes those with a given string, which defaults to “output”.
These outputs are by default generated in the current directory, but
that can be changed using the `-d` option. Had we not changed the output
prefix, SEPP would have refused to run to avoid overwriting
results of your previous run saved with “output” prefix. Try this,
and you would get the following error:

`
Output directory [a_drecotry] already contains files with prefix [output]... Terminating to avoid loss of existing files.
`


In addition to these options, you can also use the `-F` option to 
control the size of the fragment chunks used internally in SEPP. 
Lowering the size of the fragment chunks
helps control the amount of memory used by SEPP. 

#### Diameter-based decomposition
Instead of dividing until you get to a subset size, you can use the `-M` option to specify a maximum diameter for *alignment* subsets. Thus, for example, `-M 0.5` will tell SEPP to decompose until it reaches alignment subsets with a diameter of no more than 0.5. 
The diameter of a subset is the maximum total branch length on the path between any two taxa in that subset. 

Several points should be emphasized. 

* When you specify `-M`, it makes sense to set the alignment subset size `-A` to be equal to `-P` so that alignment subset sizes do not impact the stopping criteria. Otherwise, if both `-M` and `-A` are in effect, SEPP tries to honor both and may wind up honoring neither (when impossible to honor both). 
* Note that `-M` only impacts alignment subsets and not placement subsets. 
* In addition to these `-A` and `-M` options, there is another value that can impact the alignment subset decomposition: minimum subset size. 
  SEPP imposes a minimum subset size of 2 sequences by default. 
  This minimum subset size, which is always honored, can be changed by editing the main or a custom config file (both described below) and editing the following line:
  ```
  [exhaustive]
  minsubsetsize = 2
  ```
  For larger trees, increasing this minimum to something like 20 makes sense. 
 * Similar to the minimum alignment subset size, there is a minimum placement subset size as well. 
  This, however, is expressed in terms of the fraction of the maximum placement subset (i.e., `-P`). 
  By default, this value is set to 1/4 of the maximum placement subset size, but the value can be adjusted in the config file by editing
  ```
  [exhaustive]
  placementminsubsetsizefacotr = 4
   ```
* **Important:** Currently, when using the `-M` option, you need to also use `-S` followed by either `centroid` or `midpoint`. Using `centroid` would result in the normal SEPP decomposition. Using the `midpoint` option instructs SEPP to also take into account branch lengths when dividing the tree. With `midpoint`, each time the tree is being decomposed, SEPP will find the midpoint branch instead of the centroid branch. However, if the midpoint branch results in subsets that are smaller than the minimum subset size, then SEPP will revert back to using the centroid edge decomposition for that step (but continues trying the midpoint for all other decompositions).  
 Note that specifying `-S` option changes the decomposition strategy both for alignment and placement subsets. Also note that `-S` does not have to be used with `-M` and can be used with `-A`. 


Running SEPP on a larger dataset
--------------------------------------------

We are now going to run SEPP on a larger dataset. Similar to
the small “pyrg” dataset, our larger dataset is on fragments from a WGS
sample of the HMMC mock community (<http://www.hmpdacc.org/HMMC/>), but
is based on a much larger marker gene called “rpsS” (obtained from
MetaPhyler as well). The backbone alignment and tree are again estimated using
SATé. This dataset is available under
`mock/rpsS` folder and has 1,277 full length
sequences and 2101 fragments. We are going to start a
SEPP run on this dataset. Run the following command:

```
run_sepp.py -t mock/rpsS/sate.tre -r mock/rpsS/sate.tre.RAxML_info -a mock/rpsS/sate.fasta -f mock/rpsS/rpsS.even.fas -o rpsS.out.default
```

While this is running (30 seconds on my laptop), we are going to look at SEPP outputs.

SEPP output 
========================

SEPP has two outputs: an alignment of fragments to full
length sequences (or subsets of the full length sequences), and
placements of fragments on the given backbone tree. When
SEPP finishes running, it generates two or more output
files, all prefixed with a string given using `-o` option, and placed in a
directory given using `-d` option. One of these output files is
`[prefix]_placement.json` file. This is the results of
the placement algorithm, and is in a format devised by the
pplacer software. This `.json` file is a human-readable text
file. However in most cases you want to look at those results in a
visualization tool. 
The pplacer toolkit comes with a suit of software
tools called  guppy. The guppy tool reads a
`.json` file, and among other things, produces
nice visualizations of the results. We are going to first manually look
at the plain `.json` file to understand its
content, and then will use guppy to visualize results.

-   Open the file called
    `output_placement.json` using a text
    editor.

-   Notice at the top of the file there is a newick tree with edges
    labeled with numbers inside brackets.

-   Following the tree, placement results are given for each fragment.
    Everything between “{” and “}” describes the placement of a
    single fragment. Each fragment can have multiple placements, with
    different likelihoods. Each line under the “p” attribute
    indicates one placement of the fragment.     The first value gives the
    edge label, the second value gives the log likelihood, the third
    value gives probability of that placement, and the final two values
    give the position on the edge where the fragment is placed, and the
    length (as estimated using the maximum likelihood calculation in
    pplacer) of the pendant edge for the fragment.

Next we will use guppy to turn this text file to a
visualization of the results. 
SEPP comes prepackaged with guppy included. Where guppy is to be found depends
on your installation. Possible locations include: 
`~/.sepp/bundled-v*/guppy` or if you used `-c` when installing, it could
be under `[installation directory]/.sepp/bundled-v*/guppy` where `[installation directionry]`
is where you installed SEPP from. You need to find your guppy location. We
assume it's on `~/.sepp/bundled-v3.2/guppy`.

Issue the following command to generate a
tree that has fragments placed on the backbone tree based only on
*the best placement* of each fragment.

```
    ~/.sepp/bundled-v3.2/guppy tog --xml output_placement.json
```

This command generates a new file called
`out_placement.tog.xml`. This is an XML file in
NexML format (<http://www.nexml.org/>) and can be opened with
[Archaeopteryx](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx).
Run Archaeopteryx and use `File/Read Tree from
File` to open `out_placement.tog.xml`. 
Select “colorize branch” checkbox on the right hand side panel (on newer versions of
Archaeopteryx , the option may be called something else, like “visual styles/branch colors”). 
The white branches represent the backbone tree, and branches
in red correspond to the maximum likelihood placement of fragments on
the backbone tree.

By now the run we started on the larger dataset (rpsS gene) should have
finished as well (it takes about 5 minutes). Use guppy to
visualize the results of that run as well.

Please refer to pplacer and guppy documentation
at <http://matsen.github.com/pplacer/generated_rst/pplacer.html> for
more information about `.json` format and
visualization options available in guppy. In addition to
the placement file, an extended alignment file is also generated. This
extended alignment shows how fragments are aligned to the backbone
alignment. The extended alignment is a simple Fasta file, and can be
viewed in any alignment visualization tool (e.g. JalView available at
<http://www.jalview.org/>).



SEPP Obtaining Backbone Alignment and Trees <a name="sec:backbone"></a>
=======


In all examples given above, backbone alignment and trees were already
estimated. Our suggested way of obtaining backbone alignment and trees
is through PASTA, which simultaneously estimates both an
alignment and a tree based on unaligned full length sequences. Please
refer to PASTA documentation for more information on running
PASTA.

In addition to an alignment and a tree (obtained from
PASTA or otherwise), we also need to have a RAxML info file.
If your backbone tree is estimated using RAxML you already have the info
file. Otherwise, you can optimize model parameters on your backbone tree
by running the following:

```
raxml -g [backbone tree] -s [backbone alignment] -m GTRGAMMA -n some_name_you_chooose -p some_random_number
```

This will optimize GTRGAMMA model parameters on your input
alignment/tree pair and will generate a info file
(`RAxML_info.some_name_you_chooose`), that can be used with
SEPP .

**Note:** With newer versions of RAxML, you may have to manually edit the RAxML info file to fix some formatting changes. 
For example, with version 8.0.22, we had to manually remove a line that read: ` Partition: 0 with name: No Name Provided` and
only then, the info file was recognized by pplacer.
[This script](./sepp-package/buildref/reformat-info.py) may become useful for reformatting if you ran RAxML with `-f e` to get your info file and you are using DNA sequences.

For example, in the test directory, you can go to `mock/pyrg` and run:

```
raxmlHPC-SSE3 -g sate.tre -s sate.fasta  -m GTRGAMMA -n newraxmlinf -p 24222
```
which produced a file called `RAxML_info.newraxmlinf`, which, after removing the line mentioned above, can be used as input to the `-r` option in SEPP. 


SEPP On Greengenes
==============================

We have built a stand-alone version of SEPP to place 16S fragments on the greengenes dataset.

1. Refer to [this page](../sepp-package/) for setting up (super-easy)

2. One installed, placing on the greengenes, is a one line simple command. 
   To place a test fragmentary file:

   ```
   cd [sepp-installation-location]/sepp-package
   ./run-sepp.sh test.frag test-gg
   ```

This will create the `test-gg_placement.tog.tre` file, which you can see using Archaeopteryx


SEPP Miscellaneous 
===============================

The following are some other points that are worth mentioning and
testing.

-   By default, the input is assumed to be DNA. To analyze amino acid
    datasets use `-m` option (i.e.
    `-m amino` for amino acid and `-m
    rna` for RNA).

-   By default SEPP tries to use all available cores on your
    machine in each run. If you run multiple instances of
    SEPP simultaneously, or if you want it to use fewer
    cores, be sure to set the number of cpus used by
    SEPP using `-x` option. For
    example the following runs SEPP on the large dataset but
    with only 2 cores used:

    ```
    run_sepp.py -t mock/rpsS/sate.tre -r mock/rpsS/sate.tre.RAxML_info -a mock/rpsS/sate.fasta -f mock/rpsS/rpsS.even.fas -x 2
    ```

-   In addition to commandline, SEPP can be controlled
    through a configuration file, passed to SEPP using
    `-c` option. For example, to run SEPP using `config.run1`
    configuration file, use:

    ```
    run_sepp.py -c config.run1
    ```

    -   Commandline options can be specified in the configuration file under
    the section `command line`. For specifying options from command line, 
    you need to use their long format name (as show in the SEPP help invoked by
    `-h`). For example, to set input to “rpsS”
    dataset and the alignment size to 100 and number of cpus to 3 use:

    ```
	[commandline]
	alignmentSize = 100
	tree= mock/rpsS/sate.tre
	raxml = mock/rpsS/sate.tre.RAxML_info
	alignment = mock/rpsS/sate.fasta
	fragment = mock/rpsS/rpsS.even.fas
	cpu = 3
    output = config
    ```

    -   Some extra options not available in the commandline can be
    configured in various sections of the configuration file. For
    example,

        ```
        path = /some/path
        ```

        tells SEPP that pplacer binaries can be
    found under `/some/path` instead of the default location.

    -   An example config file is available as part of the distribution
    under the test directory
    (`test/unittest/data/simulated/sample.config`).
    -   A main configuration file under
    `~/.sepp/main.config` is used to store
    some basic configurations such as the location of extra programs,
    etc. If you used `-c` when installing, this file would be instead
    located under your installation directory instead of your home. 
    -   When conflicting options are given, precedence is with those
    provided through commandline, then those specified in config file
    provided using `-c` option and finally those
    specified in the main config file. To test running from the config
    file, `cd simulated` directory and run
    
        ```
        run_sepp.py -c test.config
        ```

-   Currently SEPP has a built-in checkpointing
    functionality. By default, this functionality is turned off, but can
    be turned on using `-cp [checkpoint file
    name]`, and the time interval between two consecutive
    checkpoints can be adjusted using `-cpi [checkpoint
    frequency in seconds]`. Please note that checkpointing
    options have been tested only lightly and might have unknown issues.

-   For datasets containing both fragmentary and full-legnth sequences,
    the sequences much be separated out into a backbone set on the
    full-length sequences and a query set on the fragmentary sequences.
    From there, an alignment and tree can be estimated on the backbone
    set, and SEPP can be used to insert the query sequences back into
    the backbone tree. Run `split_sequences.py -i mock/pyrg/data/mixed.fas -o split -t 150`

    This will split the sequences into two files, `split.full.fas` (all
    sequences longer than 150 bps) and `split.frag.fas` (all sequences
    shorter than or equal to 150 bps). From here, PASTA can be used to
    estimate the backbone alignment and RAxML can be used to estimate
    the backbone tree. The query sequences can be inserted into the
    tree.

### Refrences

[1]: Mirarab, S., Nguyen, N. & Warnow, T. SEPP: SATé-Enabled Phylogenetic Placement. Pacific Symp. Biocomput. 247–58 (2012).

[2]: Liu, K., Raghavan, S., Nelesen, S. M., Linder, C. R. & Warnow, T. Rapid and Accurate Large-Scale Coestimation of Sequence Alignments and Phylogenetic Trees. Science (80-. ). 324, 1561–1564 (2009).

[3]: Liu, K. et al. SATe-II: Very Fast and Accurate Simultaneous Estimation of Multiple Sequence Alignments and Phylogenetic Trees. Syst. Biol. 61, 90–106 (2011).

[4]: Eddy, S. R. A new generation of homology search tools based on probabilistic inference. Genome Inf. 23, 205–211 (2009).

[5]: Matsen, F. A., Kodner, R. B. & Armbrust, E. V. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC Bioinformatics 11, 538 (2010).
[6]:  Liu, B., Gibbons, T., Ghodsi, M. & Pop, M. MetaPhyler: Taxonomic profiling for metagenomic sequences. in Bioinformatics and Biomedicine (BIBM), 2010 IEEE International Conference on 95–100 (IEEE, 2011).

[7]: Guindon, S. et al. New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0. Syst. Biol. 59, 307–321 (2010).

[8]: Stamatakis, A. RAxML version 8: A tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics 30, 1312–1313 (2014).

[9]:  Mirarab, S. et al. PASTA: Ultra-Large Multiple Sequence Alignment for Nucleotide and Amino-Acid Sequences. J. Comput. Biol. 22, 377–386 (2015).
