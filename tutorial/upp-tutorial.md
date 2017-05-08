Introduction
===

UPP stands for `Ultra-large alignments using Phylogeny-aware Profiles`, and so is a method for the following problem:

Alignment:
- Input: A set of unaligned sequences `S`
- Output: An alignment `A` on `S`

UPP is a modification of SEPP for performing alignments of ultra-large and fragmentary datasets.  UPP operates in four steps.  In the first step, UPP partitions set `S` into a backbone set and a query set and computes an alignment and tree on the backbone set using [PASTA](https://github.com/smirarab/pasta) (Mirarab et al., RECOMB 2014 and Journal of Computational Biology 2014), which is a direct improvement to SATe (Liu et al., Science 2009 and Systematic Biology 2012).  In the next step, UPP decomposes the backbone alignment into an ensemble of profile Hidden Markov Models (HMMs).  The third step in UPP searches for the best alignment of the query sequence to each HMM.  The final step inserts the query sequence into the backbone alignment using the best scoring HMM.  Our study shows that UPP results in accurate alignments, and that ML trees estimated on the alignments are also highly accurate. UPP has good accuracy on datasets that contain fragmentary sequences. 

* UPP(Default): The default version selects 1000 sequences at random for the backbone alignment. If the dataset has at most 1000 sequences, this means that UPP(Default) is identical to PASTA. 

* UPP(Fast): We have designed a fast version of UPP that uses a backbone with at most 100 sequences. The default version uses a backbone of 1000 sequences. The fast version can produce an alignment on 10,000 sequences in less than an hour using 12 processors, and on 1,000,000 sequences in less than 12 days; the default version requires seven hours on 10,000 sequences and would take an estimated 120 days on 1,000,000 sequences. 

* The difference in accuracy between UPP(Fast) and UPP(Default) depends on the rate of evolution. Datasets with low to moderate evolutionary diameters can be analyzed well with UPP(Fast); otherwise, we recommend the use of UPP(Default). However, on large datasets, UPP(Default) will take nearly ten times as much running time.

**Custom seed alignment and tree**: If you have a pre-computed seed alignment (with phylogenetic tree), you can provide this to UPP. See Advanced Usage information below about how to do this.

**Parallel implementation**: UPP, like SEPP and TIPP, is embarrassingly parallel. See Advanced Usage information about how to do this.


#### Fragmentary Datasets

UPP can be used in default mode, which will select the backbone sequences randomly and without trying to restrict the backbone to full length sequences. However, if the dataset contains fragments, then UPP should be used in a mode that restricts the backbone to just the "full-length" sequences. To do this, you will need to provide UPP with an estimate of the full length of sequences for your locus. See Advanced Usage information about how to do this.
Note that if there are no fragmentary sequences, running UPP using this option will be equivalent of running PASTA on the entire dataset.



Developers: Nam Nguyen, Siavash Mirarab, and Tandy Warnow.

####Publication


* Nguyen, Nam-phuong D., Siavash Mirarab, Keerthana Kumar, and Tandy Warnow. 2015. “Ultra-Large Alignments Using Phylogeny-Aware Profiles.” Genome Biology 16 (1): 124. doi:[10.1186/s13059-015-0688-z](http://genomebiology.com/2015/16/1/124).

### Note and Acknowledgment: 
- UPP bundles the following two programs into its distribution:
  2. hmmer: http://hmmer.janelia.org/
- UPP uses the [Dendropy](http://pythonhosted.org/DendroPy/) package. 
- UPP uses some code from [SATe](http://phylo.bio.ku.edu/software/sate/sate.html).


Installation
===

You have two options for installing UPP. 

 - **Windows:** If you have a Windows machine, currently using the Virtual Machine (VM) image we provide is your only option. 
 - **Linux:** and **MAC:**  If you have Linux (or other \*nix systems) or MAC, you can still use VM, but downloading the code from github and installing it is what we strongly recommend. 
 
###  1. From Source Code
Current version of UPP has been developed and tested entirely on Linux and MAC. 
Windows won't work currently (future versions may or may not support Windows). 

You need to have:

- Python (version 2.7 or later)
- Java (version > 1.5)
- PASTA (version 1.0 or later)


#### Step 1: Install SEPP
UPP is a part of the SEPP distribution package.  First download and install SEPP:

1. Open a terminal and create a directory where you want to keep SEPP. e.g. `mkdir ~/sepp-code`. Go to this directory. e.g. `cd ~/sepp-code`.

2. Clone the SEPP code repository from our [github repository](https://github.com/smirarab/sepp). For example you can use `git clone https://github.com/smirarab/sepp.git`.
If you don't have git, you can directly download a [zip file from the repository](https://github.com/smirarab/sepp/archive/master.zip)
and decompress it into your desired directory. 

3. `cd sepp` (or `cd sepp-master` if you used the zip file instead of cloning the git repository)

4. Run the following command from the SEPP directory:
    ```
    python setup.py config
    ```
   or alternatively, 
    ```
    python setup.py config -c
    ```
   The first command puts all the config files under your home directory while the second version uses the current path. 
4. Then run: 
   ```
   sudo python setup.py install
   ``` 
   If you don't have root access, remove the `sudo` part and instead  use  `--user` option. Alternativley, you can `--prefix` to install in a different location, but that different location needs to be part of your `PYTHONPATH` environmental variable.


#### Step 2: Install UPP

1. install [PASTA](https://github.com/smirarab/PASTA/) and make sure the `run_pasta.py` executable is on the `PATH` variable.  

2. Run the following command from the SEPP directory:
   ```
   python setup.py upp
   ```
   or 
   ```
   python setup.py upp -c
   ```

**Note:** It's important tat you use either use `-c` for both SEPP and UPP or don't use it for both. 



**Common Problems:**
1.  UPP requires SEPP to be installed.  If UPP is not running, first check to see if UPP was installed correctly.

2.  UPP requires PASTA to be installed and the `run_pasta.py` executable to be on the path.  

3. SEPP is installed with `-c` but UPP is not.


### 2. From Virtual Machine (VM)

VM Image (mostly for Windows users) is available for [download](https://drive.google.com/open?id=0B0lcoFFOYQf8U2NZV2Z2RmRaRjQ) (3 GB).
Once the image is downloaded, you need to run it using a VM environment ([VirtualBox](https://www.virtualbox.org/) is a good option).
After you install VirtualBox, you just need to use File/import to import the phylo.ova image that you have downloaded (If your machine has less than 3GB you might want to reduce the memory to something like 512MB). Once VM is imported, you can start it from the Virtualbox. If you are asked to login, the username and passwords are (username: phylolab, password: phylolab). TIPP and UPP are already installed on the VM machine, so you can simply proceed by opening a terminal and running it.

Note that we constantly update our software.  Before running the tutorial, it's best to grab
the most updated version of the software onto the VM machine.  
This can be done by opening a terminal in the VM and typing the following commands:


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

---------
Using UPP
===

If your installation is successful, you should be able to run UPP by running the following command from any location. Open up a terminal window and type: 

```
run_upp.py -h
``` 

Running UPP with the `-h` option produces the list of options available in UPP. 

The general command for running UPP under the default settings is:

```
run_upp.py -A alignment_size -B backbone_size -M -1 -m molecule_type -s input
```

Step 1: Running a test job
---

UPP currently can only be run from the command line.  We have provided some test data files under the `test/` directory.  A good start is running a small amino acid dataset found in `test/unittest/data/upp_frag/amino.fas`.  

Let's take a look at this file:

```
>SEQ0
IVSNASCTTNCLAPLAKVINDNFGIIEGLMTTVHATTATQKTVDGPSHKDWRGGRGASQNIIPSSTGA
>SEQ1
SKIGINGFGRIGRLVLRTALEMGAQVVAVNDPFIALEYMVYMFKYDSTHGMFKGEVKVEDGALVVDGKKITVFNEMKPENIPWSKAGAEYIVESTGVFTTIEKASAHFKGGAKKVIISAPSADAPMFVCGVNLEKYSKDMKVVSNASCTTNCLAPVAKVLHENFEIVEGLMTTVHAVTATQKTVDGPSAKDWRGGRGAAQNIIPSSTGAAKAVGKVIPELDGKLTGMAFRVPTPNVSVVDLTVRLGKECSYDDIKAAMKTASEGPLQGVLGYTEDDVVSCDFTGDNRSSIFDAKAGIQLSKTFVKVVSWYDNEFGYSQRVIDLIKHMQKVDS
>SEQ10
EYMTIKVGINGFGRIGRIVFRAAQKRSDIEIVAINDLLDADYMAYMLKYDSTHGRFDGTVEVKDGHLIVNGKKIRVTAERDPANLKWDEVGVDVVAEATGLFLTDETARKHITAGAKKVVMTGPSKDNTPMFVKGANFDKYAGQDIVSNASCTTNCLAPLAKVINDNFGIIEGLMTTVHATTATQKTVDGPSHKDWRGGRGASQNIIPSSTGAAKAVGKVLPELNGKLTGMAFRVPTPNVSVVDLTVRLEKAATYEQIKAAVKAAAEGEMKGVLGYTEDDVVSTDFNGEVCTSVFDAKAGIALNDNFVKLVSWYDNETGYSNKVLDLIAHISK
>SEQ102
MTTVHAITATQKTVDGPSGKLWRDGRGAAQNIIPASTGAAKAVGKVIPELNGKLTGMAFRVPVHDVSVVDLTCRLSKEASY
....
```

You can see that this file contains both full-length and fragmentary sequences.  Let's run UPP on this dataset:

```
run_upp.py -A 10 -B 1000 -M -1 -m amino -s test/unittest/data/upp_frag/amino.fas
```

This command results in UPP building a backbone alignment and tree on any sequences that is between 75% to 125% the median length sequence; all other sequences are treated as fragmentary.  The backbone will be generated using PASTA.  The remaining sequences are then aligned using UPP.

#### Output files:
The main outputs of UPP are two alignmet files, `<prefix>_alignment.fasta` and `<prefix>_alignment_masked.fasta`.  The  `<prefix>_alignment.fasta` file is the alignment of the unaligned sequences.  The `<prefix>_alignment_masked.fasta` is the masked alignment file; non-homologous sites in the query set are removed.  Note that UPP does not report a tree.  Trees can be computed using your favorite tree estimation software.  We prefer RAxML or FastTree (and more recently, IQTree).

The secondary outputs are the backbone alignment and tree (always named as `pasta.fasta` and `pasta.fasttree`) and the list of insertion columns (named `<prefix>_insertion_columns.txt`). 



Step 2: Viewing the resulting alignment
---

### Alignment viewing software:
Many applications exist for viewing multiple sequence alignment. Some options :

1. Alignments can be viewed in any text editor. In many situations, this is sufficient. In Linux/Mac command-line, one can use ``less``, ``vim``, or any number of other text editor applications. Also GUI-enable text editors (e.g. ``TextEdit`` for MAC) would also work. Windows ``notepad`` is not a good option (neither is MS word or any other word processing tool), but `notepad++` should work. 
2. [Seaview](http://doua.prabi.fr/software/seaview) is a relatively light-weight application that does a good job of visualizing alignments
3. [JalView](http://www.jalview.org/download) has many options, but is not necessarily light-weight. 
4. [SuiteMSA](http://bioinfolab.unl.edu/~canderson/SuiteMSA/) is also a full-featured software for alignment viewing, manipulation, and more. 

In this tutorial, we will use both a text editor and the Seaview. But feel free to use your own favorite tools. We will open ``output_alignment.fasta`` using SeaView and will look at them. 


Step 3: Using an existing alignment and tree
---
Suppose that you have computed an accurate alignment and tree using your favorite methods.  You can use the alignment and tree to align the remaining sequences in your dataset.  Go to the `test/unittest/data/upp/` directory and type:

```
run_upp.py -A 10 -s query.fas -a test.fasta -t test.tree
```

This command will run UPP using the backbone alignment and tree given in the command line.

Step 4: Running UPP from a config file
---

You can also run UPP from a config file.  From the same directory, type:

```
run_upp.py -c sample.config -o config_example
```

---------
Contact
===

Post all questions, comments, requests to: <https://groups.google.com/forum/#!forum/ensemble-of-hmms>

