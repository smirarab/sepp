Introduction
===

UPP stands for `Ultra-large alignments using Phylogeny-aware Profiles`, and so is a method for the following problem:

Alignment:
- Input: A set of unaligned sequences `S`
- Output: An alignment `A` on `S`

UPP is a modification of SEPP for performing alignments of ultra-large and fragmentary datasets.  UPP operates in four steps.  In the first step, UPP partitions set `S` into a backbone set and a query set and computes an alignment and tree on the backbone set using [PASTA](https://github.com/smirarab/pasta) (Mirarab et al., RECOMB 2014 and Journal of Computational Biology 2014), which is a direct improvement to SATe (Liu et al., Science 2009 and Systematic Biology 2012).  In the next step, UPP decomposes the backbone alignment into an ensemble of profile Hidden Markov Models (HMMs).  The third step in UPP searches for the best alignment of the query sequence to each HMM.  The final step inserts the query sequence into the backbone alignment using the best scoring HMM.  Our study shows that UPP results in accurate alignments, and that ML trees estimated on the alignments are also highly accurate. UPP has good accuracy on datasets that contain fragmentary sequences. 

UPP(Default): The default version selects 1000 sequences at random for the backbone alignment. If the dataset has at most 1000 sequences, this means that UPP(Default) is identical to PASTA. 

UPP(Fast): We have designed a fast version of UPP that uses a backbone with at most 100 sequences. The default version uses a backbone of 1000 sequences. The fast version can produce an alignment on 10,000 sequences in less than an hour using 12 processors, and on 1,000,000 sequences in less than 12 days; the default version requires seven hours on 10,000 sequences and would take an estimated 120 days on 1,000,000 sequences. The difference in accuracy between UPP(Fast) and UPP(Default) depends on the rate of evolution -- datasets with low to moderate evolutionary diameters can be analyzed well with UPP(Fast); otherwise, we recommend the use of UPP(Default). However, on large datasets, UPP(Default) will take nearly ten times as much running time.

FRAGMENTARY DATASETS: UPP can be used in default mode, which will select the backbone sequences randomly and without trying to restrict the backbone to full length sequences. However, if the dataset contains fragments, then UPP should be used in a mode that restricts the backbone to just the "full-length" sequences. To do this, you will need to provide UPP with an estimate of the full length of sequences for your locus. See Advanced Usage information about how to do this.

SUPPLYING YOUR OWN SEED ALIGNMENT AND TREE: If you have a pre-computed seed alignment (with phylogenetic tree), you can provide this to UPP. See Advanced Usage information below about how to do this.

PARALLEL IMPLEMENTATION: UPP is embarrassingly parallel. See Advanced Usage information about how to do this.


Developers: Nam Nguyen, Siavash Mirarab, and Tandy Warnow.

###Publication:
Nam Nguyen, Siavash Mirarab, Keerthana Kumar, and Tandy Warnow. `Ultra-large alignments using Phylogeny Aware Profiles`. Accepted to RECOMB 2015 (Research in Computational Molecular Biology 2015) and Genome Biology.

---
Installation
===

You have two options for installing UPP. 

 - **Windows:** If you have a Windows machine, currently using the Virtual Machine (VM) image we provide is your only option. 
 - **Linux:** and **MAC:**  If you have Linux (or other \*nix systems) or MAC, you can still use VM, but downloading the code from github and installing it is what we strongly recommend. 
 
 ### 1. From Source Code
Current version of UPP has been developed and tested entirely on Linux and MAC. 
Windows won't work currently (future versions may or may not support Windows). 

You need to have:

- Python (version 2.7 or later)
- Java (version > 1.5)
- PASTA (version 1.0 or later)

**Installation of SEPP**:
UPP is a part of the SEPP distribution package.  First download and install SEPP:

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

**Installation of UPP**:

Next install [PASTA] (https://github.com/smirarab/PASTA/) and make sure the run_pasta.py executable is on the PATH variable.  Once done, do the following. 

1. Run the following command:

```
python setup.py upp
```

The last step creates a ~/.sepp/upp.config config file. Since this is specific to a user, each user that runs UPP needs to execute the last step. 

**Common Problems:**
1.  UPP requires SEPP to be installed.  If UPP is not running, first check to see if UPP was installed correctly.

2.  UPP requires PASTA to be installed and the run_pasta.py executable to be on the path.  

### 2. From Virtual Machine (VM)

VM Image (mostly for Windows users) is available for [download](http://www.cs.utexas.edu/~phylo/software/PASTA_TIPP_UPP.ova) (2.5 GB). Once the image is downloaded, you need to run it using a VM environment ([VirtualBox](https://www.virtualbox.org/) is a good option). After you install VirtualBox, you just need to use File/import to import the PASTA_TIPP_UPP.ova image that you have downloaded (If your machine has less than 3GB you might want to reduce the memory to something like 512MB). Once VM is imported, you can start it from the Virtualbox. If you are asked to login, the username and passwords are (username: phylolab, password: phylolab). TIPP and UPP are already installed on the VM machine, so you can simply proceed by opening a terminal and running it.

Email: `ensemble-of-hmms@googlegroups.com` for all issues. 

