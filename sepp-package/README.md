This pacakge enables you to run SEPP on the [Greengenes](http://greengenes.lbl.gov/cgi-bin/nph-index.cgi) reference tree with 203,452 leaves.

### Download/install

First, download [this package](https://raw.github.com/smirarab/sepp-refs/master/gg/sepp-package.tar.bz) and decompress it; 
```
wget  "https://raw.github.com/smirarab/sepp-refs/master/gg/sepp-package.tar.bz"
tar xvfj sepp-package.tar.bz
```

then run:

``` bash
cd sepp-package/sepp
python setup.py config -c
```
to configure the package. 

Anytime you move the package to a new directory, you need to run this again. 

**Note:**  this is a self-contained package and does not require usual SEPP installation steps. 
In fact, it will not use your installed SEPP even if you have it installed. 

### Running
Once configured, you can run

``` bash
./sepp-package/run-sepp.sh [input file] [prefix of the output]
```

This will add the sequences from the input file to the Greengene tree. 

A test file is provided. You can test by running:

```
./sepp-package/run-sepp.sh sepp-package/test.frag test-gg
```

### Results

* You will get a `.json` file, a `.xml` file, and a `.tre` file. 
    * The `.json` file includes full results of placement. The `.json` file can be used with the [guppy](https://matsen.github.io/pplacer/generated_rst/guppy.html) package from pplacer to create visualizations and other stuff. To learn more about these, see [SEPP's tutorial](https://github.com/smirarab/sepp/tree/master/tutorial).
    *  The `.tre` file is in the newick format and includes the best placement for each fragment. 
    * The xml file is similar, but in the NexML format. * The reference tree has taxon ids instead of names. The mapping between names and the taxonomy is given in the file `sepp-package/ref/99_otu_taxonomy.txt`. 
* You can use  [tax2tree](https://github.com/biocore/tax2tree) to decorate an existing taxonomy onto the resulting phylogenetic tree. However, with the current
version, the final newick output is already decorated. 

### Notes on installation/running

* Any version of python works
* If you already have SEPP installed and you want to simply use your current installation, instead of the package, don't despair. 
  Check out the [DEVGUIDE.md](DEVGUIDE.md) and you will get an idea of how this can be done.
  You simply need to simlink `ref` and `python` to the right place. 
* To pass any argument to SEPP, you can pass them to `run-sepp.sh`; all SEPP arguments (except the input fragment file, output prefix, and reference alignment/tree) can be provided to `run-sepp.sh`.
* By default, this package will use all available cores on your node. Pass in `-x 4` to manually control the number of threads to limit it to 4. 
  This may prove important if you have too many cores and not enough memory per core. 
* Be default, subset size and placement size are set to 1,000 and 5,000. To change, pass in `-A` and `-P` options to `run-sepp.sh`.
* **Temporaries:** (ignore tips below if they don't make sense to you)
    * To facilitate portability, we are creating the temporary directories using `mktemp`. 
      On your machine, this may not be the best strategy. 
      To improve performance, you want to use the best available partition you have for temporary directories on your machine. 
      For example, you may have an SSD partition. 
      To change the location of temporaries, set environmental variables TMPDIR and  TMPDIRSSD. See the file for examples. 
    * After you are done with your run, you may want to manually remove the temporaries or you may want to copy them to a permanent place. 
      There are two temporary directories: the main one (`sepp-temp-XXXX`) will have the alignment files, which you may want to keep (these are huge). 

### Notes on the reference alignment/tree

* To build the reference tree, see [this page](buildref/README.md).
    * We started from [this file](ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/trees/99_otus.tree) from the public Greengenes dataset. 
    * Taxonomic labels are removed (for the next few steps)
    * The reference alignment is called `gg_13_5_ssu_align_99_pfiltered.fasta` and is included within the package
    * The reference tree is resolved and branch lengths are re-estimated using RAxML:
    ```
    raxmlHPC-PTHREADS -s gg_13_5_ssu_align_99_pfiltered.fasta -m GTRCAT -n score-f -F -g 99_otus_nice.tree -T 16 -p 32323
    ```
    the use of `GTRCAT` and `-F` are not ideal, but are necessary because of the dataset size
    * The tree is rerooted and taxonomic branch labels are put back on the tree using [this script](buildref/relabel.py)


#### Acknowledgments 

* We thank Daniel MCDonald for help with this package and contributed code. 
