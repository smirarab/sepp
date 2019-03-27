To make the sepp-package, you need to
make the following structures here:

```
ref/reference-gg-raxml-bl.tre
ref/reference-gg-raxml-bl-rooted.tre
ref/reference-gg-raxml-bl-rooted-relabelled.tre
ref/RAxML_info-reference-gg-raxml-bl.info
ref/gg_13_5_ssu_align_99_pfiltered.fasta
```
You have to get these from the [this repo](https://raw.github.com/smirarab/sepp-refs/master/gg/sepp-package.tar.bz) (too big for here)

You also need:
```
sepp/dendropy: [a copy or a symlink to the latest dendropy >4.0]
```
You can get this for example by running:
```
ln -s `echo 'import dendropy; import os; print (os.path.split(dendropy.__file__)[0])'|python` sepp/
```

Then archive them all into a single file:

```
cd ..
COPYFILE_DISABLE=1 tar cvfjh sepp-package.tar.bz  --exclude "*pyc" sepp-package/
```
