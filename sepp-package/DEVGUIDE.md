To make the sepp-package, you need to
make the following structures here:

```
ref/reference-gg-raxml-bl.tre
ref/reference-gg-raxml-bl-rooted.tre
ref/reference-gg-raxml-bl-rooted-relabelled.tre
ref/RAxML_info-reference-gg-raxml-bl.info
ref/gg_13_5_ssu_align_99_pfiltered.fasta
```
You also need:
```
sepp/dendropy: [a copy or a symlink to the latest dendropy >4.0]
```

Then zip them all up

```
cd ..
rm sepp-package.zip; zip -r sepp-package.zip sepp-package/ --exclude \*.pyc
```
