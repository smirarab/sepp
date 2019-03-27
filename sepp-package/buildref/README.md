In the following, we assume you want to build an ultra-large reference tree based on a new ultra-large alignment and tree.

We use as example, `gg_13_5_pasta_99_10masked.fasta` for the alignment and `99_otus.tree` for the tree.

1. First, you need to create a nice clean version of the reference tree:
   ```
   nw_topology -bI 99_otus.tree > 99_otus_nice.tree
   ```

2. You need to run the following to resolve polytomies:  
   ```
   t=99_otus_nice.tree
   alg=gg_13_5_pasta_99_10masked.fasta
   raxmlHPC-PTHREADS -s $alg -m GTRCAT -n score-$alg -g $t -T 16 -p $RANDOM
   ```
   note that here, `alg` is the alignment and `t` is the reference tree

3. You need to run the following to re-estimate branch lengths:
   ```
   raxmlHPC-PTHREADS -s $alg -m GTRCAT -n score-bl-$alg -F -f e -t RAxML_result.score-$alg -T 16 -p $RANDOM 
   ```

4. Now, you need to reroot your reference tree
   ```
   root.sh  archaea.99.ids RAxML_result.score-bl-$alg  reference-$alg-rooted.tre
   ```

5. You need to reformat your info file to the right format
   ```
   python reformat-info.py RAxML_info.score-bl-$alg > RAxML_info.$alg
   ```

6. Optional: you can add taxonomic labels back to your reference tree if your original tree (`$t`) had those
   ```
   python relabel.py 99_otus.tree reference-$alg-rooted.tre reference-$alg-rooted-relabelled.tre
   ```
   You can perhaps also just use tax2tree directly on the final output

7. [May be needed] Depending on the steps used, sometimes, the bacteria branch is left with no length. This can be solved using
   ``` bash
   sed -i -e "s/'k__Bacteria'/'k__Bacteria':0.0/g" reference-$alg-rooted-relabelled.tre
   ```

At the end, your reference files are: 
* your original alignment, 
* your final tree (`reference-$alg-rooted.tre` or `reference-$alg-rooted-relabelled.tre` if you included the last step)
* your final info file (`RAxML_info.$alg`)
