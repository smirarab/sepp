
Introduction 
===

UPP2 is a wrapper around UPP, and is subsequently also a method for performing alignments of ultra-large and fragmentary datasets. Please find more details about [UPP on its respective repository](https://github.com/smirarab/sepp/blob/master/README.UPP.md). We have two key strategies to note: 

* UPP-hs (Hierarchical) 
* UPP (Bitscore Weighting)

For a particular query sequence and a particular HMM, we've made some changes to allow UPP to recalculate the bitscore by which it finds the "best" HMM match for a particular query sequence. This recalculation now takes into account the number of sequences summarized in the HMM. We call this the adjustedUPP strategy. 

According to this new weighting the, we also present several strategies to redesign how the ensembles of HMMs are organized, taking into account the phylogenetic information. That is, we explored several strategies that allow us to incorporate more sequences or fewer sequences, and generate overlapping HMMs that each query sequence can be tested against. 

We have also implemented a fastUPP strategy which does a fast search down the hierarchy of trees. That is, recursively from the root, we choose the child HMM which best fits our current query sequence, continuing until the leaf. 

### Inputs

The main pipeline for UPP2 can be found in the script `run_upp2_v2.py`. The expected inputs for this are as follows: 
- A sequence file of unaligned sequences 
- Aligned backbone file
- True alignment file
- Pasta backbone tree file 
- Out tag
- HMMER directory 
- SEPP bundled package directory 
- decomposition parameter
- strategies file

The command to run UPP2 therefore follows the format: 
`python run_upp2_v2.py <unaligned_file> <aligned_backbone_file> <true_aligned_file> <pasta_backbone_tree_file> <out_tag> <hmmer_dir> <sepp_bundled_package_dir> <decomp> <doResort> <strategies_file> > <outputfile>`

For a test dataset, we have provided a fragmentary dataset under the `trial/` directory. Please run the following command there. 

```
python run_upp2_v2.py trial/1000M4/R0/unaligned_frag.txt trial/1000M4/R0/pasta_backbone_align.txt trial/1000M4/R0/true_align_fragged.txt trial/1000M4/R0/pasta_backbone.tre testing "" /Users/gillianchu/.sepp/bundled-v4.3.17/ 30 False trial/1000M4/R0/strategies.txt > testing.txt
```

### Outputs
This generates a folder called `alignData`. Inside `alignData`, we have a hierarchy of folders that will be generated. Most importantly, the `alignData/tmpfiles/` folder holds the HMMs constructed during the original UPP process. It will also hold the decomposition config file generated from the user-given parameters. The `alignData/UPPOutput/` folder will hold the alignments produced by UPP. 

The main output of UPP are the new predicted alignments according to each strategy. You can find this in the folder `alignData/hmmQueryList/merged/*_strategyName_alignmentFasta.fasta`. 

We generate the SPFN and SPFP of each predicted alignment by using FastSP, which is also bundled in the repository.

### Configurations


### Running UPP2


### Notes / To-Do's: 

UPP will be run with no parallelization. A later release may change this. 

* Why is hmmSearcher set to import from stefanHMM_concurrent throwing an error? 
* does stefan_UPP have to be run first? if so, why? 

* add a check for if we just ran saveInitialSteps for the right file, if so then can skip
* make sure file only has one output.xdlk;fd file in it

* make edits for when there is no true alignment file
* make some arguments required and others optional
* if the true alignments is specified in the strategies, then the true alignment file better be provided, etc. 
* add in save_scores.py, edit scoreAlignment to return the SP-score etc., and generate a xlsx report at the end comparing everything.
* take stefan's name out of the strategies...

* add option to run upp the way upp does i.e. should understand all the same flags and commands that upp does. 

* refactor the stefanHMM_concurrent script and rename the files.

* take out reference aln parameter, move save_scores.py into the score evaluation, take that flag out. 

* take out the strategies

* take out backbone and fragmentary sequences as input, move to running UPP as a module 

* change it so that UPP-hs isn't running UPP all the way.
