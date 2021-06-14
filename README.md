
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
- decomposition parameter (optional)
- optional flag to run original UPP all the way through (or not)
- optional flag to run with adjusted bitscore (bitscore_adjust)
- optional flag to run with hierchical search (upp_hier)

The command to run UPP2 therefore follows the format: 
`python run_upp.py <unaligned_file> <hmmer_dir> <decomp> <config file> <tmpdir> <other optional flags> > <outputfile>`

Note that the aligned backbone and corresponding tree are optional inputs to UPP. For a test dataset, we have provided a fragmentary dataset under the `trial/` directory. Please run the following command there. 

```
python run_upp.py -s trial/1000M4/R0/unaligned_frag.txt -x 1 -a trial/1000M4/R0/pasta_backbone_align.txt -t trial/1000M4/R0/pasta_backbone.tre -A 30 -c trial/1000M4/R0/decomp30.config --tempdir tmpfiles2 -j True
```

### Outputs
This generates a folder called `alignData`. Inside `alignData`, we have a hierarchy of folders that will be generated. Most importantly, the `alignData/tmpfiles/` folder holds the HMMs constructed during the original UPP process. It will also hold the decomposition config file generated from the user-given parameters. The `alignData/UPPOutput/` folder will hold the alignments produced by UPP. 

The main output of UPP are the new predicted alignments according to each strategy. You can find this in the folder `tmp/output.XX/data/internalData/default-value-not-emptyp/<strategyName>/hmmQueryList/merged`. 


### Notes / To-Do's: 

UPP will be run with no parallelization. A later release may change this. 

* add a check for if we just ran saveInitialSteps for the right file, if so then can skip
* add in save_scores.py, edit scoreAlignment to return the SP-score etc., and generate a xlsx report at the end comparing everything.
* take stefan's name out of the strategies
* *split the hmm_concurrent script*
* clean up directory structure inside /tmp
* Make the ./data folder part of the /tmpfiles/ensembleData folder... integrate
* Figure out where the resulting alignments are going... folder hierarchy is being messy... 
