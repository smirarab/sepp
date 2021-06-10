import os, glob
import subprocess
import argparse

from sepp.hmm_concurrent import *
from sepp.hmm_searcher import *
from sepp.scheduler import JobPool

def create_dirs(path_name):
    if not os.path.exists(path_name): 
        subprocess.call(['mkdir', path_name])

def parse_strats(stratsfile): 
    stratlist = []
    with open(stratsfile, "r") as reader:
        lines = reader.readlines()
        for line in lines:
            stratlist.append(line.strip())
    return stratlist


def makedirstruct(dirpath): 
    print("[makedirstruct]")
    # clear the folders
    dirpath = os.path.join(dirpath, "ensembleData")
    create_dirs(dirpath)
    for firstlvl in ["trueAlignment","temporaryFileSave","subsetTrueAln", "Searcher","queryToHmm","newHMM","ML","initialHMM","hmmSeqAlign","hmmScores","hmmQueryList","fullPrediction", "seqFileNames"]:
        
        firstdir = '%s/%s' % (dirpath, firstlvl)
        if os.path.isdir(firstdir):
            subprocess.call(['rm', '-rf', firstdir])

        subprocess.call(['mkdir', firstdir])
        if firstlvl == 'fullPrediction':
            subprocess.call(['mkdir', firstdir + '/sequences/'])
        elif firstlvl == 'hmmQueryList':
            for z in ['inputQuery', 'merged', 'predictedQuery']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'hmmScores':
            for z in ['newRawOld', 'processedOld', 'rawOld', 'scoresFull', 'temporaryStorage']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'ML':
            for z in ['hmm', 'models', 'queryHMM', 'scores', 'sequences', 'temporaryStorage']:
                subprocess.call(['mkdir', firstdir + '/' + z])
            subprocess.call(['mkdir', firstdir + '/sequences/np'])
        elif firstlvl == 'newHMM':
            for z in ['columnSets', 'hmm', 'newHMMseq']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'queryToHmm':
            for z in ['original', 'withUPP']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'Searcher':
            subprocess.call(['mkdir', '%s/scoreFiles' % firstdir])
        elif firstlvl == 'trueAlignment':
            for z in ['original', 'subset']:
                subprocess.call(['mkdir', firstdir + '/' + z])


def run_upp_strats(abstract_algorithm, dirpath, decomp, strats, out_tag, trueAlign, doResort=False):
    ''' Note: abstract_algorithm can be used to do commands like below.
        The moment you call addHMMBuildJob, jobs will get enqueued and run eventually
    addHMMBuildJob(abstract_algorithm, <hmmbuild profile output file path>, <fasta file path>)
    JobPool().wait_for_all_jobs()
    '''
    return
    print("[run_upp_strats]")
    decomp = str(decomp)

    print("[START HERE] decomp:%s" % (decomp), flush=True)
    dirname = "%s/tmpfiles/%s/" % (dirpath, out_tag)
    outputdirname = os.listdir(dirname)[0]
    hmmSeqFile = '%s/%s/root/P_0/' % (dirname, outputdirname)
    fragmentfile = os.listdir(dirname + outputdirname + "/fragment_chunks/")[0]
    queryName = "%s/%s/fragment_chunks/%s" % (dirname, outputdirname, fragmentfile)
    trueAlignment = trueAlign
    predictionName = './%s/UPPoutput/%s_output_alignment.fasta' % (dirpath, out_tag)

    print("hmmseqfile: %s" % hmmSeqFile)
    print("queryName: %s" % queryName)
    print("trueAlignment: %s" % trueAlignment)
    print("predictionName: %s" % predictionName)

    setAllFileNames(hmmSeqFile, queryName, trueAlignment, predictionName)
    saveInitialSteps()

    hierchySearch()

    for strat in strats: 
        if strat != 'stefan_trueUPP':
            print("[processing %s]" % strat)
            print("[running scoresToHMMSeq]")
            scoresToHMMSeq(strat)
            print("[running buildAlignMerge, doResort is %s]" % doResort)
            buildAlignMerge(strat, doResort=doResort)
        #print("[running scoreAlignment]")
        #scoreAlignment(strat)

# def main(): 
#     parser = argparse.ArgumentParser()

#     ## take in upp tree
#     ## take in upp backbone alignment
    
#     parser.add_argument('unaligned_file', type=str, help="relative path for unaligned file")
#     parser.add_argument('aligned_file', type=str, help="relative path for aligned file")
#     parser.add_argument('true_align', type=str, help='relative path for true align file')
#     parser.add_argument('treefile', type=str, help="relative path for tree file")
#     parser.add_argument('out_tag', type=str, help="out tag")

#     parser.add_argument('hmmer_packagedir', type=str, help="hmmer package dir")
#     parser.add_argument('bundle_packagedir', type=str, help="bundled package dir")

#     parser.add_argument('decomp', type=int, help="decomp minimumsubset size, default is 10")
#     parser.add_argument('doResort', type=str, help="True/False, default is False")
#     parser.add_argument('strats', type=str, help="file with strategies to run one per line, see README for more details")
#     args = parser.parse_args()

#     unaligned_file = args.unaligned_file
#     align_file = args.aligned_file
#     treefile = args.treefile
#     trueAlign = args.true_align
#     out_tag = args.out_tag

#     bundle_packagedir = args.bundle_packagedir
#     packagedir = args.hmmer_packagedir

#     decomp = args.decomp
#     doResort = args.doResort
#     strats = parse_strats(args.strats)

#     dirpath = 'alignData' #'ensemble_data'
#     create_dirs(dirpath)

#     root_tmpdir = '%s/tmpfiles/' % dirpath
#     create_dirs(root_tmpdir)
#     outdir = '%s/UPPOutput/' % dirpath 
#     create_dirs(outdir)

#     build_upp_config(root_tmpdir, unaligned_file, bundle_packagedir, packagedir, align_file, treefile, outdir, out_tag, decomp)
    
#     #makedirstruct(dirpath)
#     #run_upp_strats(dirpath, decomp, strats, out_tag, trueAlign, doResort)

    
# if __name__ == '__main__':
#     main()
