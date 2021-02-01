#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
Created on May 13, 2014
@author: namphuon
'''
import argparse
import os
from argparse import Namespace
import sepp.config
import sepp.metagenomics
import sepp


def parse_args():
    parser = argparse.ArgumentParser(
        description='Performs various tools for TIPP.')
    parser.add_argument(
        "-g", "--gene",
        dest="gene", metavar="GENE",
        default=None,
        type=str,
        help="use GENE\'s reference package")
    parser.add_argument(
        "-o", "--output",
        dest="output", metavar="OUTPUT",
        default="output",
        type=sepp.config.valid_file_prefix,
        help="output files with prefix OUTPUT. [default: %(default)s]")
    parser.add_argument(
        "-d", "--outdir",
        dest="outdir", metavar="OUTPUT_DIR",
        default=os.path.curdir,
        type=sepp.config.valid_dir_path,
        help="output to OUTPUT_DIR directory. full-path required. "
             "[default: %(default)s]")
    parser.add_argument(
        "-i", "--input",
        dest="input", metavar="INPUT",
        default=None,
        type=str,
        help="input file. full-path required. "
             "(_classification.txt file from running run_tipp.py)")
    parser.add_argument(
        "-t", "--threshold",
        dest="threshold",  metavar="THRESHOLD",
        default=0.95,
        type=float,
        help="Threshold for classification [default: 0.95]")
    return parser.parse_args()


def profile(inputf, gene, output, prefix, threshold):
    root_p = open(os.path.join(os.path.split(
                  os.path.split(sepp.__file__)[0])[0], "home.path"))
    root_p = root_p.readlines()[0].strip()
    tipp_config_path = os.path.join(root_p, "tipp.config")
    sepp.config.set_main_config_path(tipp_config_path)
    opts = Namespace()
    sepp.config._read_config_file(open(tipp_config_path, 'r'), opts)
    (taxon_map, level_map, key_map) = sepp.metagenomics.load_taxonomy(
        "%s/%s.refpkg/all_taxon.taxonomy" % (opts.reference.path, gene))
    gene_classification = sepp.metagenomics.generate_classification(
        inputf, threshold)
    sepp.metagenomics.remove_unclassified_level(gene_classification)
    tstr = str("%d" % (threshold * 100))
    cfile = str("%s/%s.classification_%s.txt" % (output, prefix, tstr))
    sepp.metagenomics.write_classification(gene_classification, cfile)
    sepp.metagenomics.write_abundance(gene_classification, output)


def main():
    args = parse_args()
    profile(args.input, args.gene, args.outdir, args.output,
            args.threshold)


if __name__ == "__main__":
    main()
