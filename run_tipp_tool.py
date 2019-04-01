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
        '-g', '--gene', default=None, metavar='GENE',
        help='use GENE\'s reference package',
        type=str, dest='gene')
    parser.add_argument(
        '-a', '--action', default=None, metavar='ACTION', help='Run ACTION',
        required=True, type=str, dest='action')
    parser.add_argument(
        '-o', '--output', default='output', metavar='OUTPUT',
        help='OUTPUT directory', type=str, dest='output')
    parser.add_argument(
        '-p', '--prefix', default='prefix', metavar='PREFIX', help='PREFIX',
        type=str, dest='prefix')
    parser.add_argument(
        '-i', '--input', default='input', metavar='INPUT',
        help='INPUT destination',
        type=str, dest='input')
    parser.add_argument(
        '-t', '--threshold', default=0.95, metavar='THRESHOLD',
        help='threshold for classification',
        type=float, dest='threshold')

    args = parser.parse_args()
    return args


root_p = open(os.path.join(os.path.split(
    os.path.split(sepp.__file__)[0])[0], "home.path")).readlines()[0].strip()
tipp_config_path = os.path.join(root_p, "tipp.config")


def profile(input, gene, output, prefix, threshold):
    sepp.config.set_main_config_path(tipp_config_path)
    opts = Namespace()
    sepp.config._read_config_file(open(tipp_config_path, 'r'), opts)
    (taxon_map, level_map, key_map) = sepp.metagenomics.load_taxonomy(
        "%s/refpkg/%s.refpkg/all_taxon.taxonomy" % (opts.reference.path, gene))
    gene_classification = sepp.metagenomics.generate_classification(
        input, threshold)
    sepp.metagenomics.remove_unclassified_level(gene_classification)
    sepp.metagenomics.write_classification(
        gene_classification, "%s/%s.classification" % (output, prefix))
    sepp.metagenomics.write_abundance(gene_classification, output)


def main():
    args = parse_args()
    if (args.action == 'profile'):
        profile(args.input, args.gene, args.output, args.prefix,
                args.threshold)


if __name__ == "__main__":
    main()
#
