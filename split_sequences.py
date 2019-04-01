'''
Created on May 13, 2014
@author: namphuon
'''
import argparse
import os
from sepp.alignment import MutableAlignment


def parse_args():
    parser = argparse.ArgumentParser(
        description='Performs various tools for TIPP.')
    parser.add_argument(
        '-t', '--threshold', default=None, metavar='THRESHOLD',
        help='Split based on this threshold of length',
        type=int, dest='threshold')
    parser.add_argument(
        '-i', '--input', default=None, metavar='INPUT',
        help='INPUT sequence file (default=None)', required=True,
        type=str, dest='input')
    parser.add_argument(
        '-o', '--output', default='output', metavar='OUTPUT',
        help=('OUTPUT prefix, will write fragmentary file to OUTPUT.frag.fas '
              'and full-length file to OUTPUT.full.fas (default=output)'),
        type=str, dest='output')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    sequences = MutableAlignment()
    assert os.path.isfile(args.input) and os.access(args.input, os.R_OK), \
        "Input file %s does not exist\n" % args.input
    sequences.read_file_object(args.input)
    frag = MutableAlignment()
    full = MutableAlignment()

    for (key, seq) in sequences.items():
        if (len(seq) <= args.threshold):
            frag[key] = seq
        else:
            full[key] = seq
    frag.write_to_path("%s.frag.fas" % args.output)
    full.write_to_path("%s.full.fas" % args.output)


if __name__ == "__main__":
    main()
