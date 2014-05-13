'''
Created on May 13, 2014
@author: namphuon
'''
import sys,argparse
from argparse import ArgumentParser, Namespace
from sepp.alignment import MutableAlignment, ExtendedAlignment,_write_fasta

def parse_args():
  parser = argparse.ArgumentParser(description='Separates sequences out based upon length threshold.')
  parser.add_argument('-t', '--threshold', help='threshold for fragmentary sequences, inclusive',   
                   type=int,dest='threshold')
  parser.add_argument('-i', '--input', help='input sequence file',   
                   type=str,dest='input')
  parser.add_argument('-o', '--output', help='output prefix, will write fragmentary file to OUTPUT.frag.fas and full-length file to OUTPUT.full.fas',   
                   type=str,dest='output')
  args = parser.parse_args()                   
  return args

def main():
  args = parse_args()
  sequences = MutableAlignment()
  sequences.read_file_object(args.input)
  frag = MutableAlignment()
  full = MutableAlignment()
  
  for (key,seq) in sequences.items():
    if (len(seq) <= args.threshold):
      frag[key]=seq
    else:
      full[key]=seq
  frag.write_to_path("%s.frag.fas" % args.output)
  full.write_to_path("%s.full.fas" % args.output)


if __name__ == "__main__":
    main()
#


