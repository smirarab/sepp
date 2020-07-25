#!/usr/bin/env python


import itertools
from sepp.alignment import ExtendedAlignment


gencode = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'  # ,
    #    'CCA':'Z', 'CCC':'Z', 'CCG':'Z', 'CCT':'Z',
    #    'GCA':'Z', 'GCC':'Z', 'GCG':'Z', 'GCT':'Z'
}

dnacode = {'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
           'S': ['G', 'C'], 'R': ['G', 'A'], 'Y': ['T', 'C'],
           'W': ['A', 'T'], 'M': ['A', 'C'], 'K': ['G', 'T'],
           'B': ['G', 'C', 'T'], 'H': ['A', 'C', 'T'],
           'D': ['G', 'A', 'T'], 'V': ['G', 'C', 'A'],
           'N': ['A', 'C', 'G', 'T']}


def is_compatible(cd, aa):
    if aa == 'Z':
        return is_compatible(cd, 'E') or is_compatible(cd, 'Q')
    elif aa == 'B':
        return is_compatible(cd, 'N') or is_compatible(cd, 'D')
    else:
        return aa == 'X' or aa in set(gencode[''.join(a)]
                                      for a in itertools.product(
            [''.join(a) for a in itertools.product(
                dnacode[cd[0]], dnacode[cd[1]])], dnacode[cd[2]]))


def is_ambiguous(cd):
    return len(set(gencode[''.join(a)]
                   for a in itertools.product(
        [''.join(a) for a in itertools.product(
            dnacode[cd[0]], dnacode[cd[1]])],
        dnacode[cd[2]]))) > 1


def backtranslate(faa, fna):
    newfna = ExtendedAlignment(faa.fragments)
    for k, s in fna.items():
        if k in faa.keys():
            aa = faa[k].upper()
            cd = []
            i = 0
            for r in aa:
                cds = s[i:i + 3]
                if r == '-':
                    cd.append('---')
                else:
                    if is_compatible(cds, r):
                        cd.append(cds)
                        i += 3
                    else:
                        if i == 0 and (cds == 'GTG' or cds == 'TTG'):
                            cd.append(cds)
                            i += 3
                        else:
                            raise ValueError(
                                '%s at position %d of %s '
                                'does not translate to %s' % (cds, i, k, r))
            newfna[k] = ''.join(cd)
        else:
            continue
    col_lab = faa.col_labels
    for i in col_lab:
        newfna._col_labels = newfna._col_labels + [i, i, i]
    return newfna
