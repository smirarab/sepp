#!/bin/env python

import sys

raxmlinfo = '''
Base frequencies: %f %f %f %f

Inference[0]: Time %f CAT-based likelihood -0000, best rearrangement setting 5
alpha[0]: 1.000000 rates[0] ac ag at cg ct gt: %f %f %f %f %f %f


NOT conducting any final model optimizations on all 1 trees under CAT-based
model ....

Final GAMMA  likelihood: %f
'''

ls = open(sys.argv[1], 'r').readlines()

for (i, l) in enumerate(ls):
    if "called as follows" in l:
        break
i += 3

v = dict()

for l in ls[i:-1]:
    lr = l[0:-1].split(' ')
    if "<->" in l:
        v[' '.join(lr[1:4])] = float(lr[4])
    elif "freq" in l:
        v[lr[1]] = float(lr[2])
    elif "Final GAMMA" in l and "likelihood" in l:
        v['lk'] = float(lr[-1])
    elif "Overall Time for Tree Evaluation" in l:
        v['tm'] = float(lr[-1])

print('This is a RAxML_info file from a -f e run, automatically reformatted')
print(''.join(ls[0:i]))
print(raxmlinfo % (v['pi(A):'], v['pi(C):'], v['pi(G):'], v['pi(T):'],
                   v['tm'],
                   v['A <-> C:'], v['A <-> G:'], v['A <-> T:'], v['C <-> G:'],
                   v['C <-> T:'], v['G <-> T:'],
                   v['lk']))
