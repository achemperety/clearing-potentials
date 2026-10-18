"""Exhaustive potential searches for the unilateral six-vertex digraphs #03..#12 of Proposition 5.2.

Run from this directory (it imports scf_core.py and scf_groebner.py; requires passagemath-singular):
    python3 unilateral_open.py OUTFILE [first last]
Digraphs #00..#02 were decided in logs/earlier/unilateral6_groebner.txt and #03, #04 in
logs/reruns/unilateral6_groebner_03_04.txt (no potential over any field); #05..#12 are open, and each takes of the order of an
hour or more.
Each digraph is searched over all graphs U up to automorphisms that survive the pruning, in 'clearing' mode
(the weights of vertices with a clearing set fixed by Lemma 3.3), exactly as for #00..#02.
"""
import json
import sys
import time

import scf_groebner as sg
from scf_core import Digraph

out = sys.argv[1]
first = int(sys.argv[2]) if len(sys.argv) > 2 else 3
last = int(sys.argv[3]) if len(sys.argv) > 3 else 12
import os
rows = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'unilateral6_setpot13.json')))
with open(out, 'a') as fh:
    for i in range(first, last + 1):
        D = Digraph(6, [tuple(a) for a in rows[i]])
        t = time.time()
        r = sg.search(D, 'clearing', log=fh)
        r['index'] = i
        r['positive'] = [sorted(map(list, U)) for U in r['positive']]
        fh.write(json.dumps(r) + '\n')
        fh.write('#%02d done in %.0f s: %s\n' % (i, time.time() - t,
                 'no potential over any field' if not r['positive'] else 'POTENTIAL FOUND'))
        fh.flush()
