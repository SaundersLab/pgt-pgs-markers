import numpy as np
from Bio.SeqIO import parse
import sys

"""
Filtering the following tree... 

>foo
NNCCNCCCNN
>baz
CNNCCCNNCC
>bar
NNCCNNNNCN

With pct_coverage_threshold = 40 Would give the folowing result:

>foo
CCCN
>baz
NCCC
>bar
CCNC
"""

tree_input = sys.argv[1]
pct_coverage_threshold = int(sys.argv[2])
filtered_tree_input = sys.argv[3]

msa = np.array([list(r.seq) for r in parse(tree_input, 'fasta')])
samples = [r.id for r in parse(tree_input, 'fasta')]
pct_covered = 100 * (msa != 'N').mean(axis=0)
high_coverage_loci = np.flatnonzero(pct_covered >= pct_coverage_threshold)
with open(filtered_tree_input, 'wt') as f:
    for i, sample in enumerate(samples):
        f.write(f'>{sample}\n' + ''.join(msa[i, high_coverage_loci]) + '\n')
