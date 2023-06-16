from Bio.SeqIO import parse
import sys

for r in parse(sys.argv[1], 'fasta'):
    print('>' + r.id + '\n' + str(r.seq)[::3])
