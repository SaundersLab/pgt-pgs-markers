from Bio.SeqIO import parse
import pandas as pd
import sys

ref_fasta = sys.argv[1] # 'data/reference/Pgt_201_B1_AssemblyScaffolds2.fasta'
chromosome_lengths_path = sys.argv[2] # 'data/reference/chromosome_lengths.csv'

pd.DataFrame([
    {'chromosome': record.id, 'length': len(record.seq)}
    for record in parse(ref_fasta, 'fasta')
    if record.id.startswith('chr_')
]).sort_values(by='chromosome').to_csv(chromosome_lengths_path, index=None)
