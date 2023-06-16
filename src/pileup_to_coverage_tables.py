import pandas as pd
from typing import IO
import gzip
import numpy as np
import os
import argparse

def file(path, mode='rt') -> IO:
    """Create a file object from path - infers compression from extension."""
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def pileup_to_coverage_tables(
    pileup: str, out_dir: str, chromosome_lengths: str,
) -> None:
    
    os.makedirs(f'{out_dir}/coverage', exist_ok=True)

    sample = os.path.basename(pileup).split('.')[0]

    chromosome_to_length = pd.read_csv(chromosome_lengths, index_col=0).length.to_dict()

    chromosome_to_coverage = {
        chromosome: np.zeros(length, dtype='int32')
        for chromosome, length in chromosome_to_length.items()
    }

    with file(pileup) as f:
        for line in f:
            seqid, pos, _, coverage, *_ = line.split()
            if seqid not in chromosome_to_coverage:
                continue
            chromosome_to_coverage[seqid][int(pos) - 1] = int(coverage)

    for chromosome, coverage in chromosome_to_coverage.items():
        pd.DataFrame({sample: coverage}).to_csv(
            f'{out_dir}/coverage/{chromosome}.csv', index=None
        )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a table of coverage for each chromosome from a pileup file.')
    parser.add_argument('--pileup', type=str, required=True)
    parser.add_argument('--out_dir', help='Directory to create coverage subdirectory in', type=str, required=True)
    parser.add_argument('--chromosome_lengths', help='Path to csv file of chromosome lengths', type=str, required=True)
    args = vars(parser.parse_args())
    pileup_to_coverage_tables(**args)
