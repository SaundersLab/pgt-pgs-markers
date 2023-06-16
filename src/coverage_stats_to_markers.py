"""
Find regions of high Pgs coverage that are near to regions of low
Pgs coverage, such that both regions overlap a single region of
high Pgt coverage.

Arguments:
1 Path to the reference genome in FASTA format
2 A directory containing a table for each chromosome, with the columns:
    - pgs_min: minimum covergae of any Pgs sample at each position
    - pgs_max: maximum covergae of any Pgs sample at each position
    - pgt_min: minimum covergae of any Pgt sample at each position
3 Path to a spreadsheet to write markers to

Example usage:
python3 src/coverage_stats_to_markers.py \
    data/reference/Pgt_201_B1_AssemblyScaffolds2.fasta\
    data/coverage_stats \
    data/markers.xlsx
"""

from Bio.SeqIO import parse
import numpy as np
import pandas as pd
from pandas.arrays import IntervalArray
from typing import Iterable
import sys


def true_regions(
    l: Iterable, min_len=0, scale_factor=1
) -> IntervalArray:
    l = np.array([False, *l, False])
    starts = scale_factor * np.flatnonzero(l[1:] > l[:-1])
    ends = scale_factor * np.flatnonzero(l[1:] < l[:-1])
    true_intervals = IntervalArray.from_arrays(starts, ends, 'left')
    return true_intervals[true_intervals.length >= min_len]


def test_true_regions():
    l = [False, True, True, False, False, True]
    assert (true_regions(l) == IntervalArray.from_arrays(
        [1, 5], [3, 6], 'left')).all()
    assert (true_regions(l, scale_factor=10) ==
            IntervalArray.from_arrays([10, 50], [30, 60], 'left')).all()
    assert (true_regions(l, min_len=2) ==
            IntervalArray.from_arrays([1], [3], 'left')).all()
    assert (true_regions(l, min_len=20, scale_factor=10) ==
            IntervalArray.from_arrays([10], [30], 'left')).all()


def pad_interval(interval: pd.Interval, n: float) -> pd.Interval:
    return pd.Interval(interval.left - n, interval.right + n, interval.closed)


def test_pad_interval():
    padded = pad_interval(pd.Interval(20, 50, 'left'), 10)
    assert padded.left == 10
    assert padded.right == 60
    assert padded.closed == 'left'


test_true_regions()
test_pad_interval()


def find_markers_on_chromosome(
    pgs_min: np.ndarray,
    pgs_max: np.ndarray,
    pgt_min: np.ndarray,
    chromosome_ref_seq: str,
    chunk_size: int = 20,
    max_dist: int = 400,
    min_region_length: int = 150,
) -> pd.DataFrame:
    """
    Find regions of high Pgs coverage that are near to regions of low
    Pgs coverage, such that both regions overlap a single region of
    high Pgt coverage.

    Parameters
    ----------
    pgs_min : ndarray
        1-D array of the minimum coverage of any Pgs sample at each position
        on the chromosome.
    pgs_max : ndarray
        1-D array of the maximum coverage of any Pgs sample at each position
        on the chromosome.
    pgt_min : ndarray
        1-D array of the minimum coverage of any Pgt sample at each position
        on the chromosome.
    chromosome_ref_seq: str
        The reference sequence of the chromosome. Used to create the template
        sequence within which to design primers for the marker.
    chunk_size: int
        When calculating summary statistics on the pgs_min, pgs_max, and
        pgt_min values, how big should the chunks be.
    max_dist: int
        How for apart on the chromosome can regions of high and low Pgs
        coverage be when forming a marker.
    min_region_length: int
        How many bases long must a high or low Pgs coverage region be to
        form part of a marker.
    """

    assert pgs_min.size == pgs_max.size == pgt_min.size

    n_chunks = pgs_min.size // chunk_size

    # A low pgs chunk is one where no sample had any coverage at all.
    pgs_max = pgs_max[:n_chunks * chunk_size].reshape((n_chunks, chunk_size))
    pgs_low_regions = true_regions(
        pgs_max.max(axis=1) == 0,
        min_len=min_region_length,
        scale_factor=chunk_size,
    )

    # Given the minimum coverage of any pgs sample at each position
    # on the chromosome, group the values into chunks of "chunk_size" (e.g. 20).
    # If the 10th percentile of a chunk is at least 3, then we call the chunk
    # high. Contiguous high chunks (spanning at least "min_region_length" bases)
    # form a high pgs region.
    pgs_min = pgs_min[:n_chunks * chunk_size].reshape((n_chunks, chunk_size))
    pgs_high_regions = true_regions(
        np.percentile(pgs_min, 10, axis=1) >= 3,
        min_len=min_region_length,
        scale_factor=chunk_size,
    )

    # Defined the same as high pgs regions but for pgt.
    pgt_min = pgt_min[:n_chunks * chunk_size].reshape((n_chunks, chunk_size))
    pgt_high_regions = true_regions(
        np.percentile(pgt_min, 10, axis=1) >= 3,
        min_len=2 * min_region_length,
        scale_factor=chunk_size,
    )

    markers_rows = []
    for pgs_high in pgs_high_regions:
        # If a padded pgs_high interval overlaps a pgs_low interval, then the
        # distance between the original pgs_high interval and the pgs_low
        # interval is at most the padding size.
        padded_pgs_high = pad_interval(pgs_high, max_dist)
        for pgs_low in pgs_low_regions[pgs_low_regions.overlaps(padded_pgs_high)]:
            pgs_left, pgs_right = sorted([pgs_high, pgs_low])
            # A pgt high must cover at least this much of the pgs low/high pair
            minimal_marker = pad_interval(
                pd.Interval(pgs_left.right, pgs_right.left, 'left'),
                min_region_length
            )
            for pgt_high in pgt_high_regions:
                if pgt_high.left <= minimal_marker.left <= minimal_marker.right <= pgt_high.right:
                    markers_rows.append({
                        'pgs_lo': 'left' if pgs_low < pgs_high else 'right',
                        'pgs_hi_l': pgs_high.left,
                        'pgs_hi_r': pgs_high.right,
                        'pgs_lo_l': pgs_low.left,
                        'pgs_lo_r': pgs_low.right,
                        'pgt_hi_l': pgt_high.left,
                        'pgt_hi_r': pgt_high.right,
                    })

    if not markers_rows:
        return pd.DataFrame()

    markers = pd.DataFrame(markers_rows)
    # Extract the reference sequence for the left most and right most positions
    # related to this marker - we'll use these this for primer design
    markers['template_l'] = markers[markers.columns[-6:]].min(axis=1)
    markers['template_r'] = markers[markers.columns[-6:]].max(axis=1)
    markers['template'] = [
        chromosome_ref_seq[l:r]
        for l, r in markers[['template_l', 'template_r']].values
    ]
    # Convert from chromosome coordinates to template coordindates to help
    # with primer design.
    coordinate_cols = [
        'pgs_hi_l', 'pgs_hi_r', 'pgs_lo_l', 'pgs_lo_r', 'pgt_hi_l', 'pgt_hi_r'
    ]
    for col in coordinate_cols:
        markers[f'template_{col}'] = markers[col] - markers['template_l']

    return markers


def main(
    ref_path: str,
    coverage_stats_dir: str,
    markers_out: str,
):
    markers_dataframes = []
    ref = {r.id: str(r.seq) for r in parse(ref_path, 'fasta')}
    print('processing chromosome')
    # Find markers on each chromosome and then concatenate them
    for chromosome in range(1, 19):
        print(chromosome)
        chromosome_coverage_stats = pd.read_csv(
            f'{coverage_stats_dir}/chr_{chromosome}.csv')
        chromosome_ref_seq = ref[f'chr_{chromosome}']
        chromosome_markers = find_markers_on_chromosome(
            chromosome_coverage_stats.pgs_min.values,
            chromosome_coverage_stats.pgs_max.values,
            chromosome_coverage_stats.pgt_min.values,
            chromosome_ref_seq=chromosome_ref_seq,
        )
        if chromosome_markers.shape[0]:
            chromosome_markers.insert(0, 'chromosome', f'chr_{chromosome}')
            markers_dataframes.append(chromosome_markers)
    if markers_dataframes:
        markers = pd.concat(markers_dataframes)
    else:
        markers = pd.DataFrame(
            columns=[
                'chromosome', 'pgs_lo', 'pgs_hi_l', 'pgs_hi_r', 'pgs_lo_l', 'pgs_lo_r', 'pgt_hi_l', 'pgt_hi_r',
                'template_l', 'template_r', 'template', 'template_pgs_hi_l', 'template_pgs_hi_r',
                'template_pgs_lo_l', 'template_pgs_lo_r', 'template_pgt_hi_l', 'template_pgt_hi_r'
            ]
        )
    markers.to_csv(markers_out, index=None)


if __name__ == '__main__':
    ref_path = sys.argv[1]
    coverage_stats_dir = sys.argv[2]
    markers_out = sys.argv[3]
    main(ref_path, coverage_stats_dir, markers_out)
