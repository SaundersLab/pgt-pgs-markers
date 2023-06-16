import os
import re
import subprocess
from collections import abc
from typing import Callable, Dict, Iterable, List, Union
from uuid import uuid4

import numpy as np
import pandas as pd
from Bio.Seq import reverse_complement as rc


class IpcressException(Exception):
    pass

def ipcress_row_to_dict(row: str) -> dict:
    pattern = '(ipcress: )(?P<sequence_id>.*)(:filter\(unmasked\) )(?P<experiment_id>.*) '
    pattern += '(?P<product_length>\d+) (?P<primer_5>[AB]) (?P<pos_5>\d+) (?P<mismatch_5>\d+) '
    pattern += '(?P<primer_3>[AB]) (?P<pos_3>\d+) (?P<mismatch_3>\d+) (?P<description>[a-z]+)'
    match = re.search(pattern, row)
    if match is None:
        raise ValueError("Inavlid IPCRESS row " + row)
    row_str_dict = match.groupdict()
    numeric_substrings = ['pos', 'mismatch', 'length']
    numeric_fields = [
        field for field in row_str_dict 
        if any(s in field for s in numeric_substrings)
    ]
    return {
        field: int(val) if field in numeric_fields else val
        for field, val in row_str_dict.items()
    }


def alignment_block_lines_to_dict(aligment_block: List[str]) -> Dict[str, str]:
    alignment: Dict[str, str] = {}
    fwd, _, primers, _, rev = aligment_block
    alignment['fwd_align_5_to_3'] = fwd[:fwd.index(' # forward')].replace('.', '')
    alignment['rev_align_3_to_5'] = rev[:rev.index(' # revcomp')].replace('.', '')
    alignment['fwd_primer_5_to_3'] = re.findall(".*5'-(.*)-3'", primers)[0]
    alignment['rev_primer_3_to_5'] = re.findall(".*3'-(.*)-5'", primers)[0]
    return alignment


def block_to_dict(block: str) -> dict:
    lines = [line for line in block.split('\n')]
    lines = [line for line in lines if line.strip() and not line.startswith('--')]
    result = ipcress_row_to_dict(lines[11])
    product = ''.join(lines[13:])
    result['product'] = rc(product) if result['description'] == 'revcomp' else product
    return {
        **result,
        **alignment_block_lines_to_dict(lines[6:11]),
    }

def empty_ipcress_df() -> pd.DataFrame:
    cols = [
        'sequence_id', 'experiment_id', 'product_length', 'primer_5', 'pos_5',
        'mismatch_5', 'primer_3', 'pos_3', 'mismatch_3', 'description',
        'product', 'fwd_align_5_to_3', 'rev_align_3_to_5', 'fwd_primer_5_to_3',
        'rev_primer_3_to_5', 'rev_primer_5_to_3', 'rev_align_5_to_3',
        'fwd_pattern_5_to_3', 'rev_pattern_5_to_3', 'would_amplify_5',
        'would_amplify_6', 'would_amplify_7', 'would_amplify_8', 'is_intended'
    ]        
    return pd.DataFrame(columns=cols)

def ipcress_pretty_content_to_df(
    ipcress_pretty_content: str,
) -> pd.DataFrame:
    block_sep = '\nIpcress result\n--------------\n'
    result_blocks = ipcress_pretty_content.split(block_sep)[1:]
    if not result_blocks:    
        return empty_ipcress_df()
    df = pd.DataFrame(map(block_to_dict, result_blocks))
    df.columns = [col.replace(' ', '_').lower() for col in df.columns]
    df['rev_primer_5_to_3'] = df.rev_primer_3_to_5.apply(rc)
    df['rev_align_5_to_3'] = df.rev_align_3_to_5.apply(rc)
    df = df.drop(columns=[col for col in df if col.endswith('3_to_5')])
    return df

def overlap_size(a: List[int], b: List[int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def is_intended_product(
    product: pd.Series,
    intended_products: pd.DataFrame,
    ) -> bool:

    for col in ['p_name', 'sequence_id']:
        if col in intended_products:
            intended_products = intended_products.query(f'{col} == "{product[col]}"')
        if not intended_products.shape[0]:
            return False

    if 'max_product_size_diff' in intended_products:
        intended_product_lengths = intended_products.end - intended_products.start
        product_size_diffs = (intended_product_lengths - product.p_product_size).abs()
        intended_products = intended_products[product_size_diffs <= intended_products.max_product_size_diff]
        if not intended_products.shape[0]:
            return False
    if 'min_n_overlap' in intended_products:
        product_interval = [product.pos_5, product.pos_3]
        overlap_sizes = [
            overlap_size(product_interval, intended_interval)
            for intended_interval in intended_products[['start', 'end']].values
        ]
        overlap_big_enough = [
            overlap > min_overlap 
            for overlap, min_overlap in zip(overlap_sizes, intended_products.min_n_overlap)
        ]
        intended_products = intended_products[overlap_big_enough]
        if not intended_products.shape[0]:
            return False
    return True

def are_intended_products(
    products: pd.DataFrame,
    intended_products: pd.DataFrame,
) -> np.ndarray:
    if not products.shape[0]:
        return np.array([], dtype='bool')
    return products.apply(
        lambda product: is_intended_product(product, intended_products),
        axis=1
    ).values

def comparison_str(primer: str, align: str) -> str:
    if len(primer) != len(align):
        raise ValueError('primer and alignment have different lengths')
    return ''.join('|' if p == a else '.' for p, a in zip(primer, align))

def comparison_strs(primers: str, aligns: str) -> List[str]:
    return [comparison_str(primer, align) for primer, align in zip(primers, aligns)]

def would_amplify(
    product,
    mismatches_3_end=2,
    ignore_mismatches=6,
    size_of_3_end=5,
    mismatch_char='.'
) -> bool:

    for direction in ['fwd', 'rev']:
        comparion_5_to_3 = product[f'{direction}_comparison_5_to_3']
        n_snps = comparion_5_to_3.count(mismatch_char)
        n_snps_3_end = comparion_5_to_3[-size_of_3_end:].count(mismatch_char)

        if n_snps_3_end >= mismatches_3_end or n_snps >= ignore_mismatches:
            return False

    return True

def would_amplify_bias_to_unintended(product) -> bool:
    if product.is_intended:
        return would_amplify(product, size_of_3_end=7)
    return would_amplify(product)

def intended_product_in_products(
    intended_product,
    products,
) -> pd.DataFrame:
    for col in ['p_name', 'sequence_id']:
        if col in intended_product:
            products = products.query(f'{col} == "{intended_product[col]}"')
    if 'max_product_size_diff' in intended_product:
        intended_product_length = intended_product.end - intended_product.start
        product_size_diffs = (intended_product_length - products.p_product_size).abs()
        products = products[product_size_diffs <= intended_product.max_product_size_diff]
    if 'min_n_overlap' in intended_product:
        intended_interval = [intended_product.start, intended_product.end]
        overlap_sizes = [
            overlap_size(product_interval, intended_interval)
            for product_interval in products[['pos_5', 'pos_3']].values
        ]
        overlap_big_enough = [
            overlap > intended_product.min_n_overlap for overlap in overlap_sizes
        ]
        products = products[overlap_big_enough]
    return products
        
def is_intended_product_in_products(
    intended_product,
    products,
) -> bool:
    return intended_product_in_products(intended_product, products).shape[0] > 0

def add_product_summary_to_primer_pairs(
    primer_pairs,
    products,
    intended_products=None,
):
    if intended_products is None:
        intended_products = pd.DataFrame()
    primer_pairs = primer_pairs.copy()
    amplifying_products = products.query('would_amplify')
    n_unintended_products = []
    n_intended_products = []
    for p_name in primer_pairs.p_name:
        p_products = amplifying_products.query(f'p_name == "{p_name}"')
        n_unintended_products.append((~ p_products.is_intended).sum())
        n_intended_products.append(
            sum(
                is_intended_product_in_products(product, p_products)
                for _, product in intended_products.iterrows()
            )
        )
    primer_pairs['n_unintended_products'] = n_unintended_products
    if intended_products.shape[0]:
        primer_pairs['n_intended_products'] = n_intended_products
        # How many products were expected for each primer pair
        n_expected_intended_products = [
            (intended_products.p_name == p_name).sum() if 'p_name' in intended_products else intended_products.shape[0]
            for p_name in primer_pairs.p_name
        ]
        primer_pairs['pct_intended_products'] = [
            100 if n_exp == 0 else 100 * n_obs / n_exp # avoid division by 0 error
            for n_exp, n_obs in zip(n_expected_intended_products, n_intended_products)
        ]
    return primer_pairs


def ipcress(
    assembly_paths: Union[str, Iterable[str]],
    primer_pair_names: Iterable[Union[str, int]],
    fwd_primers: Iterable[str],
    rev_primers: Iterable[str],
    min_product_len: Union[int, Iterable[int]]=None,
    max_product_len: Union[int, Iterable[int]]=4000,
    max_mismatch: int=3,
    memory: int=2048,
    seed: int=12,
    intended_products: pd.DataFrame=None,
    metadata: pd.DataFrame=None,
    f_would_amplify: Callable=would_amplify,
) -> pd.DataFrame:
    
    # max product length of 2 ** 31 or greater crashes IPCRESS
    if isinstance(max_product_len, abc.Iterable):
        if max(max_product_len) >= 2 ** 31:
            raise IpcressException(f'max_product_len exceeding {2 ** 31 - 1} is not supported')
    else:
        if max_product_len >= 2 ** 31:
            raise IpcressException(f'max_product_len exceeding {2 ** 31 - 1} is not supported')
    
    if isinstance(assembly_paths, str):
        assembly_paths = [assembly_paths]
    
    ipcress_uuid = str(uuid4())
    ipcress_in_path = f'{ipcress_uuid}.txt'

    if min_product_len is None:
        min_product_len = [len(f + r) for f, r in zip(fwd_primers, rev_primers)]

    ipcress_in = pd.DataFrame({
        'primer_pair_names': primer_pair_names,
        'fwd_primers': fwd_primers,
        'rev_primers': rev_primers,
        'min_product_len': min_product_len,
        'max_product_len': max_product_len,
    })

    ipcress_cmd = [
        'ipcress', '--input', ipcress_in_path,
        '--pretty', 'true',
        '--mismatch', str(max_mismatch),
        '--products', 'true',
        '--memory', str(memory),
        '--seed', str(seed),
        '--sequence', *assembly_paths
    ]

    ipcress_in.to_csv(ipcress_in_path, sep=' ', header=None, index=None)
    try:
        ipcress_res = subprocess.run(ipcress_cmd, capture_output=True, text=True)
    except Exception as e:
        raise e
    finally:
        os.remove(ipcress_in_path)
    if ipcress_res.returncode:
        raise IpcressException(ipcress_res.stderr)
    products = ipcress_pretty_content_to_df(ipcress_res.stdout)
    # rename columns for consistency with Primer3
    products = products.rename(columns={'experiment_id': 'p_name', 'product_length': 'p_product_size'})
    products = products.sort_values(by=['p_name', 'sequence_id', 'pos_5', 'pos_3']).reset_index(drop=True)

    for direction in ['fwd', 'rev']:
        products.insert(
            list(products).index(f'{direction}_primer_5_to_3') + 1,
            f'{direction}_comparison_5_to_3',
            comparison_strs(products[f'{direction}_primer_5_to_3'], products[f'{direction}_align_5_to_3']),
        )


    if intended_products is not None:
        products['is_intended'] = are_intended_products(products, intended_products)
    else:
        products['is_intended'] = [False] * products.shape[0]

    if products.shape[0]:
        products['would_amplify'] = products.apply(lambda product: f_would_amplify(product), axis=1)
    else:
        products['would_amplify'] = []

    if metadata is not None:
        products = pd.merge(products, metadata, on='sequence_id', how='left')

    return products


def primer_pairs_to_products(
    assembly_paths: Union[str, Iterable[str]],
    primer_pairs: pd.DataFrame,
    min_product_len=None,
    max_product_len=4000,
    max_mismatch=3,
    memory=2048,
    seed=12,
    intended_products: pd.DataFrame=None,
    metadata: pd.DataFrame=None,
    f_would_amplify: Callable=would_amplify,
) -> pd.DataFrame:
    products = ipcress(
        assembly_paths=assembly_paths,
        primer_pair_names=primer_pairs.p_name.values,
        fwd_primers=primer_pairs.l_seq_5_to_3.values,
        rev_primers=primer_pairs.r_seq_5_to_3.values,
        min_product_len=min_product_len,
        max_product_len=max_product_len,
        max_mismatch=max_mismatch,
        memory=memory,
        seed=seed,
        intended_products=intended_products,
        metadata=metadata,
    )
    return products
