import matplotlib.pyplot as plt
from pymsaviz import MsaViz
from Bio.Seq import reverse_complement as rc
from Bio.SeqIO import parse
import pandas as pd
from matplotlib.patches import Rectangle
from Bio import pairwise2
import json
from collections import Counter
import os
from typing import List, Set
results_dir = 'results/samples_msas'
figures_dir = 'figures/samples_msas'
os.makedirs(results_dir, exist_ok=True)
os.makedirs(figures_dir, exist_ok=True)

data_dir = 'data'
primers = pd.read_csv(f'{data_dir}/primers/marker_primers.csv')

metadata = pd.read_csv(f'{data_dir}/metadata/metadata.csv', index_col=0)
metadata = metadata[metadata.Host != 'Oat']
accession_to_name = metadata['Sample Name'].to_dict()
metadata['pgs'] = metadata.Host == 'Rye'
metadata = metadata.sort_values(by=['pgs', 'Sample Name'])

base_to_color = {
    'A': '#B18A51', 'R': '#83BFF1', 'N': '#FFFFFF', 'D': '#29A578', 'C': '#F85604',
    'Q': '#7295AE', 'E': '#2DA0A1', 'G': '#B1C23C', 'H': '#2E94F9', 'I': '#F27663',
    'L': '#DF6E75', 'K': '#7FC3D7', 'M': '#FA9DB0', 'F': '#F9559D', 'P': '#4FA32A',
    'S': '#B4BD9B', 'T': '#D2B576', 'W': '#F92CED', 'Y': '#C96ECF', 'V': '#FA997B',
    'B': '#FFFFFF', 'X': '#FFFFFF', 'Z': '#FFFFFF', 'J': '#FFFFFF', 'O': '#FFFFFF',
    'U': '#FFFFFF', '-': '#FFFFFF'
}

fontsize = 6
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.size'] = fontsize

primers_dict = {}
primer_color = 'C2'

chromosomes = ['chr_8', 'chr_10', 'chr_12', 'chr_18']
for chromosome in chromosomes:
    
    # convert sample accessions to sample names
    accessions = f'{results_dir}/{chromosome}_msa.fasta'
    names = f'{results_dir}/{chromosome}_msa_with_sample_names.fasta'
    records = {r.id: str(r.seq) for r in parse(accessions, 'fasta')}
    with open(names, 'wt') as f:
        for accession, name in zip(metadata.index, metadata['Sample Name']):
            f.write(f'>{name}\n{records[accession].upper()}\n')
    mv = MsaViz(
        names, show_seq_char=False, show_label=(chromosome in {'chr_8', 'chr_12'}),
    )
    mv.set_plot_params(x_unit_size=.01, y_unit_size=.15, ticks_interval=100)

    for record in mv.msa:
        if record.id == 'PGT21-0':
            ref_seq = record.seq
            break

    mv.set_custom_color_func(
        lambda row, col, char, msa: '#EEE' if mv.consensus_seq[col] == char else base_to_color[char]
    )
    top = len(mv.msa)

    fig = mv.plotfig()
    ax = fig.get_axes()[0]
    ax.tick_params('both', length=1, width=1, which='major', labelsize=fontsize, pad=1)

    for _, row in primers[primers.chromosome == chromosome].iterrows():
        direction = row.direction
        seq = rc(row.seq) if direction == 'Reverse' else row.seq

        # where does the primer sequence match the MSA?
        for a in pairwise2.align.localxs(ref_seq, seq, -1, -1):
            start = a.start
            end = a.end
            break

        ax.add_patch(
            Rectangle((start, 0), (end - start), top, color=primer_color, alpha=.03, lw=0)
        )
        
        offset = .07 # hack to account for line thickness
        marker_kwargs = {'color': 'g', 'clip_on': False, 'lw': .9, 'solid_capstyle': 'butt'}
        arrow_color = 'g'
        arrow_far = 7
        arrow_near = .6
        arrow_head_width = .25
        top = len(mv.msa)
        if direction == 'Forward':
            ax.plot([start, end], [top + offset] * 2, **marker_kwargs)
            ax.plot([end - arrow_far, end - arrow_near], [len(mv.msa) + arrow_head_width + offset, len(mv.msa) + offset], **marker_kwargs)
            ax.plot([start, end], [-offset] * 2, clip_on=False, color=arrow_color, lw=.85, solid_capstyle='butt')
            ax.plot([end - arrow_far, end - arrow_near], [-arrow_head_width - offset, -offset], **marker_kwargs)
        else:
            ax.plot([start, end], [top + offset] * 2, **marker_kwargs)
            ax.plot([start + arrow_far, start + arrow_near], [len(mv.msa) + arrow_head_width + offset, len(mv.msa) + offset], **marker_kwargs)
            ax.plot([start, end], [-offset] * 2, **marker_kwargs)
            ax.plot([start + arrow_far, start + arrow_near], [-arrow_head_width - offset, -offset], **marker_kwargs)
        
        # do the text here instead of with add_text_annotation so we only have to loop over
        # the primers once (after we already have the axes)
        ax.text(
            (start + end) / 2, top + .5, {'Forward': 'FWD', 'Reverse': 'REV'}[direction],
            ha='center', fontsize=fontsize, color='g'
        )

        # add data about the primers
        d = row.to_dict()
        d['original_seq'] = d['seq']
        d['seq'] = seq
        d['msa'] = {
            record.id: str(record.seq[start:end]) for record in mv.msa
        }
        primers_dict[row['name']] = d
        unique_bases: List[Set[str]] = list(map(set, seq))
        pgs_samples = {'DE-04', 'DE-05', 'DE-06', 'DE-07'}
        d['n_diffs'] = {}
        for sample, bases in d['msa'].items():
            if (d['amplifies'] == 'Pgt') and (sample in pgs_samples):
                continue
            if len(bases) != len(d['seq']):
                bases = bases.replace('-', '')
            if len(bases) != len(d['seq']):
                print(sample)
                continue
            for i, base in enumerate(bases):
                unique_bases[i].add(base)
            d['n_diffs'][sample] = sum(a != b for a, b in zip(seq, bases))
        if d['n_diffs']:
            d['max_n_diffs'] = max(d['n_diffs'].values())
        d['unique_bases'] = list(map(''.join, unique_bases))
        d['n_variable_sites'] = len([bases for bases in d['unique_bases'] if len(bases) != 1])

    ax.text(
        sum(ax.get_xlim()) / 2,
        top + .75,
        f'M{chromosomes.index(chromosome) + 1}',
        ha='center',
    )
    for t in ax.texts:
        if t.get_text() in {'M1', 'M2', 'M3', 'M4'}:
            t.set_fontsize(7)
        else:
            t.set_fontsize(fontsize)

    fig.set_size_inches(h=3, w=3 + .25 * (chromosome in {'chr_8', 'chr_12'}))
    plt.savefig(f'{figures_dir}/{chromosome}_samples.pdf', bbox_inches='tight', transparent=True)
    plt.close()

    
with open(f'{results_dir}/primer_msa.json', 'wt') as f:
    json.dump(primers_dict, f, indent=4)


n_snps = []
for k, d in primers_dict.items():
    n_snps += list(d['n_diffs'].values())

counts = Counter(n_snps)
df = pd.DataFrame([{'n_mismatches': n, 'n_occurences': counts[n]} for n in sorted(counts)])
df['n_occurences_cum'] = df.n_occurences.cumsum()
df['pct_occurences'] = 100 * df.n_occurences / df.n_occurences.sum()
df['pct_occurences_cum'] = df.pct_occurences.cumsum()
df.to_csv(f'{results_dir}/mismatch_summary.csv', index=None)
