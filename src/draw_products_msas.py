from ipcress import ipcress
import pandas as pd
import glob
from Bio.Seq import reverse_complement as rc
import matplotlib.pyplot as plt
from Bio.SeqIO import parse
import os
import subprocess
from pymsaviz import MsaViz
from matplotlib.patches import Rectangle
from Bio import pairwise2
from collections import Counter
from matplotlib.patches import FancyBboxPatch


results_dir = 'results/assembly_products'
figures_dir = 'figures/assembly_products'
os.makedirs(results_dir, exist_ok=True)
os.makedirs(figures_dir, exist_ok=True)

data_dir = 'data'
assemblies_dir = f'{data_dir}/reference/rust_assembly_sequences'

assembly_paths = glob.glob(f'{assemblies_dir}/*')

# Load the primers
primers = pd.read_csv(f'{data_dir}/primers/primers.csv')
primer_details = primers.name.str.extract(
    'AB_(C)(?P<chromosome>\d+)(?P<direction>[FR])(?P<amplifies>W|COM|)_*(?P<version>2)*'
)
primer_details = primer_details[
    (primer_details.chromosome.isin({'10', '12', '18'})) | (primer_details.version == '2')
]
primer_details = primer_details[['chromosome', 'direction', 'amplifies']]
primers = primers.join(primer_details, how='inner')


# # Remove the contigs and keep only scaffolds from GCA_000149925.1
# # The assembly contains scaffolds, but also individual records for the contigs that
# # comprise the scaffolds. If we run ipcress on this assembly we'll get twice the
# # expected products (one on the scaffold and one on the contig).
# with open(f'{assemblies_dir}/GCA_000149925.1_superconts_only.fa', 'wt') as f:
#     for record in parse(f'{assemblies_dir}/GCA_000149925.1.fa', 'fasta'):
#         if 'supercont' in record.description:
#             f.write(f'>{record.description}\n{str(record.seq)}\n')
# os.rename(
#     f'{assemblies_dir}/GCA_000149925.1_superconts_only.fa',
#     f'{assemblies_dir}/GCA_000149925.1.fa'
# )

# convert assembly accession to a the name for that isolate 
assembly_to_name = {
    'GCF_000149925.1': 'CRL-75-36-700-3',
    'GCA_002762355.2': '99KS76A-1',
    'GCA_008522505.1': '21-0',
    'GCA_008520325.1': 'Ug99',
    'GCA_903797515.1': 'UK-01',
}
# the ipcress results contain the seqid of a product, but not the assembly
# (we give it with multiple assemblies in one run) so we need a way to get
# the assembly from the sequence_id
seqid_to_assembly = {}
seqid_to_description = {}
for path in assembly_paths:
    assembly = path.split('/')[-1][:-3]
    for record in parse(path, 'fasta'):
        if record.id in seqid_to_assembly:
            raise ValueError('sequence_ids are not unique')
        seqid_to_assembly[record.id] = assembly_to_name[assembly]
        seqid_to_description[record.id] = record.description.split(maxsplit=1)[1]

# Find all the products from each pair of primers
for chromosome in ['8', '10', '12', '18']:
    products_tables = []
    chromosome_primers = primers[primers.chromosome == chromosome]
    for fname, forward_primer in chromosome_primers[chromosome_primers.direction == 'F'].iterrows():
        for rname, reverse_primer in chromosome_primers[chromosome_primers.direction == 'R'].iterrows():
            products_tables.append(
                ipcress(
                    assembly_paths=assembly_paths,
                    primer_pair_names=[f'{fname}-{rname}'],
                    fwd_primers=[forward_primer.seq],
                    rev_primers=[reverse_primer.seq],
                    seed=12,
                    max_mismatch=4,
                )
            )
    products = pd.concat(products_tables).reset_index(drop=True)
    # Assume that any matches with 2 or more SNPs at the 3' end of the primer
    # will not amplify - so we drop those products
    products = products[
         (products['fwd_comparison_5_to_3'].str.slice(-5).str.count('\.') < 2) &
         (products['rev_comparison_5_to_3'].str.slice(-5).str.count('\.') < 2)
    ]
    products = products.reset_index(drop=True)
    
    # For a given marker, one amplimer (the Pgt and Pgs product) will be entirely
    # inside of another larger amplimer (the Pgt product). So we have 2 products, one
    # entirely covered by the other, to create the MSA we only use the containing product.
    contained_products = []
    # For each sequence_id (i.e. contig, scaffold, or chromosome depending on the assembly)
    for record, record_products in products.groupby('sequence_id'):
        # get the string of the longest product on this sequence_id
        longest_product = sorted(record_products['product'], key=len)[-1]
        # record the products that are inside the longest product
        for index, product in zip(record_products.index, record_products['product']):
            if (product != longest_product) and (product in longest_product):
                contained_products.append(index)
    
    # To create the MSA we want the products that are not contained 
    # in a larger product on the same sequence_id
    non_contained_products = products.drop(contained_products)
    non_contained_products.to_csv(f'{results_dir}/chr_{chromosome}_products.csv', index=None)
    names = non_contained_products.sequence_id.values
    sequences = [
        rc(p) if desc == 'revcomp' else p 
        for p, desc in non_contained_products[['product', 'description']].values
    ]
    
    with open(f'{results_dir}/chr_{chromosome}_products.fasta', 'wt') as f:
        for name, sequence in zip(names, sequences):
            f.write(f'>{name}\n{"".join(sequence)}\n')

# create an MSA for each set of non-contained products
for path in glob.glob(f'{results_dir}/chr_*_products.fasta'):
    chromosome = path.split('/')[-1].rsplit('_', 1)[0]
    subprocess.run(
        ['muscle', '-align', path, '-output', f'{results_dir}/{chromosome}_msa.fasta']
    )



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


def seqid_to_label(seqid: str) -> str:
    """
    Generate a short meanigful label from a sequence id
    """
    assembly = seqid_to_assembly[record.id]
    description = seqid_to_description[record.id]
    for s in (
        'Puccinia graminis f. sp. tritici ',
        ', whole genome shotgun sequence.',
        'isolate RKQQC ',
        'strain Ug99 ',
        'CRL 75-36-700-3 ',
        'genome assembly, contig: ',
        'isolate 21-0 ',
        ' genomic scaffold',
        '_pilon5',
        ', whole genome shotgun sequence',
        'omosome ',
    ):
        description = description.replace(s, '')
    return assembly + ' ' + description

fontsize = 6
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.size'] = fontsize

# convert the sequence ids to labels
chromosomes = ['chr_8', 'chr_10', 'chr_12', 'chr_18']
for chromosome in chromosomes:
    accessions = f'{results_dir}/{chromosome}_msa.fasta'
    labelled_records = []
    for record in parse(accessions, 'fasta'):
        labelled_records.append((seqid_to_label(record.id), str(record.seq)))
    labelled_records = sorted(labelled_records, key=lambda l: l[0].lower())
    labelled = f'{results_dir}/{chromosome}_msa_labelled.fasta'
    with open(labelled, 'wt') as f:
        for label, seq in labelled_records:
            f.write(f'>{label}\n{seq}\n')

name_to_color = {
    '21-0': 'C0',
    '99KS76A-1': 'C1',
    'CRL-75-36-700-3': 'C6',
    'Ug99': 'C8',
    'UK-01': 'C5',
}
# Which marker to put next to the contig/scaffole/chromosome in the plot
name_to_symbol = {
    '21-0': '*',
    '99KS76A-1': 's',
    'CRL-75-36-700-3': 'o',
    'Ug99': '^',
    'UK-01': 'P',
}


def get_text_bbox(ax, t):
    # get the inverse of the transformation from data coordinates to pixels
    bb = t.get_window_extent()
    bb_datacoords = bb.transformed(ax.transData.inverted())
    return bb_datacoords

primer_color = 'C2'

n_snps = []
for chromosome in chromosomes:
    labelled = f'{results_dir}/{chromosome}_msa_labelled.fasta'

    mv = MsaViz(labelled, show_seq_char=False, show_label=True, label_type='description')
    mv.set_plot_params(x_unit_size=.01, y_unit_size=.15, ticks_interval=100)
    ref_seq = mv.consensus_seq
    # show agreement with consensus as #EEE
    mv.set_custom_color_func(
        lambda row, col, char, msa: '#EEE' if ref_seq[col] == char else base_to_color[char]
    )

    fig = mv.plotfig()
    top = len(mv.msa)
    fig.set_size_inches(h=top / 12, w=2.8)
    ax = fig.get_axes()[0]

    # make the ticks short and move the numbers right up to the MSA grid
    ax.tick_params('both', length=1, width=1, which='major', labelsize=fontsize, pad=1)

    for _, row in primers[primers.chromosome == chromosome].iterrows():
        direction = row.direction
        seq = rc(row.seq) if direction == 'Reverse' else row.seq

        # where does the primer sequence match the MSA?
        for a in pairwise2.align.localxs(ref_seq, seq, -1, -1):
            start = a.start
            end = a.end
            break        

        # highlight the primer region
        ax.add_patch(
            Rectangle((start, 0), (end - start), top, color=primer_color, alpha=.03, lw=0)
        )
        
        # offset is a hack to make sure the primer arrow
        # doesn't cover the MSA (it accounts for line thickness)
        offset = .07 

        marker_kwargs = {'color': primer_color, 'clip_on': False, 'lw': .9, 'solid_capstyle': 'butt'}
        # extremes of the arrow head
        arrow_far = 7
        arrow_near = .6
        arrow_head_width = .25
        
        # Draw an arrow at the top and bottom of the MSA grid for each primer
        ax.plot([start, end], [top + offset] * 2, **marker_kwargs)
        ax.plot([start, end], [-offset] * 2, **marker_kwargs)
        if direction == 'Forward':
            ax.plot([end - arrow_far, end - arrow_near], [top + arrow_head_width + offset, top + offset], **marker_kwargs)
            ax.plot([end - arrow_far, end - arrow_near], [-arrow_head_width - offset, -offset], **marker_kwargs)
        else:
            ax.plot([start + arrow_far, start + arrow_near], [top + arrow_head_width + offset, top + offset], **marker_kwargs)
            ax.plot([start + arrow_far, start + arrow_near], [-arrow_head_width - offset, -offset], **marker_kwargs)
        
        # adding the text here instead of with add_text_annotation means we only
        # have to loop over the primers once (after we already have the axes)
        ax.text(
            (start + end) / 2,
            top + .5, {'Forward': 'FWD', 'Reverse': 'REV'}[direction],
            ha='center',
            fontsize=fontsize,
            color=primer_color,
        )

        for record in mv.msa:
            n_snps.append(sum(a != b for a, b in zip(seq, record.seq[start:end])))

    # label which marker this is
    ax.text(
        sum(ax.get_xlim()) / 2,
        top + .75,
        f'M{chromosomes.index(chromosome) + 1}',
        ha='center',
    )

    # get the width of a space character (used for calculations later)
    space = ax.text(0, 0, ' ', fontsize=fontsize)
    bbox = get_text_bbox(ax, space)
    space_width = bbox.xmax - bbox.xmin
    space.remove()

    for t in ax.texts:
        # if this is a label for a product then set it to the appropriate color
        # for the assembly it came from
        t.set_color(name_to_color.get(t.get_text().split()[0], t.get_color()))

        if t.get_text() in {'M1', 'M2', 'M3', 'M4'}:
            t.set_fontsize(7)
        else:
            t.set_fontsize(fontsize)

        # if this text is the label for a product (i.e "assembly_name some description of sequence")  
        if ' ' in t.get_text():
            assembly_name, label = t.get_text().split(maxsplit=1)
            t.set_text(label)
            bbox = get_text_bbox(ax, t)
            # put a marker to the left of the label
            ax.scatter(
                bbox.xmin - 1.7 * space_width,
                .05 + (bbox.ymin + bbox.ymax) / 2,
                clip_on=False,
                zorder=100,
                s=11,
                marker=name_to_symbol[assembly_name],
                color=name_to_color[assembly_name],
            )

    plt.savefig(f'{figures_dir}/{chromosome}_assembly_products.pdf', bbox_inches='tight', transparent=True)
    plt.close()

# create a legend for the figure (which marker and color for which assembly)
fig, ax = plt.subplots(1, 1, figsize=(3.9, .22))
ax.set_xticks([])
ax.set_yticks([])
x = .03
ax.set_xlim(0, 1)
legend_to_color = {'Assembly:': 'k', **name_to_color}
space = ax.text(0, 0, ' ', fontsize=fontsize)
bbox = get_text_bbox(ax, space)
space_width = bbox.xmax - bbox.xmin
space.remove()
for i, (name, color) in enumerate(legend_to_color.items()):
    t = ax.text(x, .45, name, color=color, fontsize=fontsize, va='center', ha='left')
    plt.gcf().canvas.draw()
    if i == 0:
        t.set_fontweight('bold')
    else:
        ax.scatter(x - 1.7 * space_width, .5, s=11, marker=name_to_symbol[name], color=color)
    bbox = get_text_bbox(ax, t)
    x = bbox.xmax + .05

ax.set_ylim(0, 1)
ax.spines[['right', 'top', 'left', 'bottom']].set_visible(False)

size = fig.get_size_inches()
aspect = size[0] / size[1]
# use the same aesthetic as the default matplotlib legend
ax.add_patch(
    FancyBboxPatch(
        (.015,.255), .97,.5,
        boxstyle="round,pad=0.008", mutation_aspect=aspect, ec='#D5D5D5', fc='none'
    )
)

plt.savefig(f'{figures_dir}/assembly_products_legend.pdf', bbox_inches='tight', transparent=True)
plt.close()


# summarise the snps between the assembly and the primers
counts = Counter(n_snps)
df = pd.DataFrame([{'n_mismatches': n, 'n_occurences': counts[n]} for n in sorted(counts)])
df['n_occurences_cum'] = df.n_occurences.cumsum()
df['pct_occurences'] = 100 * df.n_occurences / df.n_occurences.sum()
df['pct_occurences_cum'] = df.pct_occurences.cumsum()
df.to_csv(f'{results_dir}/assembly_products_mismatch_summary.csv', index=None)
