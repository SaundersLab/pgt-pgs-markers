import colorsys
import sys
from typing import Any, Dict, List, Tuple, Union

import matplotlib
import matplotlib.colors as mc
import matplotlib.image as mpimg
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.Phylo import Newick, draw, read
from matplotlib.figure import Figure


def get_parent(tree: Newick.Tree, child: Union[str, Newick.Clade]) -> Newick.Clade:
    if isinstance(child, str):
        child = list(tree.find_clades(child))[0]
    return tree.get_path(child)[-2]


def load_metadata(path):
    metadata = pd.read_csv(path).set_index('Sample Accession')
    # Abbreviate because there's little space in the metadata column of the plot
    metadata['Country'] = metadata['Country'].str.replace(
        'South Africa', 'S. Africa'
    )
    # Distinction between wheat and rye is more important than wheat and bread wheat
    metadata['Host'] = metadata['Host'].str.replace('Bread wheat', 'Wheat')
    return metadata


def load_tree(path: str) -> Newick.Tree:
    tree = read(path, 'newick')
    tree.root_at_midpoint()
    tree.ladderize(reverse=True)
    return tree


def get_boundaries(l: List[Tuple[str, float, float]]) -> List[list]:
    """
    Find the boundaries of contiguous labels in a metadata column
    from the value and position of the texts.

    Examples
    --------
    >>> get_boundaries([
        ('A', .5, 1),
        ('A', 20, 2),
        ('A', 9, 3),
        ('B', .1, 4),
        ('C', 5, 5),
        ('D', 5, 6),
        ('D', 5, 7),
        ('F', 5, 8),
    ])
    [['A', 1, 3], ['B', 4, 4], ['C', 5, 5], ['D', 6, 7], ['F', 8, 8]]
    """
    curr_val = l[0][0]
    curr_y = l[0][2]
    boundaries = [[curr_val, curr_y]]
    for val, x, y in l[1:]:
        if val != curr_val:
            boundaries[-1].append(curr_y)
            boundaries.append([val, y])
        curr_val = val
        curr_y = y
    boundaries[-1].append(curr_y)
    return boundaries


def color_hls(color):
    try:
        c = mc.cnames[color]
    except:
        c = color
    return colorsys.rgb_to_hls(*mc.to_rgb(c))


def remove_axis_decorations(ax: plt.Axes):
    for side in ['top', 'left', 'bottom', 'right']:
        ax.spines[side].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(None)
    ax.set_ylabel(None)


def set_branch_thickness(tree_ax: plt.Axes, thickness: float):
    for collection in tree_ax.collections:
        paths = collection.get_paths()
        path = paths[0]
        extents = path.get_extents()
        if extents.x1 == 0:
            # The very first vertical line would be
            # split in half when we set xlim to (0, ...) so
            # compensate for that here doubling its thickness
            collection.set_linewidths([thickness * 2])
        else:
            collection.set_linewidths([thickness])


def replace_confidence_numbers_with_markers(
    tree_ax: plt.Axes, confidence_threshold=80,
):
    """
    Phylo.draw puts confidence numbers as text at the midpoint of the branch.
    For a visually simpler summary, replace that text with a marker at the
    end of the branch if the confidence is below some threshold. 
    """

    # Create a mapping from midpoint of the branch (i.e. where the condifidence
    # number will be draw by Phylo.draw) to the right most position of the branch
    # (i.e. where we will put the confidence marker)
    midpoint_to_line_end = {}
    for collection in tree_ax.collections:
        paths = collection.get_paths()
        path = paths[0]
        extents = path.get_extents()
        if extents.height or (not extents.width):
            continue
        xmid = round(extents.x0 + extents.width / 2, 5)
        midpoint_to_line_end[(xmid, round(extents.y0, 5))] = (extents.x1, extents.y1)

    bootstrap_marker_props: Dict[str, Any] = dict(
        marker='o',
        zorder=10,
        # color=plt.get_cmap('Set1')(4),
        color=plt.get_cmap('Dark2')(7),
        edgecolors='k',
        s=15,
        lw=.5,
    )

    for t in tree_ax.texts:
        # When Phylo.draw creates the node labels they all begin with spaces.
        # This is not true of confidence numbers. So skip any text starting with ' '
        if not t.get_text().startswith(' '):
            try:
                # double check that this text is indeed a confidence score
                confidence = int(t.get_text())
                if confidence < 0 or confidence > 100:
                    raise ValueError('Not confidence')

                # If the confidence is of concern (i.e. below the threshold) then
                # put a marker at the right end of the branch
                if confidence < confidence_threshold:
                    pos: Tuple[Any, Any] = (
                        round(t.get_position()[0], 5),
                        round(t.get_position()[1], 5),
                    )                    
                    # tuple(map(lambda v: round(v, 5), t.get_position()))
                    tree_ax.scatter(
                        *midpoint_to_line_end[pos],
                        **bootstrap_marker_props,
                    )
                # Hide the confidence number
                t.set_visible(False)
            except:
                # This was not actually a confidence number
                pass

    tree_ax.scatter(
        [],
        [],
        label=f'Bootstrap < {confidence_threshold}',
        **bootstrap_marker_props
    )
    tree_ax.legend(loc='center left')


def draw_field_table(
    meta_axes,
    table_fields,
    field_to_val_to_color,
    metadata,
    sample_id_positions,
    col_margin,
    row_margin=.12,
    # forma specialis should be written italicised
    should_italicise_field=lambda field: field == 'fsp',
):
    for field_index, (field, field_ax) in enumerate(zip(table_fields, meta_axes)):
        val_to_color = field_to_val_to_color[field]
        value_positions = [
            (metadata.loc[sample_id, field], x, y)
            for sample_id, x, y
            in sample_id_positions
        ]
        value_boundaries = get_boundaries(value_positions)
        for val, start, end in value_boundaries:
            color = val_to_color[val]
            hls = color_hls(color)
            # light foreground for dark background (vice versa)
            font_color = '#BBB' if hls[1] < .4 else '#000'
            mid = (start + end) / 2
            height = end - start + 1

            # For balance, shift words without descenders down very slightly
            has_descenders = any(letter in str(val) for letter in 'gjpqyQ')
            field_ax.text(
                0,
                mid + (.08 if not has_descenders else .02),
                ('$\it{' + str(val) +
                 '}$') if should_italicise_field(field) else str(val),
                color=font_color,
                ha='center',
                va='center'
            )

            if field_index:
                width = 2 - 2 * col_margin
                rect_x = -1 + col_margin
            else:
                # The very first column is the label for a clade and so do not
                # insert a horizontal margin to the left
                width = 2 - col_margin
                rect_x = -1
            field_ax.add_patch(patches.Rectangle(
                (rect_x, start - .5 + row_margin),
                width, height - 2 * row_margin, color=color)
            )
        field_ax.set_xlim(-1, 1)


def add_scale_bar(
    ax: plt.Axes,
    lw=.75,
    scale_bar_width=.01,
    scale_bar_left=0,
    scale_bar_y=0,
    scale_bar_end_bar_height=.8,
):
    scale_bar_right = scale_bar_left + scale_bar_width
    scale_bar_kwargs = {'lw': .75, 'color': 'k'}
    ax.hlines(y=scale_bar_y, xmin=scale_bar_left,
              xmax=scale_bar_right, color='k', lw=lw)
    # draw the end bars
    for x in [scale_bar_left, scale_bar_right]:
        ax.vlines(
            x=x,
            ymin=scale_bar_y - scale_bar_end_bar_height / 2,
            ymax=scale_bar_y + scale_bar_end_bar_height / 2,
            color='k',
            # double the thickness to account for half the line being
            # cut off if the scale bar starts at 0
            lw=lw * (2 if x == 0 else 1)
        )
    ax.text(
        scale_bar_left + scale_bar_width / 2,
        scale_bar_y + .4,
        scale_bar_width,
        va='top',
        ha='center'
    )


def get_rightmost_end_of_any_node_label(fig: Figure, tree_ax: plt.Axes) -> float:
    # Set the xlim to match the right most position of any node label
    rightmost_text_end = 0
    for t in tree_ax.texts:
        bbox = t.get_tightbbox(fig.canvas.get_renderer()).transformed(
            tree_ax.transData.inverted())
        rightmost_text_end = max(rightmost_text_end, bbox.x1)
    return rightmost_text_end


def change_sample_ids_to_names(tree_ax: plt.Axes, metadata: pd.DataFrame):
    for t in tree_ax.texts:
        if not t.get_text().startswith(' '):
            continue
        sample_id = t.get_text()[1:]
        x, y = t.get_position()
        t.set_text(' ' + metadata.loc[sample_id, 'Sample Name'])
        t.set_position((x, y + .1))


def draw_tree_with_metadata(
    tree: Newick.Tree,
    metadata: pd.DataFrame,
    field_to_val_to_color: Dict[str, Dict[str, str]],
    width=1,
    col_margin=.05,
    row_margin=.12,
    confidence_threshold=80,
    figsize=(6.8, 6),
):
    # fields to display with colors next to the tree
    table_fields = [col for col in metadata if col != 'Sample Name']
    n_table_fields = len(table_fields)
    fig, axes = plt.subplots(
        1, 1 + n_table_fields,
        gridspec_kw={'width_ratios': [6] + [1] * n_table_fields},
        figsize=figsize,
    )
    tree_ax = axes[0]
    meta_axes = axes[1:]
    plt.subplots_adjust(wspace=0, hspace=0)
    draw(tree, axes=tree_ax, do_show=False, show_confidence=True)

    # Create a list of tuples [(sample_id, x, y), ..., (sample_id, x, y)]
    # using the fact that Phylo.draw positions a text element where the
    # text is the sample name and a space, at the position of the node
    sample_id_positions = [
        (t.get_text()[1:], ) + t.get_position()
        for t in tree_ax.texts
        if t.get_text().startswith(' ')
    ]
    # Sort by y position
    sample_id_positions = sorted(
        sample_id_positions,
        key=lambda sample_id_x_y: sample_id_x_y[-1]
    )

    draw_field_table(
        meta_axes,
        table_fields,
        field_to_val_to_color,
        metadata,
        sample_id_positions,
        col_margin,
        row_margin=row_margin
    )
    set_branch_thickness(tree_ax, .75)

    ylim = tree_ax.get_ylim()
    ylim = (
        ylim[0] + 2,  # Make space for the scale bar at the bottom
        ylim[1] - .5
    )
    for ax in axes:
        ax.set_ylim(ylim)

    replace_confidence_numbers_with_markers(tree_ax, confidence_threshold)

    # fill the clades of the first field of the metadata
    # TODO: should probably throw an error if we the different values
    #       for the field create overlapping rectangles
    field = metadata.columns[1]
    tree_ax_width = (tree_ax.get_xlim()[1] - tree_ax.get_xlim()[0])
    clade_fill_extra_left = .01 * tree_ax_width
    for val, val_sample_accessions in metadata.reset_index().groupby(field)['Sample Accession'].unique().to_dict().items():
        val_min_y = 2 ** 32
        val_max_y = - 2 ** 32
        for sample_id, x, y in sample_id_positions:
            if sample_id in val_sample_accessions:
                val_min_y = min(val_min_y, y)
                val_max_y = max(val_max_y, y)
        mrca = tree.is_monophyletic(
            [list(tree.find_clades(s))[0] for s in val_sample_accessions])
        left = sum([c.branch_length for c in tree.get_path(mrca)]
                   ) - clade_fill_extra_left
        color = field_to_val_to_color[field][val]
        tree_ax.add_patch(patches.Rectangle(
            (left, val_min_y + row_margin - .5),
            width=50 * tree_ax_width,  # has to fill all the way up to the metadata field
            height=(val_max_y - val_min_y + 1 - 2 * row_margin),
            color=color,
        ))

        # # If you want to make the labels and branches in a clade a different
        # # shade of the background color then unccoment this but it gives poor contrast
        # color_to_dark_color = {
        #     plt.get_cmap('Pastel1')(2): plt.get_cmap('Set1')(2),
        #     plt.get_cmap('Pastel2')(5): plt.get_cmap('Dark2')(5),
        #     plt.get_cmap('Pastel1')(1): plt.get_cmap('Set1')(1),
        # }
        # for t in tree_ax.texts:
        #     if t.get_text()[1:] in val_sample_accessions:
        #         t.set_color(color_to_dark_color[color])
        # for collection in tree_ax.collections:
        #     paths = collection.get_paths()
        #     path = paths[0]
        #     extents = path.get_extents()
        #     if ((val_min_y - 1) < extents.y0 <= extents.y1 < val_max_y + 1) and (left <= extents.x0):
        #         collection.set_color(color_to_dark_color[color])

    # # to center the scale bar in the tree axis:
    # scale_bar_left = .5 * (tree_ax.get_xlim()[0] + tree_ax.get_xlim()[1] - scale_bar_width)
    add_scale_bar(tree_ax, scale_bar_y=ax.get_ylim()[0] - 1.5)

    for ax in axes:
        remove_axis_decorations(ax)

    change_sample_ids_to_names(tree_ax, metadata)
    tree_ax.set_xlim(0, get_rightmost_end_of_any_node_label(fig, tree_ax))

    return fig, axes


if __name__ == '__main__':

    metadata_path = sys.argv[1]  # 'data/metadata/metadata.csv'
    tree_path = sys.argv[2]  # 'results/tree/msa_80pct_covered.raxml.support'
    out_prefix = sys.argv[3]  # 'figures/tree'

    pastels = (
        [plt.get_cmap('Pastel1')(i / 9) for i in [1]] +
        [plt.get_cmap('Pastel2_r')(i / 8) for i in [2]] +
        [plt.get_cmap('Pastel1')(i / 9) for i in [2]] +
        [plt.get_cmap('Pastel1')(i / 9) for i in range(9) if i not in [1, 2]][::-1] +
        [plt.get_cmap('Pastel2_r')(i / 8) for i in range(8) if i not in [2]]
    )

    metadata = load_metadata(metadata_path)
    tree = load_tree(tree_path)

    metadata['fsp'] = metadata['Host'].apply(
        lambda s: {'Wheat': 'Pgt', 'Rye': 'Pgs', 'Oat': 'Pga'}[s]
    )

    i = 0
    field_to_val_to_color = {}
    for field in ['fsp', 'Country']:
        val_to_color = {}
        for val in sorted(set(metadata[field])):
            val_to_color[val] = pastels[i]
            i += 1
        field_to_val_to_color[field] = val_to_color

    new_colors = {
        'Mexico': plt.get_cmap('Pastel2')(7),
        'Uruguay': plt.get_cmap('Pastel1')(4),
        'Iran': plt.get_cmap('Pastel2')(2),
        'Sweden': plt.get_cmap('Pastel1')(0),
    }
    for country, color in new_colors.items():
        field_to_val_to_color['Country'][country] = color

    fontsize = 7
    matplotlib.rc(
        'font',
        **{
            'size': fontsize,
            'family': 'sans-serif',
            'sans-serif': ['Helvetica']
        }
    )

    figsize=(7, 6)
    fig, axes = draw_tree_with_metadata(
        tree,
        metadata[['Sample Name', 'fsp', 'Country']],
        field_to_val_to_color,
        1,
        confidence_threshold=70,
        figsize=figsize
    )
    # Save figure twice, the second version will be changed into grayscale
    for with_color in [True, False]:
        plt.savefig(
            f'{out_prefix}{".svg" if with_color else "_gray.png"}',
            bbox_inches='tight',
            pad_inches=0,
        )
    plt.close()
    # make a grayscale version to check it's clear without color
    def rgb2gray(rgb):
        return np.dot(rgb[..., :3], [0.2989, 0.5870, 0.1140])
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    remove_axis_decorations(ax)
    ax.imshow(rgb2gray(mpimg.imread(f'{out_prefix}_gray.png')), cmap='gray')
    plt.savefig(f'{out_prefix}_gray.png', bbox_inches='tight', pad_inches=0)
