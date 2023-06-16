import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
from matplotlib.patches import Rectangle
import matplotlib.ticker as ticker

def primer_poly_path(direction, start, stop, y, height=1, point_length=12):
    """
    Draws an optinally pointy rectangle, could represent a primer
    or even an amplimer
    """

    if direction == 'F':
        return [
            [start, y],
            [start, y + height],
            [stop - point_length, y + height],
            [stop, y + height / 2],
            [stop - point_length, y],
            [start - .05, y],
        ]
    return [
        [start + point_length, y],
        [start, y + height / 2],
        [start + point_length, y + height],
        [stop, y + height],
        [stop, y],
        [start + point_length, y],
    ]

if __name__ == '__main__':

    primers = pd.read_csv('data/primers/marker_primers.csv')
    markers = pd.read_csv('data/primers/markers_used.csv')

    primers['chromosome'] = pd.Categorical(primers.chromosome, categories=['chr_8', 'chr_10', 'chr_12', 'chr_18'])
    markers['chromosome'] = pd.Categorical(markers.chromosome, categories=['chr_8', 'chr_10', 'chr_12', 'chr_18'])

    coverage_stats = {
        chromosome: pd.read_csv(f'results/coverage_stats/{chromosome}.csv')
        for chromosome in set(markers.chromosome)
    }

    primers['stop'] = primers.start + primers.seq.str.len()
    primers['amplimers'] = primers['name'].apply(lambda s: 1 if 'W' in s else (2 if 'COM' in s else 3))

    groupby = primers.groupby('chromosome')
    primer_boundaries = pd.DataFrame({'l': groupby.start.min(), 'r': groupby.stop.max()})


    fontsize = 8
    matplotlib.rc(
        'font',
        **{
            'size': fontsize,
            'family': 'sans-serif',
            'sans-serif': ['Helvetica']
        }
    )

    pgs_color = plt.get_cmap('Pastel2')(5 / 8)
    pgs_ec = plt.get_cmap('Dark2')(5 / 8)
    pgt_color = plt.get_cmap('Pastel1')(2 / 9)
    pgt_ec = plt.get_cmap('Set1')(2 / 9)
    primer_color = 'w'
    primer_ec = plt.get_cmap('Dark2')(8 / 8)

    fig = plt.figure(figsize=(7, 2.5))
    gs = gridspec.GridSpec(3, 1, figure=fig)
    plot_gs = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[:2, :])
    legend_gs = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[2, :])
    plots_gs = gridspec.GridSpecFromSubplotSpec(7, 4, subplot_spec=plot_gs[:, :], hspace=0, wspace=.08)

    lw = .75

    for i, ((chromosome, chromosome_marker), (_, chromosome_primers)) in enumerate(zip(markers.groupby('chromosome'), primers.groupby('chromosome'))):
        
        # Show up to 20 bases either side of the Pgt amplimer
        l = primer_boundaries.loc[chromosome, 'l'] - 20
        r = primer_boundaries.loc[chromosome, 'r'] + 20

        amplimer_ax = fig.add_subplot(plots_gs[:2, i])
        amplimer_ax.set_yticks([])
        amplimer_ax.set_xticklabels([])
        amplimer_ax.tick_params(axis='both', which='both', length=0)
        coverage_ax = fig.add_subplot(plots_gs[2:, i])
        coverage_ax.set_ylim(0, 29)
        coverage_ax.set_xlabel(f'Position on Chr{chromosome.split("_")[-1]}')
        
        # Add a label to only the left-most marker axis
        if i == 0:
            amplimer_ax.set_ylabel('Amplimers', rotation=0, ha='right', va='center')
            coverage_ax.set_ylabel('Depth', rotation=0, ha='right')
        else:
            # The y axis is shared among all markers
            coverage_ax.set_yticklabels([])
            coverage_ax.tick_params(axis='y', which='both', length=0)

        amplimer_ax.set_title(f'M{i + 1}', fontsize=fontsize)# don't use a higher fontsize

        coverage_ax.plot(
            coverage_stats[chromosome].loc[l:r].pgt_min,
            lw=lw,
            color=pgt_ec,
        )
        coverage_ax.fill_between(
            x=list(range(l, r + 1)),
            y1=coverage_stats[chromosome].loc[l:r].pgs_min,
            y2=coverage_stats[chromosome].loc[l:r].pgs_max,
            lw=lw,
            fc=pgs_color,
            ec=pgs_ec,
        )

        both_primers = chromosome_primers[chromosome_primers.amplifies == 'Pgt & Pgs'].sort_values(
            by='start').reset_index(drop=True)
        left_both_primer = both_primers.loc[0]
        right_both_primer = both_primers.loc[1]
        pgt_primer = chromosome_primers[chromosome_primers.amplifies == 'Pgt'].reset_index(
            drop=True).loc[0]
        x = left_both_primer.start
        width = (both_primers.loc[1].start + len(right_both_primer.seq)) - x
        
        
        amplimer_height = .5
        amplimer_props = {
            'height': amplimer_height,
            'lw': lw,
            'zorder': 2, # zorder 2 ensures amplimers behind primers
            'ec': primer_ec
        }
        amplimer_ax.add_patch(Rectangle(xy=(x, 1 + amplimer_height / 2 ), width=width, color=pgt_ec, **amplimer_props))
        amplimer_ax.add_patch(Rectangle(xy=(x, 2 + amplimer_height / 2), width=width, color=pgs_ec, **amplimer_props))
        if pgt_primer.start < x:
            width = (both_primers.loc[1].start + len(right_both_primer.seq)) - pgt_primer.start
            amplimer_ax.add_patch(Rectangle(xy=(pgt_primer.start, amplimer_height / 2), width=width, color=pgt_ec, **amplimer_props))
        else:
            width = (pgt_primer.start + len(pgt_primer.seq)) - left_both_primer.start
            amplimer_ax.add_patch(Rectangle(xy=(left_both_primer.start, amplimer_height / 2), width=width, color=pgt_ec, **amplimer_props))
        
        # Draw the primers
        n_amplimers = 1
        for y, height in zip([.1, 1.1, .1], [.8, 1.8, 2.8]):
            primer = chromosome_primers[chromosome_primers.amplimers == n_amplimers].iloc[0]
            amplimer_ax.add_patch(
                Rectangle(
                    xy=(primer.start, y),
                    width=len(primer.seq),
                    height=height,
                    color=primer_color,
                    ec=primer_ec, lw=lw, zorder=2,
                )
            )
            n_amplimers += 1

        amplimer_ax.set_ylim(-.2, 3.2)

        # avoid having to write an offset or e because that'd be too crowded
        coverage_ax.ticklabel_format(style='plain', useOffset=False, scilimits=(0, 0), axis='x',)

        for ax in [amplimer_ax, coverage_ax]:
            ax.set_xlim(l, r)
            ax.grid(alpha=.4, lw=lw)
            for side in ['top', 'left', 'bottom', 'right']:
                amplimer_ax.spines[side].set_linewidth(.75)

            # 300 is the only whole number which avoids collisions with axis numbers
            ax.xaxis.set_major_locator(ticker.MultipleLocator(300))
    

    # Create the legend
    legend_ax = fig.add_subplot(legend_gs[:, :])
    legend_ax.set_yticks([])
    legend_ax.set_xticks([])
    legend_ax.scatter([], [], s=0, label='Amplimers:')
    # Use primer_poly_path for the amplimer in the legend so we can get the height correcct 
    legend_ax.scatter([], [], color=pgs_ec, ec=primer_ec, lw=lw, marker=primer_poly_path('F', -1, 1, -.05, .4, point_length=0), label='$\it{Pgs}$ Amplimer', s=350)
    legend_ax.scatter([], [], color=pgt_ec, ec=primer_ec, lw=lw, marker=primer_poly_path('F', -1, 1, -.05, .4, point_length=0), label='$\it{Pgt}$ Amplimer', s=350)
    legend_ax.scatter([], [], s=0, label='Depth:')
    legend_ax.fill_between([], [], [], fc=pgs_color, ec=pgs_ec, label='$\it{Pgs}$ Depth Range', lw=lw)
    legend_ax.plot([], [], color=pgt_ec, label='$\it{Pgt}$ Minimum Depth', lw=lw)
    legend_ax.scatter([], [], color=primer_color, ec=primer_ec, lw=lw, marker=primer_poly_path('F', 0, 1, -.8, 2, point_length=0), label='Primer', s=55)
    for side in ['top', 'left', 'bottom', 'right']:
        legend_ax.spines[side].set_visible(False)
    legend_ax.legend(ncol=4)
    # Put the legend handles in the correct order
    handles, labels = legend_ax.get_legend_handles_labels()
    order = [1, 4, 2, 5, 3, 0, 6]
    legend = legend_ax.legend(
        [handles[idx] for idx in order],[labels[idx] for idx in order],
        ncol=4,
        loc='upper center',
        # Centering the legend explicitly because of spacing around text
        bbox_to_anchor=(.47, 1),
        framealpha=0, # Don't draw a box around the legend
    )
    # Make the right side of the "Depth:" label line up with the "Amplimers:" label
    renderer = fig.canvas.get_renderer()
    amplimer_right = [t.get_window_extent(renderer).x1 for t in legend.get_texts() if t.get_text() == 'Amplimers:'][0]
    depth_right = [t.get_window_extent(renderer).x1 for t in legend.get_texts() if t.get_text() == 'Depth:'][0]
    for t in legend.get_texts():
        if t.get_text() == 'Depth:':
            t.set_position((t.get_position()[0] + (amplimer_right - depth_right), t.get_position()[1]))

    # Because the primer marker is much narrower than the other symbols, there's a large gap
    # to the text - get rid of that gap
    for t in legend.get_texts():
        if t.get_text() == 'Primer':
            t.set_position((t.get_position()[0] - 6, t.get_position()[1]))

    fig.tight_layout()
    plt.savefig('figures/markers.svg', bbox_inches='tight', pad_inches=0)


