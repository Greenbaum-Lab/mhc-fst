'''
FST against the number of variants behind it, the check that the size matched
null behaves.

The null interval should be wide where a gene holds few variants and close in
as the gene holds more, and the genes called outside the null should be the
ones leaving that band. A band that fails to narrow, or calls that sit inside
it, means the normalisation is not doing what the Manhattan assumes.
'''

import matplotlib.pyplot as plt

from temporal_scan.manhattan import DIRECTION_COLORS, save

BAND_COLOR = '#333333'
GENOME_COLOR = '#27ae60'
INSIDE_COLOR = '#b4b4b4'
POINT_SIZE = 5
CALLED_SIZE = 14
FIGURE_SIZE = (7, 5)


def draw_null_band(axes, genes):
	ordered = genes.sort_values('n_variants')
	axes.plot(ordered['n_variants'], ordered['null_ci_low'], color=BAND_COLOR, linewidth=1)
	axes.plot(ordered['n_variants'], ordered['null_ci_high'], color=BAND_COLOR, linewidth=1, label='null interval')
	axes.axhline(ordered['genome_fst'].iloc[0], color=GENOME_COLOR, linewidth=1, linestyle='--', label='genome wide FST')


def draw_genes(axes, genes):
	inside = genes[genes['outside_null'] != 1]
	axes.scatter(inside['n_variants'], inside['fst'], s=POINT_SIZE, color=INSIDE_COLOR, linewidths=0, label=f'inside null ({len(inside)})')
	for direction, color in DIRECTION_COLORS.items():
		block = genes[(genes['outside_null'] == 1) & (genes['direction'] == direction)]
		axes.scatter(block['n_variants'], block['fst'], s=CALLED_SIZE, color=color, linewidths=0, label=f'{direction} ({len(block)})')


def plot_diagnostic(genes, group, path):
	'''
	The genes of one group against the number of variants each was measured
	from, with the null interval that number buys drawn over them.
	'''
	figure, axes = plt.subplots(figsize=FIGURE_SIZE)
	draw_genes(axes, genes)
	draw_null_band(axes, genes)
	axes.set_xscale('log')
	axes.set_xlabel('variants in the gene')
	axes.set_ylabel('FST')
	axes.set_title(f'{group}: gene FST against gene size')
	axes.legend(fontsize=8, frameon=False)
	save(figure, path)
