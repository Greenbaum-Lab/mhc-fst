'''
Manhattan plots of where every gene falls inside its own size matched null.

The height of a gene is its percentile in that null, so a gene near one has
changed more between the oldest and the newest time bin than random variant
sets of the same size, and a gene near zero has changed less. Genes called
outside the null are coloured by the side they leave it on.
'''

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from temporal_scan.scan_results import AUTOSOMES, HIGH, LOW

CHROMOSOME_COLORS = ['#b4b4b4', '#8a8a8a']
DIRECTION_COLORS = {HIGH: '#c0392b', LOW: '#2980b9'}
GROUP_MARKERS = ['o', '^', 's', 'D']
BACKGROUND_SIZE = 4
CALLED_SIZE = 16
FIGURE_SIZE = (14, 4.5)
PERCENTILE_LIMITS = (-0.03, 1.03)


def draw_background(axes, genes):
	for index, chromosome in enumerate(AUTOSOMES):
		block = genes[genes['chrom'] == chromosome]
		axes.scatter(block['position'], block['null_percentile'], s=BACKGROUND_SIZE, color=CHROMOSOME_COLORS[index % len(CHROMOSOME_COLORS)], linewidths=0)


def draw_called(axes, genes, marker):
	for direction, color in DIRECTION_COLORS.items():
		block = genes[(genes['outside_null'] == 1) & (genes['direction'] == direction)]
		axes.scatter(block['position'], block['null_percentile'], s=CALLED_SIZE, color=color, marker=marker, linewidths=0, label=f'{direction} ({len(block)})')


def label_axes(axes, centres, title, handles=None):
	axes.set_xticks(list(centres.values()))
	axes.set_xticklabels(list(centres.keys()), fontsize=7)
	axes.set_xlabel('chromosome')
	axes.set_ylabel('percentile in the size matched null')
	axes.set_ylim(*PERCENTILE_LIMITS)
	axes.set_title(title)
	axes.legend(handles=handles, fontsize=8, frameon=False, loc='center left', bbox_to_anchor=(1.0, 0.5))


def marker_handle(label, marker, color):
	return Line2D([], [], linestyle='none', marker=marker, color=color, label=label)


def save(figure, path):
	figure.tight_layout()
	figure.savefig(path, dpi=200)
	plt.close(figure)


def plot_group(genes, centres, group, path):
	'''
	Every usable gene of one group, with the genes called outside their null
	standing out by the side of the null they fall on.
	'''
	block = genes[genes['group'] == group]
	figure, axes = plt.subplots(figsize=FIGURE_SIZE)
	draw_background(axes, block)
	draw_called(axes, block, GROUP_MARKERS[0])
	label_axes(axes, centres, f'{group}: oldest against newest time bin, {len(block)} genes')
	save(figure, path)


def plot_shared(shared, centres, path):
	'''
	Only the genes called outside their null in every group, one marker per
	group, so that a gene answering the same way in both can be seen at once.
	'''
	figure, axes = plt.subplots(figsize=FIGURE_SIZE)
	handles = [marker_handle(direction, GROUP_MARKERS[0], color) for direction, color in DIRECTION_COLORS.items()]
	for index, group in enumerate(sorted(shared['group'].unique())):
		block = shared[shared['group'] == group]
		marker = GROUP_MARKERS[index % len(GROUP_MARKERS)]
		colors = [DIRECTION_COLORS[direction] for direction in block['direction']]
		axes.scatter(block['position'], block['null_percentile'], s=CALLED_SIZE * 2, marker=marker, color=colors, linewidths=0)
		handles.append(marker_handle(f'{group} ({len(block)})', marker, '#444444'))
	label_axes(axes, centres, f'genes called outside the null in every group, {shared["gene"].nunique()} genes', handles)
	save(figure, path)
