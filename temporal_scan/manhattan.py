'''
Manhattan plots of how far every gene falls outside its own size matched null.

The height of a gene is the two sided tail of its percentile in that null, so
ordinary genes sit near zero and only the genes departing from the null rise.
The score is signed: a gene that changed more between the oldest and the newest
time bin than random variant sets of the same size rises above the line, and a
gene that changed less falls below it.
'''

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from temporal_scan.scan_results import AUTOSOMES, HIGH, LOW

CHROMOSOME_COLORS = ['#b4b4b4', '#8a8a8a']
DIRECTION_COLORS = {HIGH: '#c0392b', LOW: '#2980b9'}
GROUP_MARKERS = ['o', '^', 's', 'D']
BACKGROUND_SIZE = 4
CALLED_SIZE = 16
FIGURE_SIZE = (14, 5.0)
LIMIT_MARGIN = 1.05


def draw_background(axes, genes):
	for index, chromosome in enumerate(AUTOSOMES):
		block = genes[genes['chrom'] == chromosome]
		axes.scatter(block['position'], block['score'], s=BACKGROUND_SIZE, color=CHROMOSOME_COLORS[index % len(CHROMOSOME_COLORS)], linewidths=0)


def draw_called(axes, genes, marker):
	for direction, color in DIRECTION_COLORS.items():
		block = genes[(genes['outside_null'] == 1) & (genes['direction'] == direction)]
		axes.scatter(block['position'], block['score'], s=CALLED_SIZE, color=color, marker=marker, linewidths=0, label=f'{direction} ({len(block)})')


def label_axes(axes, centres, title, handles=None):
	axes.set_xticks(list(centres.values()))
	axes.set_xticklabels(list(centres.keys()), fontsize=7)
	axes.set_xlabel('chromosome')
	axes.set_ylabel('signed -log10 two sided tail against the null')
	axes.axhline(0.0, color='#666666', linewidth=0.8)
	limit = max(abs(value) for value in axes.get_ylim()) * LIMIT_MARGIN
	axes.set_ylim(-limit, limit)
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
		axes.scatter(block['position'], block['score'], s=CALLED_SIZE * 2, marker=marker, color=colors, linewidths=0)
		handles.append(marker_handle(f'{group} ({len(block)})', marker, '#444444'))
	label_axes(axes, centres, f'genes called outside the null in every group, {shared["gene"].nunique()} genes', handles)
	save(figure, path)
