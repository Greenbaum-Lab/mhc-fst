'''
The finished figure: one panel per focal locus, in a column for the trend
expected of it, drawn from the results tables alone.

	python plot_final_figure.py --results-dir .

It reads `fst_time_series.csv` and `gene_background.csv` and writes
`fst_final_individuals.png` and `fst_final_snp_blocks.png` beside them.
'''

import argparse
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

FOCAL_LOCI = [
	{'label': 'SLC24A5', 'phenotype': 'pigmentation', 'trend': 'high'},
	{'label': 'SLC45A2', 'phenotype': 'pigmentation', 'trend': 'high'},
	{'label': 'HERC2/OCA2', 'phenotype': 'pigmentation', 'trend': 'high'},
	{'label': 'EDAR', 'phenotype': 'ectodermal morphology', 'trend': 'high'},
	{'label': 'FADS1_FADS2', 'phenotype': 'lipid metabolism', 'trend': 'low'},
	{'label': 'TLR6_TLR1_TLR10', 'phenotype': 'innate immunity', 'trend': 'low'},
	{'label': 'MARK3', 'phenotype': 'infectious disease immunity', 'trend': 'low'},
	{'label': 'IL23R', 'phenotype': 'immune signalling', 'trend': 'low'},
	{'label': 'IL1RL1', 'phenotype': 'immune signalling', 'trend': 'low'},
]
PERIODS = [
	{'period': 'Neolithic Chalcolithic', 'start_bp': 8500, 'end_bp': 4500},
	{'period': 'Bronze Age', 'start_bp': 4500, 'end_bp': 3150},
	{'period': 'Iron Age Historical', 'start_bp': 3150, 'end_bp': 500},
]
TREND_ORDER = ['high', 'low']
UNCERTAINTIES = ['individuals', 'snp_blocks']
GENOME_WIDE_TARGET = 'genome_wide'
FOCAL_COLOR = '#c0392b'
GENE_BACKGROUND_COLOR = '#2e7d32'
BACKGROUND_COLOR = 'black'
PERIOD_COLORS = ['#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3']
PERIOD_ALPHA = 0.13
BAND_ALPHA = 0.2
AXIS_MARGIN = 0.05
PANEL_WIDTH = 5.0
PANEL_HEIGHT = 2.2
TIME_LABEL = 'Thousand years before present'
FST_LABEL = '$F_{ST}$'


def bin_midpoints(rows):
	return (rows['time_start'] + rows['time_end']) / 2000.0


def shade_periods(axis):
	'''
	Each period as a band across the panel, so a locus is read against the
	archaeological sequence rather than against bare years.
	'''
	for position, period in enumerate(PERIODS):
		axis.axvspan(
			period['end_bp'] / 1000.0, period['start_bp'] / 1000.0,
			color=PERIOD_COLORS[position % len(PERIOD_COLORS)], alpha=PERIOD_ALPHA, linewidth=0)


def draw_series(axis, rows, low_column, high_column, color, line_style):
	'''
	One series as a line with its interval as a shaded band.
	'''
	years = bin_midpoints(rows)
	axis.plot(years, rows['fst'], color=color, linestyle=line_style, linewidth=1.4)
	axis.fill_between(years, rows[low_column], rows[high_column], color=color, alpha=BAND_ALPHA, linewidth=0)


def series_limits(frames, low_column, high_column):
	'''
	The lowest and highest value any of these series reaches, bands included.
	'''
	low = min(frame[low_column].min() for frame in frames)
	high = max(frame[high_column].max() for frame in frames)
	margin = (high - low) * AXIS_MARGIN
	return low - margin, high + margin


def column_limits(table, loci, background_rows, gene_background, uncertainty):
	'''
	One pair of limits for a whole column, so the loci expected to do the same
	thing are read on the same scale.
	'''
	frames = [table[table['target'] == locus['label']] for locus in loci]
	low, high = series_limits(frames + [background_rows], f'ci_low_{uncertainty}', f'ci_high_{uncertainty}')
	gene_low, gene_high = series_limits([gene_background], 'ci_low', 'ci_high')
	return min(low, gene_low), max(high, gene_high)


def time_limits(table):
	years = bin_midpoints(table)
	return years.max(), years.min()


def draw_panel(axis, locus, target_rows, background_rows, gene_background, uncertainty):
	shade_periods(axis)
	draw_series(axis, background_rows, f'ci_low_{uncertainty}', f'ci_high_{uncertainty}', BACKGROUND_COLOR, '--')
	draw_series(axis, gene_background, 'ci_low', 'ci_high', GENE_BACKGROUND_COLOR, '--')
	draw_series(axis, target_rows, f'ci_low_{uncertainty}', f'ci_high_{uncertainty}', FOCAL_COLOR, '-')
	axis.set_title(f'{locus["label"]}\n{locus["phenotype"]}', fontsize=10, linespacing=1.3)
	axis.set_ylabel(FST_LABEL, fontsize=10)
	axis.tick_params(labelsize=9)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)


def fill_column(column_axes, loci, table, background_rows, gene_background, uncertainty):
	limits = column_limits(table, loci, background_rows, gene_background, uncertainty)
	for axis, locus in zip(column_axes, loci):
		target_rows = table[table['target'] == locus['label']].sort_values('time_start')
		draw_panel(axis, locus, target_rows, background_rows, gene_background, uncertainty)
		axis.set_xlim(time_limits(table))
		axis.set_ylim(limits)
	for axis in column_axes[len(loci):]:
		axis.axis('off')
	column_axes[len(loci) - 1].set_xlabel(TIME_LABEL, fontsize=10)


def add_headers(figure, axes, columns):
	for column, (trend, loci) in enumerate(columns.items()):
		box = axes[0, column].get_position()
		figure.text(box.x0 + box.width / 2.0, 0.99, f'expected {trend} ({len(loci)})',
		            ha='center', va='top', fontsize=14, fontweight='bold')


def add_legend(figure):
	handles = [
		Patch(facecolor=PERIOD_COLORS[position % len(PERIOD_COLORS)], alpha=PERIOD_ALPHA, label=period['period'])
		for position, period in enumerate(PERIODS)
	]
	handles.append(Line2D([], [], color=FOCAL_COLOR, label='locus'))
	handles.append(Line2D([], [], color=BACKGROUND_COLOR, linestyle='--', label='genome wide'))
	handles.append(Line2D([], [], color=GENE_BACKGROUND_COLOR, linestyle='--', label='mean over all annotated genes'))
	figure.legend(handles=handles, loc='lower center', ncol=len(handles), frameon=False, fontsize=9)


def loci_by_trend():
	return {
		trend: [locus for locus in FOCAL_LOCI if locus['trend'] == trend]
		for trend in TREND_ORDER
	}


def build_figure(table, gene_background, uncertainty):
	'''
	One column per expected trend, one panel per locus, all on one grid so
	every panel is the same size.
	'''
	columns = loci_by_trend()
	row_count = max(len(loci) for loci in columns.values())
	background_rows = table[table['target'] == GENOME_WIDE_TARGET].sort_values('time_start')
	figure, axes = plt.subplots(
		row_count, len(columns), squeeze=False,
		figsize=(PANEL_WIDTH * len(columns), PANEL_HEIGHT * row_count))
	for column, loci in enumerate(columns.values()):
		fill_column(axes[:, column], loci, table, background_rows, gene_background, uncertainty)
	add_legend(figure)
	figure.tight_layout(rect=[0, 0.03, 1, 0.97])
	add_headers(figure, axes, columns)
	return figure


def main():
	parser = argparse.ArgumentParser(description='The finished trend grid, one figure per jackknife')
	parser.add_argument('--results-dir', default='.')
	parser.add_argument('--output-dir', default='.')
	parser.add_argument('--dpi', type=int, default=300)
	args = parser.parse_args()
	results_dir = pathlib.Path(args.results_dir)
	table = pd.read_csv(results_dir / 'fst_time_series.csv')
	gene_background = pd.read_csv(results_dir / 'gene_background.csv').sort_values('time_start')
	for uncertainty in UNCERTAINTIES:
		figure = build_figure(table, gene_background, uncertainty)
		figure.savefig(pathlib.Path(args.output_dir) / f'fst_final_{uncertainty}.png', dpi=args.dpi)
	plt.show()


if __name__ == '__main__':
	main()
