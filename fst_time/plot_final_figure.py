'''
The finished figure: one panel per focal locus, in a column for the trend
expected of it, drawn from the results tables alone.

The loci drawn and their order are the list below, which is a choice of the
figure and not of the run, so a measured locus can be left out of the page
without measuring anything again. Both jackknives are drawn, one figure each.

	python plot_final_figure.py --results-dir .

It reads `fst_time_series.csv` and `gene_background.csv` and writes
`fst_final_individuals.png` and `fst_final_snp_blocks.png` beside them.

Every panel is placed by hand, in inches, so all of them are the same size
whatever they are labelled and however many loci a column holds.
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
PANEL_WIDTH_INCHES = 4.0
PANEL_HEIGHT_INCHES = 1.45
LEFT_MARGIN_INCHES = 0.75
RIGHT_MARGIN_INCHES = 0.3
TITLE_INCHES = 0.6
TICKS_INCHES = 0.35
HEADER_INCHES = 0.5
TIME_LABEL_INCHES = 0.4
LEGEND_INCHES = 0.4
TIME_LABEL = 'Thousand years before present'
FST_LABEL = '$F_{ST}$'


def bin_midpoints(rows):
	return (rows['time_start'] + rows['time_end']) / 2000.0


def column_pitch():
	return LEFT_MARGIN_INCHES + PANEL_WIDTH_INCHES + RIGHT_MARGIN_INCHES


def row_pitch():
	return TITLE_INCHES + PANEL_HEIGHT_INCHES + TICKS_INCHES


def figure_size(column_count, row_count):
	return (
		column_count * column_pitch(),
		HEADER_INCHES + row_count * row_pitch() + TIME_LABEL_INCHES + LEGEND_INCHES)


def panel_axes(figure, row, column, row_count):
	'''
	The box of one panel, given in inches and turned into figure fractions, so
	every panel of the page holds exactly the same area.
	'''
	width, height = figure.get_figwidth(), figure.get_figheight()
	left = column * column_pitch() + LEFT_MARGIN_INCHES
	bottom = LEGEND_INCHES + TIME_LABEL_INCHES + (row_count - 1 - row) * row_pitch() + TICKS_INCHES
	return figure.add_axes([
		left / width, bottom / height, PANEL_WIDTH_INCHES / width, PANEL_HEIGHT_INCHES / height])


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


def draw_panel(axis, locus, target_rows, background_rows, gene_background, uncertainty, limits):
	shade_periods(axis)
	draw_series(axis, background_rows, f'ci_low_{uncertainty}', f'ci_high_{uncertainty}', BACKGROUND_COLOR, '--')
	draw_series(axis, gene_background, 'ci_low', 'ci_high', GENE_BACKGROUND_COLOR, '--')
	draw_series(axis, target_rows, f'ci_low_{uncertainty}', f'ci_high_{uncertainty}', FOCAL_COLOR, '-')
	axis.set_title(f'{locus["label"]}\n{locus["phenotype"]}', fontsize=9, linespacing=1.3)
	axis.set_ylabel(FST_LABEL, fontsize=9)
	axis.set_ylim(limits)
	axis.tick_params(labelsize=8)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)


def time_limits(table):
	'''
	The axis runs from the oldest bin to the most recent one, ending where the
	data ends.
	'''
	years = bin_midpoints(table)
	return years.max(), years.min()


def column_limits(table, loci, background_rows, gene_background, uncertainty):
	'''
	One pair of limits for a whole column, so the loci expected to do the same
	thing are read on the same scale.
	'''
	frames = [table[table['target'] == locus['label']] for locus in loci]
	low, high = series_limits(frames + [background_rows], f'ci_low_{uncertainty}', f'ci_high_{uncertainty}')
	gene_low, gene_high = series_limits([gene_background], 'ci_low', 'ci_high')
	return min(low, gene_low), max(high, gene_high)


def fill_column(figure, column, row_count, loci, table, background_rows, gene_background, uncertainty):
	limits = column_limits(table, loci, background_rows, gene_background, uncertainty)
	for row, locus in enumerate(loci):
		axis = panel_axes(figure, row, column, row_count)
		target_rows = table[table['target'] == locus['label']].sort_values('time_start')
		draw_panel(axis, locus, target_rows, background_rows, gene_background, uncertainty, limits)
		axis.set_xlim(time_limits(table))


def add_column_labels(figure, column, trend, locus_count):
	'''
	The trend a column stands for above it and the time axis below it, once
	each, while the years themselves stay on every panel.
	'''
	height = figure.get_figheight()
	center = (column * column_pitch() + LEFT_MARGIN_INCHES + PANEL_WIDTH_INCHES / 2.0) / figure.get_figwidth()
	figure.text(center, 1.0 - HEADER_INCHES / height * 0.7, f'expected {trend} ({locus_count})',
	            ha='center', va='center', fontsize=13, fontweight='bold')
	figure.text(center, (LEGEND_INCHES + TIME_LABEL_INCHES * 0.45) / height, TIME_LABEL,
	            ha='center', va='center', fontsize=10)


def add_legend(figure):
	handles = [
		Patch(facecolor=PERIOD_COLORS[position % len(PERIOD_COLORS)], alpha=PERIOD_ALPHA, label=period['period'])
		for position, period in enumerate(PERIODS)
	]
	handles.append(Line2D([], [], color=FOCAL_COLOR, label='locus'))
	handles.append(Line2D([], [], color=BACKGROUND_COLOR, linestyle='--', label='genome wide'))
	handles.append(Line2D([], [], color=GENE_BACKGROUND_COLOR, linestyle='--', label='mean over all annotated genes'))
	figure.legend(handles=handles, loc='center', bbox_to_anchor=(0.5, LEGEND_INCHES / figure.get_figheight() * 0.5),
	              ncol=len(handles), frameon=False, fontsize=8)


def loci_by_trend():
	return {
		trend: [locus for locus in FOCAL_LOCI if locus['trend'] == trend]
		for trend in TREND_ORDER
	}


def build_figure(table, gene_background, uncertainty):
	'''
	One column per expected trend, one panel per locus, every panel of a
	column on the same axes.
	'''
	columns = loci_by_trend()
	row_count = max(len(loci) for loci in columns.values())
	background_rows = table[table['target'] == GENOME_WIDE_TARGET].sort_values('time_start')
	figure = plt.figure(figsize=figure_size(len(columns), row_count))
	for column, (trend, loci) in enumerate(columns.items()):
		fill_column(figure, column, row_count, loci, table, background_rows, gene_background, uncertainty)
		add_column_labels(figure, column, trend, len(loci))
	add_legend(figure)
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
