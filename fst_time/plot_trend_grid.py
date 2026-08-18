'''
Every locus on one page, in a column for the trend expected of it.

Reads only the results table, the gene background and the period list, so it
can sit next to them and run away from the pipeline and the cluster.

	python plot_trend_grid.py
	python plot_trend_grid.py --include LCT MHC ERAP2
	python plot_trend_grid.py --exclude G6PD EPAS1

Every input defaults to the working directory, and the figure is written there
too.
'''

import json
import argparse
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

GENOME_WIDE_TARGET = 'genome_wide'
TREND_COLUMNS = ['high', 'low', 'neutral']
UNCERTAINTIES = ['snp_blocks', 'individuals']
FOCAL_COLOR = '#c0392b'
GENE_BACKGROUND_COLOR = '#2e7d32'
BACKGROUND_COLOR = 'black'
PERIOD_COLORS = ['#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3']
PERIOD_ALPHA = 0.13
BAND_ALPHA = 0.2
HEADER_INCHES = 0.55
FOOTER_INCHES = 0.80
PANEL_INCHES = 1.7


def load_periods(periods_path):
	with open(periods_path) as periods_file:
		return json.load(periods_file)


def bin_midpoints(rows):
	return (rows['time_start'] + rows['time_end']) / 2000.0


def shade_periods(axis, periods):
	'''
	Each period as a band across the panel, so a locus is read against the
	archaeological sequence rather than against bare years.
	'''
	for position, period in enumerate(periods):
		axis.axvspan(
			period['end_bp'] / 1000.0, period['start_bp'] / 1000.0,
			color=PERIOD_COLORS[position % len(PERIOD_COLORS)], alpha=PERIOD_ALPHA, linewidth=0)


def draw_series(axis, rows, uncertainty, color, line_style):
	'''
	One series as a line with its jackknife interval as a shaded band.
	'''
	years = bin_midpoints(rows)
	axis.plot(years, rows['fst'], color=color, linestyle=line_style, linewidth=1.4)
	axis.fill_between(
		years, rows[f'ci_low_{uncertainty}'], rows[f'ci_high_{uncertainty}'],
		color=color, alpha=BAND_ALPHA, linewidth=0)


def draw_gene_background(axis, gene_background):
	'''
	The mean over every gene of the annotation, the same in every panel, with
	error bars spanning the genes.
	'''
	axis.errorbar(
		bin_midpoints(gene_background), gene_background['fst'],
		yerr=[gene_background['fst'] - gene_background['ci_low'],
		      gene_background['ci_high'] - gene_background['fst']],
		color=GENE_BACKGROUND_COLOR, linewidth=1.2, marker='s', markersize=3,
		capsize=2.5, elinewidth=0.9)


def draw_expected_time(axis, time_bp):
	'''
	The time the expectation changes, for the loci that have one.
	'''
	if pd.isna(time_bp):
		return
	axis.axvline(time_bp / 1000.0, color='0.15', linestyle=':', linewidth=1.3)


def draw_panel(axis, target_rows, background_rows, gene_background, periods, uncertainty):
	shade_periods(axis, periods)
	draw_series(axis, background_rows, uncertainty, BACKGROUND_COLOR, '--')
	draw_gene_background(axis, gene_background)
	draw_series(axis, target_rows, uncertainty, FOCAL_COLOR, '-')
	draw_expected_time(axis, target_rows.iloc[0]['time_bp'])
	if target_rows['n_variants'].max() == 0:
		axis.annotate('no variants', xy=(0.5, 0.5), xycoords='axes fraction',
		              ha='center', fontsize=9, color='0.4')


def style_panel(axis, target_rows):
	first = target_rows.iloc[0]
	axis.set_title(f'{first["target"]}\n{first["phenotype"]}', fontsize=9, linespacing=1.3)
	axis.tick_params(labelsize=8)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)


def chosen_targets(table, included, excluded):
	'''
	The loci to draw, kept to those named where any are named, and without the
	ones left out.
	'''
	focal = table[table['target'] != GENOME_WIDE_TARGET]
	if included:
		focal = focal[focal['target'].isin(included)]
	return focal[~focal['target'].isin(excluded)]


def trend_targets(table, included, excluded):
	'''
	The loci of each trend column, in the order the table holds them.
	'''
	focal = chosen_targets(table, included, excluded)
	return {
		trend: list(dict.fromkeys(focal[focal['trend'] == trend]['target']))
		for trend in TREND_COLUMNS
	}


def fill_column(column_axes, targets, table, background_rows, gene_background, periods, uncertainty):
	'''
	One panel per locus, with the years read off the lowest panel the column
	fills rather than the lowest panel of the grid.
	'''
	for axis, target in zip(column_axes, targets):
		target_rows = table[table['target'] == target].sort_values('time_start')
		draw_panel(axis, target_rows, background_rows, gene_background, periods, uncertainty)
		style_panel(axis, target_rows)
	for axis in column_axes[len(targets):]:
		axis.set_visible(False)
	if targets:
		column_axes[len(targets) - 1].tick_params(labelbottom=True)


def add_column_headers(figure, axes, targets_by_trend, top):
	for column_position, trend in enumerate(TREND_COLUMNS):
		box = axes[0, column_position].get_position()
		figure.text(
			box.x0 + box.width / 2.0, top + (1.0 - top) * 0.35,
			f'expected {trend} ({len(targets_by_trend[trend])})',
			ha='center', fontsize=13, fontweight='bold')


def add_legend(figure, periods):
	handles = [
		Patch(facecolor=PERIOD_COLORS[position % len(PERIOD_COLORS)], alpha=PERIOD_ALPHA,
		      label=period['period'].replace('_', ' '))
		for position, period in enumerate(periods)
	]
	handles.append(Line2D([], [], color=FOCAL_COLOR, label='locus'))
	handles.append(Line2D([], [], color=BACKGROUND_COLOR, linestyle='--', label='genome wide'))
	handles.append(Line2D([], [], color=GENE_BACKGROUND_COLOR, marker='s', markersize=3,
	                      label='mean over all annotated genes'))
	figure.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.0),
	              ncol=len(handles), frameon=False, fontsize=8)


def finish_layout(figure, axes, targets_by_trend, periods):
	'''
	Reserve the height the headers and the legend were given when the figure
	was sized, so a grid of any number of rows keeps the same margins and the
	panels keep the same height.
	'''
	height = figure.get_figheight()
	bottom, top = FOOTER_INCHES / height, 1.0 - HEADER_INCHES / height
	figure.tight_layout(rect=[0.03, bottom, 1, top])
	add_column_headers(figure, axes, targets_by_trend, top)
	add_legend(figure, periods)
	figure.supxlabel('Thousand years before present', fontsize=11, y=bottom * 0.42)
	figure.supylabel('$F_{ST}$', fontsize=12)


def build_figure(table, gene_background, periods, uncertainty, included, excluded):
	'''
	One column per expected trend, one panel per locus, all on a shared time
	axis running from oldest to most recent.
	'''
	targets_by_trend = trend_targets(table, included, excluded)
	row_count = max(len(targets) for targets in targets_by_trend.values())
	figure, axes = plt.subplots(
		row_count, len(TREND_COLUMNS), squeeze=False, sharex=True,
		figsize=(4.2 * len(TREND_COLUMNS), PANEL_INCHES * row_count + HEADER_INCHES + FOOTER_INCHES))
	background_rows = table[table['target'] == GENOME_WIDE_TARGET].sort_values('time_start')
	for column_position, trend in enumerate(TREND_COLUMNS):
		fill_column(axes[:, column_position], targets_by_trend[trend], table,
		            background_rows, gene_background, periods, uncertainty)
	axes[0, 0].set_xlim(table['time_end'].max() / 1000.0, 0)
	finish_layout(figure, axes, targets_by_trend, periods)
	return figure


def main():
	parser = argparse.ArgumentParser(description='All loci on one page, in a column per expected trend')
	parser.add_argument('--results', default='fst_time_series.csv')
	parser.add_argument('--periods', default='time_periods.json')
	parser.add_argument('--gene-background', default='gene_background.csv')
	parser.add_argument('--output-dir', default='.')
	parser.add_argument('--uncertainty', choices=UNCERTAINTIES, default=UNCERTAINTIES[0])
	parser.add_argument('--include', nargs='*', default=[], metavar='LOCUS')
	parser.add_argument('--exclude', nargs='*', default=[], metavar='LOCUS')
	args = parser.parse_args()
	table = pd.read_csv(args.results)
	gene_background = pd.read_csv(args.gene_background).sort_values('time_start')
	figure = build_figure(
		table, gene_background, load_periods(args.periods), args.uncertainty,
		args.include, args.exclude)
	figure.savefig(pathlib.Path(args.output_dir) / f'fst_trend_grid_{args.uncertainty}.png', dpi=200)


if __name__ == '__main__':
	main()
