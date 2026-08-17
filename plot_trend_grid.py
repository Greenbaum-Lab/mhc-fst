'''
Every locus on one page, in a column for the trend expected of it.

Reads only the results table and the period list, so it can sit next to them
and run away from the pipeline and the cluster.

	python plot_trend_grid.py
	python plot_trend_grid.py --exclude G6PD EPAS1

Both inputs default to the working directory, and the figure is written there
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
from statistics import NormalDist

GENOME_WIDE_TARGET = 'genome_wide'
TREND_COLUMNS = ['high', 'low', 'neutral']
UNCERTAINTIES = ['snp_blocks', 'individuals']
FOCAL_COLOR = '#c0392b'
MEAN_COLOR = '#1f4e79'
GENE_BACKGROUND_COLOR = '#2e7d32'
BACKGROUND_COLOR = 'black'
PERIOD_COLORS = ['#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3']
PERIOD_ALPHA = 0.13
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


def draw_gene_background(axis, gene_background):
	'''
	The mean over every gene of the annotation, with error bars spanning the
	genes, which is the background a focal locus is compared against.
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
	axis.plot(bin_midpoints(background_rows), background_rows['fst'],
	          color=BACKGROUND_COLOR, linestyle='--', linewidth=1.0)
	draw_gene_background(axis, gene_background)
	axis.plot(bin_midpoints(target_rows), target_rows['fst'], color=FOCAL_COLOR, linewidth=1.5)
	axis.fill_between(
		bin_midpoints(target_rows), target_rows[f'ci_low_{uncertainty}'], target_rows[f'ci_high_{uncertainty}'],
		color=FOCAL_COLOR, alpha=0.2, linewidth=0)
	draw_expected_time(axis, target_rows.iloc[0]['time_bp'])
	if target_rows['n_variants'].max() == 0:
		axis.annotate('no variants', xy=(0.5, 0.5), xycoords='axes fraction',
		              ha='center', fontsize=9, color='0.4')


def mean_across_targets(table, targets, confidence_level):
	'''
	The mean FST of a set of loci at each time bin, and how widely they spread
	around it. The interval is the spread between loci, not the uncertainty of
	any single one of them, so it answers whether the loci agree.
	'''
	rows = table[table['target'].isin(targets)]
	summary = rows.groupby(['time_start', 'time_end'])['fst'].agg(['mean', 'std', 'count']).reset_index()
	standard_error = (summary['std'] / np.sqrt(summary['count'])).fillna(0.0)
	quantile = NormalDist().inv_cdf(1.0 - (1.0 - confidence_level) / 2.0)
	summary['fst'] = summary['mean']
	summary['ci_low'] = summary['mean'] - quantile * standard_error
	summary['ci_high'] = summary['mean'] + quantile * standard_error
	return summary.sort_values('time_start')


def draw_mean_panel(axis, summary, background_rows, gene_background, periods):
	'''
	The loci of one column averaged, with error bars across them.
	'''
	shade_periods(axis, periods)
	axis.plot(bin_midpoints(background_rows), background_rows['fst'],
	          color=BACKGROUND_COLOR, linestyle='--', linewidth=1.0)
	draw_gene_background(axis, gene_background)
	axis.errorbar(
		bin_midpoints(summary), summary['fst'],
		yerr=[summary['fst'] - summary['ci_low'], summary['ci_high'] - summary['fst']],
		color=MEAN_COLOR, linewidth=1.6, marker='o', markersize=3.5, capsize=3, elinewidth=1.0)


def style_mean_panel(axis, summary):
	locus_count = int(summary['count'].max())
	axis.set_title(
		f'mean of {locus_count} loc{"us" if locus_count == 1 else "i"}\nbars span loci',
		fontsize=9, linespacing=1.3)
	axis.tick_params(labelsize=8)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)


def style_panel(axis, target_rows):
	first = target_rows.iloc[0]
	axis.set_title(f'{first["target"]}\n{first["phenotype"]}', fontsize=9, linespacing=1.3)
	axis.tick_params(labelsize=8)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)


def trend_targets(table, excluded):
	'''
	The loci of each trend column, in the order the table holds them, without
	the ones left out.
	'''
	focal = table[(table['target'] != GENOME_WIDE_TARGET) & ~table['target'].isin(excluded)]
	return {
		trend: list(dict.fromkeys(focal[focal['trend'] == trend]['target']))
		for trend in TREND_COLUMNS
	}


def fill_column(column_axes, targets, table, background_rows, gene_background, periods, uncertainty, confidence_level):
	'''
	The loci of one column averaged in the top panel and drawn one per panel
	below it, with the years read off the lowest panel the column fills rather
	than the lowest panel of the grid.
	'''
	if targets:
		summary = mean_across_targets(table, targets, confidence_level)
		draw_mean_panel(column_axes[0], summary, background_rows, gene_background, periods)
		style_mean_panel(column_axes[0], summary)
	for axis, target in zip(column_axes[1:], targets):
		target_rows = table[table['target'] == target].sort_values('time_start')
		draw_panel(axis, target_rows, background_rows, gene_background, periods, uncertainty)
		style_panel(axis, target_rows)
	for axis in column_axes[len(targets) + 1:]:
		axis.set_visible(False)
	if targets:
		column_axes[len(targets)].tick_params(labelbottom=True)


def add_column_headers(figure, axes, targets_by_trend, top):
	for column_position, trend in enumerate(TREND_COLUMNS):
		box = axes[0, column_position].get_position()
		figure.text(
			box.x0 + box.width / 2.0, top + (1.0 - top) * 0.35,
			f'expected {trend} ({len(targets_by_trend[trend])})',
			ha='center', fontsize=13, fontweight='bold')


def add_period_legend(figure, periods):
	handles = [
		Patch(facecolor=PERIOD_COLORS[position % len(PERIOD_COLORS)], alpha=PERIOD_ALPHA,
		      label=period['period'].replace('_', ' '))
		for position, period in enumerate(periods)
	]
	handles.append(Line2D([], [], color=BACKGROUND_COLOR, linestyle='--', label='genome wide'))
	handles.append(Line2D([], [], color=GENE_BACKGROUND_COLOR, marker='s', markersize=3,
	                      label='mean over all annotated genes'))
	figure.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.0),
	              ncol=len(periods), frameon=False, fontsize=9)


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
	add_period_legend(figure, periods)
	figure.supxlabel('Thousand years before present', fontsize=11, y=bottom * 0.42)
	figure.supylabel('$F_{ST}$', fontsize=12)


def build_figure(table, gene_background, periods, uncertainty, excluded, confidence_level):
	'''
	One column per expected trend, the loci of that column averaged in its top
	panel and drawn one per panel below, all on a shared time axis running
	from oldest to most recent.
	'''
	targets_by_trend = trend_targets(table, excluded)
	row_count = 1 + max(len(targets) for targets in targets_by_trend.values())
	figure, axes = plt.subplots(
		row_count, len(TREND_COLUMNS), squeeze=False, sharex=True,
		figsize=(4.2 * len(TREND_COLUMNS), PANEL_INCHES * row_count + HEADER_INCHES + FOOTER_INCHES))
	background_rows = table[table['target'] == GENOME_WIDE_TARGET].sort_values('time_start')
	for column_position, trend in enumerate(TREND_COLUMNS):
		fill_column(axes[:, column_position], targets_by_trend[trend], table,
		            background_rows, gene_background, periods, uncertainty, confidence_level)
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
	parser.add_argument('--exclude', nargs='*', default=[], metavar='LOCUS')
	parser.add_argument('--confidence-level', type=float, default=0.95)
	args = parser.parse_args()
	table = pd.read_csv(args.results)
	gene_background = pd.read_csv(args.gene_background).sort_values('time_start')
	figure = build_figure(
		table, gene_background, load_periods(args.periods), args.uncertainty,
		args.exclude, args.confidence_level)
	figure.savefig(pathlib.Path(args.output_dir) / f'fst_trend_grid_{args.uncertainty}.png', dpi=200)


if __name__ == '__main__':
	main()
