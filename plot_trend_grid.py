'''
Every locus on one page, in a column for the trend expected of it.

Reads only the results table and the period list, so it can sit next to them
and run away from the pipeline and the cluster.

	python plot_trend_grid.py

Both inputs default to the working directory, and the figure is written there
too.
'''

import json
import argparse
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

GENOME_WIDE_TARGET = 'genome_wide'
TREND_COLUMNS = ['high', 'low', 'neutral']
UNCERTAINTIES = ['snp_blocks', 'individuals']
FOCAL_COLOR = '#c0392b'
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


def draw_expected_time(axis, time_bp):
	'''
	The time the expectation changes, for the loci that have one.
	'''
	if pd.isna(time_bp):
		return
	axis.axvline(time_bp / 1000.0, color='0.15', linestyle=':', linewidth=1.3)


def draw_panel(axis, target_rows, background_rows, periods, uncertainty):
	shade_periods(axis, periods)
	axis.plot(bin_midpoints(background_rows), background_rows['fst'],
	          color=BACKGROUND_COLOR, linestyle='--', linewidth=1.0)
	axis.plot(bin_midpoints(target_rows), target_rows['fst'], color=FOCAL_COLOR, linewidth=1.5)
	axis.fill_between(
		bin_midpoints(target_rows), target_rows[f'ci_low_{uncertainty}'], target_rows[f'ci_high_{uncertainty}'],
		color=FOCAL_COLOR, alpha=0.2, linewidth=0)
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


def trend_targets(table):
	'''
	The loci of each trend column, in the order the table holds them.
	'''
	focal = table[table['target'] != GENOME_WIDE_TARGET]
	return {
		trend: list(dict.fromkeys(focal[focal['trend'] == trend]['target']))
		for trend in TREND_COLUMNS
	}


def fill_column(column_axes, targets, table, background_rows, periods, uncertainty):
	'''
	One panel per locus, and the years read off the lowest panel the column
	fills rather than the lowest panel of the grid.
	'''
	for axis, target in zip(column_axes, targets):
		target_rows = table[table['target'] == target].sort_values('time_start')
		draw_panel(axis, target_rows, background_rows, periods, uncertainty)
		style_panel(axis, target_rows)
	for axis in column_axes[len(targets):]:
		axis.set_visible(False)
	column_axes[len(targets) - 1].tick_params(labelbottom=True)


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
	figure.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.0),
	              ncol=len(periods), frameon=False, fontsize=9)


def finish_layout(figure, axes, targets_by_trend, periods):
	'''
	Reserve a fixed height for the headers and the legend, so a grid of any
	number of rows keeps the same margins.
	'''
	height = figure.get_figheight()
	bottom, top = FOOTER_INCHES / height, 1.0 - HEADER_INCHES / height
	figure.tight_layout(rect=[0.03, bottom, 1, top])
	add_column_headers(figure, axes, targets_by_trend, top)
	add_period_legend(figure, periods)
	figure.supxlabel('Thousand years before present', fontsize=11, y=bottom * 0.42)
	figure.supylabel('$F_{ST}$', fontsize=12)


def build_figure(table, periods, uncertainty):
	'''
	One column per expected trend, one panel per locus, all on a shared time
	axis running from oldest to most recent.
	'''
	targets_by_trend = trend_targets(table)
	row_count = max(len(targets) for targets in targets_by_trend.values())
	figure, axes = plt.subplots(
		row_count, len(TREND_COLUMNS), squeeze=False, sharex=True,
		figsize=(4.2 * len(TREND_COLUMNS), PANEL_INCHES * row_count))
	background_rows = table[table['target'] == GENOME_WIDE_TARGET].sort_values('time_start')
	for column_position, trend in enumerate(TREND_COLUMNS):
		fill_column(axes[:, column_position], targets_by_trend[trend], table,
		            background_rows, periods, uncertainty)
	axes[0, 0].set_xlim(table['time_end'].max() / 1000.0, 0)
	finish_layout(figure, axes, targets_by_trend, periods)
	return figure


def main():
	parser = argparse.ArgumentParser(description='All loci on one page, in a column per expected trend')
	parser.add_argument('--results', default='fst_time_series.csv')
	parser.add_argument('--periods', default='time_periods.json')
	parser.add_argument('--output-dir', default='.')
	parser.add_argument('--uncertainty', choices=UNCERTAINTIES, default=UNCERTAINTIES[0])
	args = parser.parse_args()
	table = pd.read_csv(args.results)
	figure = build_figure(table, load_periods(args.periods), args.uncertainty)
	figure.savefig(pathlib.Path(args.output_dir) / f'fst_trend_grid_{args.uncertainty}.png', dpi=200)


if __name__ == '__main__':
	main()
