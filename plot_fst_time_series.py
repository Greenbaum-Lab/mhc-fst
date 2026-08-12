import argparse
import pathlib
import pandas as pd
import matplotlib.pyplot as plt

from variant_masks import GENOME_WIDE_TARGET

BACKGROUND_COLOR = 'black'
BACKGROUND_LABEL = 'genome wide'
INTERVALS = ['basic', 'percentile']


def bin_midpoints(rows):
	return (rows['time_start'] + rows['time_end']) / 2000.0


def plot_target(axis, rows, label, color, line_style, interval):
	'''
	One target as a line with its bootstrap interval as a shaded band.
	'''
	years = bin_midpoints(rows)
	axis.plot(years, rows['fst'], label=label, color=color, linestyle=line_style)
	axis.scatter(years, rows['fst'], color=color, s=18)
	axis.fill_between(
		years, rows[f'ci_low_{interval}'], rows[f'ci_high_{interval}'],
		color=color, alpha=0.15, linewidth=0)


def annotate_sample_counts(axis, background_rows, interval):
	'''
	Individuals of each population per bin, the same for every target.
	'''
	for _, row in background_rows.iterrows():
		axis.annotate(
			f'{row["n_samples_a"]}/{row["n_samples_b"]}',
			xy=((row['time_start'] + row['time_end']) / 2000.0, row[f'ci_low_{interval}']),
			xytext=(0, -12), textcoords='offset points',
			fontsize=8, ha='center', color='0.4')


def draw_series(axis, rows, interval):
	focal_targets = sorted(set(rows['target']) - {GENOME_WIDE_TARGET})
	colors = plt.cm.Dark2.colors
	for position, target in enumerate(focal_targets):
		target_rows = rows[rows['target'] == target].sort_values('time_start')
		plot_target(axis, target_rows, target, colors[position % len(colors)], '-', interval)
	background_rows = rows[rows['target'] == GENOME_WIDE_TARGET].sort_values('time_start')
	plot_target(axis, background_rows, BACKGROUND_LABEL, BACKGROUND_COLOR, '--', interval)
	annotate_sample_counts(axis, background_rows, interval)


def style_axis(axis, title):
	axis.set_xlabel('Thousand years before present')
	axis.set_ylabel('$F_{ST}$')
	axis.set_title(title, fontsize=11)
	axis.invert_xaxis()
	axis.axhline(0.0, color='0.8', linewidth=0.8)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)
	axis.legend(frameon=False, fontsize=9)


def plot_panel(rows, title, output_path, interval):
	figure, axis = plt.subplots(figsize=(8, 5))
	draw_series(axis, rows, interval)
	style_axis(axis, title)
	figure.tight_layout()
	figure.savefig(output_path, dpi=200)
	plt.close(figure)


def panel_title(rows, interval):
	first = rows.iloc[0]
	return (f'{first["polygon_a"]} vs {first["polygon_b"]}, '
	        f'{first["estimator"].replace("_", " ")}, SNP set: {first["filter_mode"]}, '
	        f'{interval} bootstrap interval')


def plot_all_panels(table, output_dir, interval):
	'''
	One figure per SNP set alternative and estimator, so the effect of both
	choices is visible side by side.
	'''
	for (filter_mode, estimator), rows in table.groupby(['filter_mode', 'estimator']):
		output_path = output_dir / f'fst_time_series_{filter_mode}_{estimator}.png'
		plot_panel(rows, panel_title(rows, interval), output_path, interval)


def main():
	parser = argparse.ArgumentParser(description='Plot FST across time for focal regions and the genome wide background')
	parser.add_argument('--results', required=True)
	parser.add_argument('--output-dir', required=True)
	parser.add_argument('--interval', choices=INTERVALS, default='basic')
	args = parser.parse_args()
	output_dir = pathlib.Path(args.output_dir)
	output_dir.mkdir(parents=True, exist_ok=True)
	plot_all_panels(pd.read_csv(args.results), output_dir, args.interval)


if __name__ == '__main__':
	main()
