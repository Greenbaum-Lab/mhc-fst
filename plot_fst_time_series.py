import argparse
import pathlib
import pandas as pd
import matplotlib.pyplot as plt

from variant_masks import GENOME_WIDE_TARGET

BACKGROUND_COLOR = 'black'
BACKGROUND_LABEL = 'genome wide'
UNCERTAINTIES = {
	'individuals': 'jackknife over individuals',
	'snp_blocks': 'jackknife over blocks of variants',
}


def bin_midpoints(rows):
	return (rows['time_start'] + rows['time_end']) / 2000.0


def plot_target(axis, rows, label, color, line_style, uncertainty):
	'''
	One target as a line with its jackknife interval as a shaded band.
	'''
	years = bin_midpoints(rows)
	axis.plot(years, rows['fst'], label=label, color=color, linestyle=line_style)
	axis.scatter(years, rows['fst'], color=color, s=18)
	axis.fill_between(
		years, rows[f'ci_low_{uncertainty}'], rows[f'ci_high_{uncertainty}'],
		color=color, alpha=0.15, linewidth=0)


def annotate_sample_counts(axis, background_rows, uncertainty):
	'''
	Individuals of each population per bin, the same for every target.
	'''
	for _, row in background_rows.iterrows():
		axis.annotate(
			f'{row["n_samples_a"]}/{row["n_samples_b"]}',
			xy=((row['time_start'] + row['time_end']) / 2000.0, row[f'ci_low_{uncertainty}']),
			xytext=(0, -12), textcoords='offset points',
			fontsize=8, ha='center', color='0.4')


def draw_series(axis, rows, uncertainty):
	focal_targets = sorted(set(rows['target']) - {GENOME_WIDE_TARGET})
	colors = plt.cm.Dark2.colors
	for position, target in enumerate(focal_targets):
		target_rows = rows[rows['target'] == target].sort_values('time_start')
		plot_target(axis, target_rows, target, colors[position % len(colors)], '-', uncertainty)
	background_rows = rows[rows['target'] == GENOME_WIDE_TARGET].sort_values('time_start')
	plot_target(axis, background_rows, BACKGROUND_LABEL, BACKGROUND_COLOR, '--', uncertainty)
	annotate_sample_counts(axis, background_rows, uncertainty)


def style_axis(axis, title):
	axis.set_xlabel('Thousand years before present')
	axis.set_ylabel('$F_{ST}$')
	axis.set_title(title, fontsize=11)
	axis.invert_xaxis()
	axis.axhline(0.0, color='0.8', linewidth=0.8)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)
	axis.legend(frameon=False, fontsize=9)


def plot_panel(rows, title, output_path, uncertainty):
	figure, axis = plt.subplots(figsize=(8, 5))
	draw_series(axis, rows, uncertainty)
	style_axis(axis, title)
	figure.tight_layout()
	figure.savefig(output_path, dpi=200)
	plt.close(figure)


def panel_title(rows, locus, uncertainty):
	first = rows.iloc[0]
	return (f'{locus}, {first["polygon_a"]} vs {first["polygon_b"]}, '
	        f'{first["estimator"].replace("_", " ")}\n{UNCERTAINTIES[uncertainty]}')


def locus_panels(table):
	'''
	The rows of each figure: one locus against the genome wide background, for
	every estimator.
	'''
	loci = sorted(set(table['locus'].fillna('')) - {''})
	for estimator, rows in table.groupby('estimator'):
		background = rows[rows['target'] == GENOME_WIDE_TARGET]
		for locus in loci:
			yield locus, estimator, pd.concat([rows[rows['locus'] == locus], background])


def plot_all_panels(table, output_dir):
	'''
	One figure per locus, estimator and jackknife, each holding the spans of
	that locus and the genome wide background.
	'''
	for locus, estimator, rows in locus_panels(table):
		for uncertainty in UNCERTAINTIES:
			output_path = output_dir / f'fst_{locus}_{estimator}_{uncertainty}.png'
			plot_panel(rows, panel_title(rows, locus, uncertainty), output_path, uncertainty)


def main():
	parser = argparse.ArgumentParser(description='Plot FST across time for focal loci and the genome wide background')
	parser.add_argument('--results', required=True)
	parser.add_argument('--output-dir', required=True)
	args = parser.parse_args()
	output_dir = pathlib.Path(args.output_dir)
	output_dir.mkdir(parents=True, exist_ok=True)
	plot_all_panels(pd.read_csv(args.results), output_dir)


if __name__ == '__main__':
	main()
