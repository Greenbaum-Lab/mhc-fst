import argparse
import pathlib
import pandas as pd
import matplotlib.pyplot as plt

from variant_masks import GENOME_WIDE_TARGET

BACKGROUND_COLOR = 'black'
BACKGROUND_LABEL = 'genome wide'


def bin_midpoints(rows):
	return (rows['time_start'] + rows['time_end']) / 2000.0


def plot_target(axis, rows, label, color, line_style):
	'''
	One target as a line with its jackknife interval as a shaded band.
	'''
	years = bin_midpoints(rows)
	axis.plot(years, rows['fst'], label=label, color=color, linestyle=line_style)
	axis.scatter(years, rows['fst'], color=color, s=18)
	axis.fill_between(years, rows['ci_low'], rows['ci_high'], color=color, alpha=0.15, linewidth=0)


def annotate_sample_counts(axis, background_rows):
	'''
	Individuals of each population per bin, the same for every target.
	'''
	for _, row in background_rows.iterrows():
		axis.annotate(
			f'{row["n_samples_a"]}/{row["n_samples_b"]}',
			xy=((row['time_start'] + row['time_end']) / 2000.0, row['ci_low']),
			xytext=(0, -12), textcoords='offset points',
			fontsize=8, ha='center', color='0.4')


def draw_series(axis, rows):
	focal_targets = sorted(set(rows['target']) - {GENOME_WIDE_TARGET})
	colors = plt.cm.Dark2.colors
	for position, target in enumerate(focal_targets):
		target_rows = rows[rows['target'] == target].sort_values('time_start')
		plot_target(axis, target_rows, target, colors[position % len(colors)], '-')
	background_rows = rows[rows['target'] == GENOME_WIDE_TARGET].sort_values('time_start')
	plot_target(axis, background_rows, BACKGROUND_LABEL, BACKGROUND_COLOR, '--')
	annotate_sample_counts(axis, background_rows)


def style_axis(axis, title):
	axis.set_xlabel('Thousand years before present')
	axis.set_ylabel('$F_{ST}$')
	axis.set_title(title, fontsize=11)
	axis.invert_xaxis()
	axis.axhline(0.0, color='0.8', linewidth=0.8)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)
	axis.legend(frameon=False, fontsize=9)


def plot_panel(rows, title, output_path):
	figure, axis = plt.subplots(figsize=(8, 5))
	draw_series(axis, rows)
	style_axis(axis, title)
	figure.tight_layout()
	figure.savefig(output_path, dpi=200)
	plt.close(figure)


def panel_title(rows, gene):
	first = rows.iloc[0]
	return (f'{gene}, {first["polygon_a"]} vs {first["polygon_b"]}, '
	        f'{first["estimator"].replace("_", " ")}, SNP set: {first["filter_mode"]}')


def gene_panels(table):
	'''
	The rows of each figure: one gene against the genome wide background, for
	every SNP set alternative and estimator.
	'''
	genes = sorted(set(table['gene'].fillna('')) - {''})
	for (filter_mode, estimator), rows in table.groupby(['filter_mode', 'estimator']):
		background = rows[rows['target'] == GENOME_WIDE_TARGET]
		for gene in genes:
			panel_rows = pd.concat([rows[rows['gene'] == gene], background])
			yield gene, filter_mode, estimator, panel_rows


def plot_all_panels(table, output_dir):
	'''
	One figure per gene, SNP set alternative and estimator, each holding the
	spans of that gene and the genome wide background.
	'''
	for gene, filter_mode, estimator, rows in gene_panels(table):
		output_path = output_dir / f'fst_{gene}_{filter_mode}_{estimator}.png'
		plot_panel(rows, panel_title(rows, gene), output_path)


def main():
	parser = argparse.ArgumentParser(description='Plot FST across time for focal regions and the genome wide background')
	parser.add_argument('--results', required=True)
	parser.add_argument('--output-dir', required=True)
	args = parser.parse_args()
	output_dir = pathlib.Path(args.output_dir)
	output_dir.mkdir(parents=True, exist_ok=True)
	plot_all_panels(pd.read_csv(args.results), output_dir)


if __name__ == '__main__':
	main()
