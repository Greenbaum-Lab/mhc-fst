import argparse
import pathlib
import pandas as pd
import matplotlib.pyplot as plt

from variant_masks import GENOME_WIDE_TARGET

FOCAL_COLOR = '#c0392b'
GENE_BACKGROUND_COLOR = '#2e7d32'
BACKGROUND_COLOR = 'black'
BACKGROUND_LABEL = 'genome wide'
GENE_BACKGROUND_LABEL = 'mean over all annotated genes'
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


def plot_gene_background(axis, gene_background):
	'''
	The mean over every gene of the annotation, with error bars spanning the
	genes rather than the uncertainty of any one of them.
	'''
	axis.errorbar(
		bin_midpoints(gene_background), gene_background['fst'],
		yerr=[gene_background['fst'] - gene_background['ci_low'],
		      gene_background['ci_high'] - gene_background['fst']],
		color=GENE_BACKGROUND_COLOR, linewidth=1.2, marker='s', markersize=4,
		capsize=3, elinewidth=1.0, label=GENE_BACKGROUND_LABEL)


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


def style_axis(axis, title):
	axis.set_xlabel('Thousand years before present')
	axis.set_ylabel('$F_{ST}$')
	axis.set_title(title, fontsize=11)
	axis.invert_xaxis()
	axis.axhline(0.0, color='0.8', linewidth=0.8)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)
	axis.legend(frameon=False, fontsize=9)


def panel_title(target_rows):
	'''
	The locus, what it is studied for, and the trend expected of it.
	'''
	first = target_rows.iloc[0]
	return f'{first["target"]}, {first["phenotype"]}, expected {first["trend"]}'


def plot_panel(target_rows, background_rows, gene_background, output_path, uncertainty):
	figure, axis = plt.subplots(figsize=(8, 5))
	plot_target(axis, target_rows, target_rows.iloc[0]['target'], FOCAL_COLOR, '-', uncertainty)
	plot_target(axis, background_rows, BACKGROUND_LABEL, BACKGROUND_COLOR, '--', uncertainty)
	plot_gene_background(axis, gene_background)
	annotate_sample_counts(axis, background_rows, uncertainty)
	style_axis(axis, panel_title(target_rows))
	figure.tight_layout()
	figure.savefig(output_path, dpi=200)
	plt.close(figure)


def plot_all_panels(table, gene_background, output_dir):
	'''
	One figure per locus and jackknife, each holding that locus, the genome
	wide background and the mean over all annotated genes.
	'''
	background_rows = table[table['target'] == GENOME_WIDE_TARGET].sort_values('time_start')
	for target in sorted(set(table['target']) - {GENOME_WIDE_TARGET}):
		target_rows = table[table['target'] == target].sort_values('time_start')
		for uncertainty in UNCERTAINTIES:
			plot_panel(target_rows, background_rows, gene_background,
			           output_dir / f'fst_{target}_{uncertainty}.png', uncertainty)


def main():
	parser = argparse.ArgumentParser(description='Plot FST across time for focal loci and the genome wide background')
	parser.add_argument('--results', required=True)
	parser.add_argument('--gene-background', required=True)
	parser.add_argument('--output-dir', required=True)
	args = parser.parse_args()
	output_dir = pathlib.Path(args.output_dir)
	output_dir.mkdir(parents=True, exist_ok=True)
	gene_background = pd.read_csv(args.gene_background).sort_values('time_start')
	plot_all_panels(pd.read_csv(args.results), gene_background, output_dir)


if __name__ == '__main__':
	main()
