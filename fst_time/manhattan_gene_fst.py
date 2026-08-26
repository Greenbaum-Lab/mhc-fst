'''
A Manhattan plot of gene FST, one pair of figures per time bin.

Significance here is a plain quantile of the genes themselves: the lowest and
the highest three percent of the bin stand out, whatever their spread. This
asks only where a gene sits among the rest, and unlike the stratified tail of
manhattan_genes.py it does not account for how many variants a gene holds.
Every gene with a finite FST enters, whatever the number of variants behind it.

Two figures are drawn from the same genes. The first is FST as measured. The
second is FST over the genome wide FST of the same bin, which places a gene
against the background its own bin sets and so makes bins comparable to one
another.

The FST of genes crowds towards its floor, which leaves the low genes on top
of one another on a plain axis. Two other axes are offered for reading that
end. A symlog axis is linear near zero and logarithmic further out, so it
keeps FST units while giving the crowded end room. A rankit axis replaces each
gene by the normal score of its rank, which spreads the genes evenly at the
cost of the axis no longer being FST. Both keep the order of the genes exactly,
so the same genes fall in the extreme quantiles whichever axis is drawn.
'''

import argparse
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statistics import NormalDist

from fst_time.manhattan_genes import usable_genes, genomic_axis
from fst_time.variant_masks import GENOME_WIDE_TARGET

CHROMOSOME_COLORS = ['#9a9a9a', '#c8c8c8']
EXTREME_COLORS = {'low': '#2980b9', 'high': '#c0392b'}
THRESHOLD_COLOR = '#c0392b'
MEASURES = {'fst': 'FST', 'ratio': 'FST / genome wide FST'}
LINEAR, SYMLOG, RANKIT = 'linear', 'symlog', 'rankit'
SCALES = [LINEAR, SYMLOG, RANKIT]
SYMLOG_THRESHOLD = 0.01
PLOTTED_COLUMN = 'plotted'
NO_VARIANT_FLOOR = 0


def genome_wide_fst(results):
	'''
	The genome wide FST of every time bin, keyed by the bin.
	'''
	rows = results[results['target'] == GENOME_WIDE_TARGET]
	return rows.set_index(['time_start', 'time_end'])['fst']


def add_ratio(genes, genome_wide_value):
	'''
	Each gene against the genome wide FST of its own time bin.
	'''
	genes = genes.copy()
	genes['ratio'] = genes['fst'] / genome_wide_value
	return genes


def normal_scores(values):
	'''
	The normal score of each value's rank, which spreads a crowded bulk evenly
	over the axis while keeping the order of the genes exactly.
	'''
	quantiles = pd.Series(values).rank(method='average') / (len(values) + 1.0)
	return np.array([NormalDist().inv_cdf(quantile) for quantile in quantiles])


def plotted_values(genes, measure, scale):
	'''
	The column the figure draws, which is the measure itself unless the axis
	is a rankit and asks for normal scores instead.
	'''
	values = genes[measure].to_numpy()
	if scale == RANKIT:
		values = normal_scores(values)
	return genes.assign(**{PLOTTED_COLUMN: values})


def label_extremes(genes, measure, fraction):
	'''
	Mark every gene low, high or ordinary by where it falls against the two
	quantiles of the bin, and return those quantiles beside it.
	'''
	values = genes[PLOTTED_COLUMN].to_numpy()
	low, high = np.quantile(values, fraction), np.quantile(values, 1.0 - fraction)
	genes = genes.copy()
	genes[f'side_{measure}'] = np.where(
		values <= low, 'low', np.where(values >= high, 'high', 'ordinary'))
	return genes, low, high


def draw_points(axis, placed, measure):
	'''
	Every gene as a point, the colour alternating by chromosome, with the
	extreme genes of the bin coloured by the side they fall on.
	'''
	sides = placed[f'side_{measure}']
	for position, chromosome in enumerate(sorted(set(placed['chromosome'].astype(str)), key=int)):
		rows = placed[(placed['chromosome'].astype(str) == chromosome) & (sides == 'ordinary')]
		axis.scatter(rows['cumulative_position'], rows[PLOTTED_COLUMN], s=6, linewidths=0,
		             color=CHROMOSOME_COLORS[position % len(CHROMOSOME_COLORS)])
	for side, color in EXTREME_COLORS.items():
		rows = placed[sides == side]
		axis.scatter(rows['cumulative_position'], rows[PLOTTED_COLUMN], s=10, linewidths=0,
		             color=color, label=side)


def top_genes(placed, count):
	'''
	The genes furthest out on each side.
	'''
	ordered = placed.sort_values(PLOTTED_COLUMN)
	return pd.concat([ordered.head(count), ordered.tail(count)])


def annotate_genes(axis, top):
	for _, row in top.iterrows():
		axis.annotate(
			row['gene'], xy=(row['cumulative_position'], row[PLOTTED_COLUMN]),
			xytext=(0, 5), textcoords='offset points', rotation=90,
			fontsize=6, ha='center', va='bottom')


def axis_label(measure, scale):
	'''
	What the y axis holds, which is the measure unless the rankit axis has
	replaced it by a rank.
	'''
	if scale == RANKIT:
		return f'{MEASURES[measure]}, normal score of rank'
	return MEASURES[measure]


def style_axis(axis, placed, measure, scale, title):
	ticks = placed.groupby(placed['chromosome'].astype(str))['cumulative_position'].mean()
	ticks = ticks.reindex(sorted(ticks.index, key=int))
	axis.set_xticks(ticks.to_numpy())
	axis.set_xticklabels(ticks.index, fontsize=8)
	axis.set_xlabel('Chromosome')
	axis.set_ylabel(axis_label(measure, scale))
	if scale == SYMLOG:
		axis.set_yscale('symlog', linthresh=SYMLOG_THRESHOLD)
	axis.set_title(title, fontsize=11)
	axis.legend(loc='upper right', markerscale=2, frameon=False)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)
	axis.margins(x=0.01, y=0.12)


def plot_measure(genes, measure, scale, fraction, annotate_count, title, output_path):
	'''
	One measure of one time bin as a Manhattan plot, the extreme genes named
	and the two quantiles they are called against drawn across the figure.
	'''
	placed, low, high = label_extremes(plotted_values(genes, measure, scale), measure, fraction)
	figure, axis = plt.subplots(figsize=(16, 5.5))
	draw_points(axis, placed, measure)
	for threshold in (low, high):
		axis.axhline(threshold, color=THRESHOLD_COLOR, linestyle='--', linewidth=1)
	annotate_genes(axis, top_genes(placed, annotate_count))
	style_axis(axis, placed, measure, scale, title)
	figure.tight_layout()
	figure.savefig(output_path, dpi=200)
	plt.close(figure)
	return placed.drop(columns=PLOTTED_COLUMN)


def scan_bin(per_gene, genome_wide, time_bin, arguments, output_dir):
	'''
	Draw both measures of one time bin and return the genes behind them.
	'''
	time_start, time_end = time_bin
	rows = per_gene[(per_gene['time_start'] == time_start) & (per_gene['time_end'] == time_end)]
	genes = add_ratio(genomic_axis(usable_genes(rows, NO_VARIANT_FLOOR)), genome_wide[time_bin])
	title = (f'{time_start}-{time_end} years before present, {len(genes)} genes, '
	         f'extreme {arguments.fraction:.0%} on each side, {arguments.scale} axis')
	for measure in MEASURES:
		genes = plot_measure(
			genes, measure, arguments.scale, arguments.fraction, arguments.annotate, title,
			output_dir / f'manhattan_gene_{measure}_{arguments.scale}_{time_start}_{time_end}.png')
	return genes


def main():
	parser = argparse.ArgumentParser(description='Manhattan plot of gene FST and of gene FST over the genome wide FST')
	parser.add_argument('--per-gene', default='results/fst_per_gene.csv')
	parser.add_argument('--results', default='results/fst_time_series.csv')
	parser.add_argument('--output-dir', default='results/manhattan')
	parser.add_argument('--annotate', type=int, default=25)
	parser.add_argument('--fraction', type=float, default=0.03)
	parser.add_argument('--scale', choices=SCALES, default=LINEAR)
	arguments = parser.parse_args()
	output_dir = pathlib.Path(arguments.output_dir)
	output_dir.mkdir(parents=True, exist_ok=True)
	per_gene = pd.read_csv(arguments.per_gene)
	genome_wide = genome_wide_fst(pd.read_csv(arguments.results))
	bins = per_gene[['time_start', 'time_end']].drop_duplicates().sort_values('time_start')
	scanned = [
		scan_bin(per_gene, genome_wide, tuple(time_bin), arguments, output_dir)
		for time_bin in bins.to_numpy()
	]
	pd.concat(scanned).to_csv(output_dir / 'gene_extremes.csv', index=False)
	print(f'{len(bins)} time bins, two figures in each')


if __name__ == '__main__':
	main()
