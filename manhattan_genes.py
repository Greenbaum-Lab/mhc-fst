'''
A Manhattan plot of every gene, one per time bin.

FST on its own does not say whether a gene stands out, because a gene measured
from twenty variants scatters far more than one measured from six hundred.
Each gene is ranked against the genes holding a similar number of variants,
so that gene size does not decide which genes reach the top.

The null is the bulk of the genes themselves. Within a stratum the centre and
the spread are taken robustly, from the median and the median absolute
deviation, so that the outliers being looked for do not inflate the very
spread they are measured against. A gene's p value is the normal tail beyond
its distance from that centre, which is what lets a gene reach further than
one over the number of genes and lets Benjamini Hochberg call anything at all.

Ranking genes against their own empirical distribution instead, as a fraction
of genes above them, cannot work here: those p values are uniform by
construction, so every q value comes back at one whatever the data hold. The
empirical fraction is still written out, as `percentile`, because it says
plainly where a gene sits among the rest.

Two assumptions come with this. Most genes are taken to be ordinary, which is
what makes the bulk a null. And the tail is taken to be normal, which the
skew of an FST distribution only approximately obeys, so a p value here is a
way of ordering candidates rather than a precise probability.
'''

import math
import argparse
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from manhattan import assign_cumulative_positions
from outlier_scan import empirical_one_sided_p, benjamini_hochberg

MAD_TO_DEVIATION = 1.4826
SMALLEST_TAIL = np.finfo(float).tiny
CHROMOSOMES = [str(number) for number in range(1, 23)]
CHROMOSOME_COLORS = ['#3b6ea5', '#8a8a8a']
THRESHOLD_COLOR = '#c0392b'
SIDES = ['high', 'low']


def usable_genes(per_gene, min_variants):
	'''
	Autosomal genes holding enough variants to be worth a number of their own.
	'''
	genes = per_gene[per_gene['chromosome'].astype(str).isin(CHROMOSOMES)]
	genes = genes[(genes['n_variants'] >= min_variants) & np.isfinite(genes['fst'])]
	return genes.copy()


def normal_upper_tail(deviations):
	'''
	The probability of exceeding each deviation under a standard normal, held
	above zero so that a gene far enough out still has a finite score.
	'''
	tails = np.array([0.5 * math.erfc(value / math.sqrt(2.0)) for value in deviations])
	return np.maximum(tails, SMALLEST_TAIL)


def robust_null(values):
	'''
	Centre and spread of the ordinary genes, taken from the median and the
	median absolute deviation so that outliers do not widen the null they are
	measured against. Where half the genes of a stratum share one value the
	deviation is zero and the plain standard deviation stands in for it.
	'''
	centre = np.median(values)
	spread = MAD_TO_DEVIATION * np.median(np.abs(values - centre))
	return centre, spread if spread > 0.0 else np.std(values)


def stratum_tails(values):
	'''
	One sided p values of each gene against the robust null of its stratum,
	and the fraction of the stratum it sits above.
	'''
	centre, spread = robust_null(values)
	deviations = (values - centre) / spread
	percentile_high, _ = empirical_one_sided_p(values)
	return normal_upper_tail(deviations), normal_upper_tail(-deviations), percentile_high


def stratified_p_values(genes, strata):
	'''
	Rank every gene against the genes holding a similar number of variants, so
	that gene size does not decide which genes reach the top.
	'''
	genes = genes.copy()
	genes['stratum'] = pd.qcut(genes['n_variants'], strata, labels=False, duplicates='drop')
	for column in ('p_high', 'p_low', 'percentile'):
		genes[column] = np.nan
	for stratum in genes['stratum'].unique():
		rows = genes['stratum'] == stratum
		p_high, p_low, percentile = stratum_tails(genes.loc[rows, 'fst'].to_numpy())
		genes.loc[rows, 'p_high'] = p_high
		genes.loc[rows, 'p_low'] = p_low
		genes.loc[rows, 'percentile'] = 1.0 - percentile
	genes['q_high'] = benjamini_hochberg(genes['p_high'].to_numpy())
	genes['q_low'] = benjamini_hochberg(genes['p_low'].to_numpy())
	return genes


def genomic_axis(genes):
	'''
	Every gene on one axis, chromosome after chromosome, ordered within each.
	'''
	frames = []
	for chromosome in CHROMOSOMES:
		rows = genes[genes['chromosome'].astype(str) == chromosome].copy()
		if rows.empty:
			continue
		rows['position'] = (rows['start'] + rows['end']) / 2.0
		frames.append(rows.sort_values('position'))
	return assign_cumulative_positions(frames)


def significance_threshold(placed, side, q_threshold):
	'''
	The largest p value Benjamini Hochberg still calls at this threshold, or
	nan where it calls nothing.
	'''
	called = placed[placed[f'q_{side}'] <= q_threshold]
	return called[f'p_{side}'].max() if not called.empty else np.nan


def top_genes(placed, side, count):
	'''
	The genes furthest into the tail, most extreme FST breaking any ties.
	'''
	ordered = placed.sort_values([f'p_{side}', 'fst'], ascending=[True, side == 'low'])
	return ordered.head(count)


def add_scores(placed, side, max_score):
	'''
	The plotted height of each gene, held at the cap so that one gene far out
	in the tail does not flatten every other. Genes at the cap are marked, as
	their true height is off the top of the axis.
	'''
	scores = -np.log10(placed[f'p_{side}'])
	placed = placed.copy()
	placed['score'] = np.minimum(scores, max_score)
	placed['capped'] = scores > max_score
	return placed


def draw_points(axis, placed):
	'''
	Every gene as a point, chromosome by chromosome, the capped ones as
	triangles.
	'''
	positions = placed['cumulative_position'].to_numpy()
	scores = placed['score'].to_numpy()
	capped = placed['capped'].to_numpy()
	for position, chromosome in enumerate(sorted(set(placed['chromosome'].astype(str)), key=int)):
		rows = (placed['chromosome'].astype(str) == chromosome).to_numpy()
		color = CHROMOSOME_COLORS[position % len(CHROMOSOME_COLORS)]
		axis.scatter(positions[rows & ~capped], scores[rows & ~capped], s=6, linewidths=0, color=color)
		axis.scatter(positions[rows & capped], scores[rows & capped], s=28, marker='^', linewidths=0, color=color)


def annotate_genes(axis, top):
	for _, row in top.iterrows():
		axis.annotate(
			row['gene'], xy=(row['cumulative_position'], row['score']),
			xytext=(0, 5), textcoords='offset points', rotation=90,
			fontsize=6, ha='center', va='bottom')


def style_axis(axis, placed, side, title):
	ticks = placed.groupby(placed['chromosome'].astype(str))['cumulative_position'].mean()
	ticks = ticks.reindex(sorted(ticks.index, key=int))
	axis.set_xticks(ticks.to_numpy())
	axis.set_xticklabels(ticks.index, fontsize=8)
	axis.set_xlabel('Chromosome')
	axis.set_ylabel(f'$-\\log_{{10}}$ p, FST {side} for gene size')
	axis.set_title(title, fontsize=11)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)
	axis.margins(x=0.01)


def plot_bin(placed, side, q_threshold, annotate_count, title, output_path):
	'''
	One time bin as a Manhattan plot, with the most extreme genes named.
	'''
	figure, axis = plt.subplots(figsize=(16, 5.5))
	draw_points(axis, placed)
	threshold = significance_threshold(placed, side, q_threshold)
	if np.isfinite(threshold) and threshold > 0.0:
		axis.axhline(-np.log10(threshold), color=THRESHOLD_COLOR, linestyle='--', linewidth=1)
	top = top_genes(placed, side, annotate_count)
	annotate_genes(axis, top)
	style_axis(axis, placed, side, title)
	axis.set_ylim(top=float(placed['score'].max()) * 1.35)
	figure.tight_layout()
	figure.savefig(output_path, dpi=200)
	plt.close(figure)
	return top


def scan_bin(per_gene, time_start, time_end, arguments, output_dir):
	'''
	Rank the genes of one time bin, draw them and return their table.
	'''
	rows = per_gene[(per_gene['time_start'] == time_start) & (per_gene['time_end'] == time_end)]
	placed = genomic_axis(stratified_p_values(usable_genes(rows, arguments.min_variants), arguments.strata))
	placed = add_scores(placed, arguments.side, arguments.max_score)
	title = (f'{time_start}-{time_end} years before present, {len(placed)} genes, '
	         f'{arguments.strata} strata by variant count, '
	         f'{int(placed["capped"].sum())} above the axis')
	top = plot_bin(
		placed, arguments.side, arguments.q_threshold, arguments.annotate, title,
		output_dir / f'manhattan_{arguments.side}_{time_start}_{time_end}.png')
	return placed, top


def main():
	parser = argparse.ArgumentParser(description='Manhattan plot of every gene, one figure per time bin')
	parser.add_argument('--per-gene', default='results/fst_per_gene.csv')
	parser.add_argument('--output-dir', default='results/manhattan')
	parser.add_argument('--min-variants', type=int, default=10)
	parser.add_argument('--strata', type=int, default=10)
	parser.add_argument('--annotate', type=int, default=50)
	parser.add_argument('--q-threshold', type=float, default=0.05)
	parser.add_argument('--max-score', type=float, default=30.0)
	parser.add_argument('--side', choices=SIDES, default='high')
	arguments = parser.parse_args()
	output_dir = pathlib.Path(arguments.output_dir)
	output_dir.mkdir(parents=True, exist_ok=True)
	per_gene = pd.read_csv(arguments.per_gene)
	bins = per_gene[['time_start', 'time_end']].drop_duplicates().sort_values('time_start')
	scanned = [scan_bin(per_gene, start, end, arguments, output_dir) for start, end in bins.to_numpy()]
	pd.concat([placed for placed, _ in scanned]).to_csv(output_dir / 'gene_significance.csv', index=False)
	pd.concat([top for _, top in scanned]).to_csv(output_dir / 'top_genes.csv', index=False)
	print(f'{len(bins)} time bins, {arguments.annotate} genes named in each')


if __name__ == '__main__':
	main()
