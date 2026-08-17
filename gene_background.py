'''
The per gene results and the background they make.

Every gene of the annotation is measured, and every gene is written out, so a
later question about any set of genes is answered from the table rather than
from the genotypes again. The background drawn beside a focal locus is the
mean over genes, with an interval that is the spread between them.
'''

import numpy as np
import pandas as pd
from statistics import NormalDist

from fst_core import ratio_of_averages

GENE_TABLE_COLUMNS = ['gene', 'chromosome', 'start', 'end', 'time_start', 'time_end', 'fst', 'n_variants']
BACKGROUND_COLUMNS = ['time_start', 'time_end', 'fst', 'ci_low', 'ci_high', 'standard_error', 'n_genes']


def gene_fst(gene_accumulators):
	'''
	FST of every gene at every time bin, from the sums of the whole sample.
	'''
	return ratio_of_averages(
		gene_accumulators['sum_a'], gene_accumulators['sum_b'], gene_accumulators['sum_c'])


def gene_table(context, gene_accumulators):
	'''
	One row per gene and time bin, holding the gene's own FST and how many
	variants it rests on, so any set of genes can be selected later.
	'''
	values = gene_fst(gene_accumulators)
	counts = gene_accumulators['variant_count']
	names, chromosomes, starts, ends = zip(*context['genes'])
	frames = []
	for bin_position, time_bin in enumerate(context['time_bins']):
		frames.append(pd.DataFrame({
			'gene': names,
			'chromosome': chromosomes,
			'start': starts,
			'end': ends,
			'time_start': time_bin['time_start'],
			'time_end': time_bin['time_end'],
			'fst': values[bin_position],
			'n_variants': counts[bin_position].astype(int),
		}))
	return pd.concat(frames, ignore_index=True)[GENE_TABLE_COLUMNS]


def background_row(time_bin, values, confidence_level):
	'''
	The mean over genes at one time bin, with the interval that the genes
	themselves spread over.
	'''
	usable = values[~np.isnan(values)]
	standard_error = float(np.std(usable, ddof=1) / np.sqrt(len(usable))) if len(usable) > 1 else np.nan
	quantile = NormalDist().inv_cdf(1.0 - (1.0 - confidence_level) / 2.0)
	mean_value = float(np.mean(usable)) if len(usable) else np.nan
	return {
		'time_start': time_bin['time_start'],
		'time_end': time_bin['time_end'],
		'fst': mean_value,
		'ci_low': mean_value - quantile * standard_error,
		'ci_high': mean_value + quantile * standard_error,
		'standard_error': standard_error,
		'n_genes': len(usable),
	}


def background_table(context, gene_accumulators, min_gene_variants, confidence_level):
	'''
	The mean over genes at each time bin, counting only the genes holding
	enough variants to be worth a number of their own.
	'''
	values = gene_fst(gene_accumulators)
	counted = np.where(gene_accumulators['variant_count'] >= min_gene_variants, values, np.nan)
	rows = [
		background_row(time_bin, counted[bin_position], confidence_level)
		for bin_position, time_bin in enumerate(context['time_bins'])
	]
	return pd.DataFrame(rows, columns=BACKGROUND_COLUMNS)


def save_gene_sums(output_path, context, gene_accumulators):
	'''
	The sums each gene's FST is formed from, so a gene can be recombined with
	others without the genotypes.
	'''
	names, chromosomes, starts, ends = zip(*context['genes'])
	np.savez(
		output_path,
		gene=np.array(names),
		chromosome=np.array(chromosomes),
		start=np.array(starts),
		end=np.array(ends),
		time_start=np.array([time_bin['time_start'] for time_bin in context['time_bins']]),
		time_end=np.array([time_bin['time_end'] for time_bin in context['time_bins']]),
		**gene_accumulators)
