'''
Reading the temporal FST scan table and preparing it for plotting.

A gene holding a handful of variants scatters far more than one holding
hundreds, so FST itself cannot be read across genes. The scan already answers
that by drawing a null from random variant sets of the same size as the gene,
and `null_percentile` is where the gene falls inside that null. That is the
quantity carried here, and the side of the null a gene falls on is its
direction.
'''

import numpy as np
import pandas as pd

AUTOSOMES = [str(number) for number in range(1, 23)]
HIGH = 'high'
LOW = 'low'
MIDDLE_OF_NULL = 0.5


def load_results(path):
	results = pd.read_csv(path, sep='\t')
	results['chrom'] = results['chrom'].astype(str)
	return results


def usable_genes(results, min_variants):
	'''
	Autosomal genes holding more than min_variants variants and carrying a
	finite position inside their own null.
	'''
	genes = results[results['chrom'].isin(AUTOSOMES)]
	genes = genes[genes['n_variants'] > min_variants]
	genes = genes[np.isfinite(genes['null_percentile'])]
	return genes.copy()


def add_direction(genes):
	'''
	Whether a gene sits in the upper or the lower half of its own null.
	'''
	genes = genes.copy()
	genes['direction'] = np.where(genes['null_percentile'] >= MIDDLE_OF_NULL, HIGH, LOW)
	return genes


def genomic_axis(genes):
	'''
	Gene midpoints laid end to end along the autosomes, returned with the
	centre of every chromosome for the ticks of the axis.
	'''
	chromosome_ends = genes.groupby('chrom')['end'].max()
	offsets = {}
	centres = {}
	offset = 0.0
	for chromosome in AUTOSOMES:
		if chromosome not in chromosome_ends.index:
			continue
		offsets[chromosome] = offset
		centres[chromosome] = offset + chromosome_ends[chromosome] / 2.0
		offset += float(chromosome_ends[chromosome])
	genes = genes.copy()
	genes['position'] = (genes['start'] + genes['end']) / 2.0 + genes['chrom'].map(offsets)
	return genes, centres


def significant_in_every_group(genes):
	'''
	The genes called outside their null in every group of the table.
	'''
	called = genes[genes['outside_null'] == 1]
	groups_per_gene = called.groupby('gene')['group'].nunique()
	shared = groups_per_gene[groups_per_gene == genes['group'].nunique()].index
	return called[called['gene'].isin(shared)].copy()


def prepare(path, min_variants):
	'''
	The usable genes of a scan table, with their direction and their place
	along the genome, and the chromosome centres of that axis.
	'''
	genes = add_direction(usable_genes(load_results(path), min_variants))
	return genomic_axis(genes)
