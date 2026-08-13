'''
Reading blocks of a PLINK BED fileset by sample and by variant.
'''

import numpy as np
from bed_reader import open_bed

AUTOSOMES = [str(number) for number in range(1, 23)]


def open_genotypes(bed_prefix, threads):
	return open_bed(f'{bed_prefix}.bed', num_threads=threads)


def autosomal_variants(source):
	'''
	File index, chromosome and base position of the autosomal variants.
	'''
	chromosome = np.array([str(value).replace('chr', '') for value in source.chromosome])
	position = np.asarray(source.bp_position)
	variant_index = np.where(np.isin(chromosome, AUTOSOMES))[0]
	return variant_index, chromosome[variant_index], position[variant_index]


def select_sample_indices(source, sample_ids):
	'''
	Source row index of each requested sample, in the order requested.
	Raises KeyError if any requested sample is absent from the source, so a
	population is never silently computed from part of its roster.
	'''
	index_by_id = {str(iid): position for position, iid in enumerate(source.iid)}
	missing = [sample_id for sample_id in sample_ids if sample_id not in index_by_id]
	if missing:
		raise KeyError(f'Samples not found in source BED: {sorted(missing)}')
	return np.array([index_by_id[sample_id] for sample_id in sample_ids])


def read_genotypes(source, sample_indices, variant_indices):
	'''
	Genotypes of the requested samples and variants as float32, with nan for
	missing calls, laid out as samples by variants.
	'''
	return source.read(index=np.s_[sample_indices, variant_indices], dtype='float32')


def allele_statistics(genotypes, weights):
	'''
	Allele count, allele number and heterozygote count per variant for every
	column of weights. Weighting the individuals is what makes a leave-one-out
	sample, so one matrix product yields every sample at once.
	'''
	called = np.isfinite(genotypes).astype(np.float32)
	filled = np.nan_to_num(genotypes, nan=0.0)
	heterozygous = (genotypes == 1.0).astype(np.float32)
	return filled.T @ weights, 2.0 * (called.T @ weights), heterozygous.T @ weights
