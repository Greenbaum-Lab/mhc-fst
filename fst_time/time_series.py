'''
FST of focal regions and of the genome across a series of time bins, for two
populations, with a jackknife over individuals and over blocks of variants.

The genotypes are read once. Each chunk of variants yields, for every time
bin, the Weir & Cockerham components of every leave-one-out sample of the
individuals, which accumulate into one sum per region. The same chunk also
accumulates the components of the whole sample per block of variants, so that
deleting a block later is a subtraction rather than another pass, and per gene
of the annotation, which gives a background made of genes. Focal regions and
the genome wide background share one set of leave-one-out samples per bin, so
their intervals are comparable.
'''

import numpy as np

from fst_time.fst_core import weir_cockerham_components
from fst_time.jackknife import paired_jackknife_weights, column_count, POINT_ESTIMATE_COLUMN
from fst_time.gene_regions import load_gene_spans, build_regions, locus_gene_names
from fst_time.time_populations import load_populations, overlapping_time_bins
from fst_time.variant_masks import region_masks, block_indices, GENOME_WIDE_TARGET
from fst_time.annotation_genes import (
	load_all_genes,
	gene_membership,
	chunk_pairs,
	empty_gene_accumulators,
	add_genes,
)
from fst_time.genotype_source import (
	open_genotypes,
	autosomal_variants,
	select_sample_indices,
	read_genotypes,
	allele_statistics,
)

ACCUMULATOR_NAMES = ['sum_a', 'sum_b', 'sum_c', 'variant_count']


def chunk_bounds(variant_count, chunk_size):
	return [(start, min(start + chunk_size, variant_count)) for start in range(0, variant_count, chunk_size)]


def union_sample_ids(time_bins):
	'''
	Every sample used by any bin, and the row each one occupies in a read
	block, so the genotypes are read once and shared by all bins.
	'''
	sample_ids = sorted({
		sample_id
		for time_bin in time_bins
		for sample_id in time_bin['samples_a'] + time_bin['samples_b']
	})
	return sample_ids, {sample_id: row for row, sample_id in enumerate(sample_ids)}


def sample_rows(row_by_sample_id, sample_ids):
	return np.array([row_by_sample_id[sample_id] for sample_id in sample_ids])


def bin_weights(time_bins):
	'''
	One leave-one-out weight matrix per population per bin.
	'''
	return [
		paired_jackknife_weights(len(time_bin['samples_a']), len(time_bin['samples_b']))
		for time_bin in time_bins
	]


def bin_columns(time_bins):
	'''
	Columns each bin fills, the widest of which sizes the accumulators.
	'''
	return [column_count(len(time_bin['samples_a']), len(time_bin['samples_b'])) for time_bin in time_bins]


def empty_accumulators(bin_count, target_count, width):
	return {name: np.zeros((bin_count, target_count, width)) for name in ACCUMULATOR_NAMES}


def add_selection(accumulators, index, components, selected):
	'''
	Add the variants of one region into the running sums, dropping the
	variants the estimator marked unusable. A bin fills only as many columns
	as it has individuals, so each sum is written into the leading columns.
	'''
	component_a, component_b, component_c = components
	columns = component_a.shape[1]
	accumulators['sum_a'][index][:columns] += np.nansum(component_a[selected], axis=0, dtype=np.float64)
	accumulators['sum_b'][index][:columns] += np.nansum(component_b[selected], axis=0, dtype=np.float64)
	accumulators['sum_c'][index][:columns] += np.nansum(component_c[selected], axis=0, dtype=np.float64)
	accumulators['variant_count'][index][:columns] += np.count_nonzero(~np.isnan(component_a[selected]), axis=0)


def add_blocks(block_accumulators, index, components, blocks, selected, block_count):
	'''
	Add the variants of one region into per block sums of the whole sample, so
	that the estimate without a block is the total minus that block.
	'''
	component_a, component_b, component_c = components
	column = POINT_ESTIMATE_COLUMN
	for name, values in (
		('sum_a', component_a[selected, column]),
		('sum_b', component_b[selected, column]),
		('sum_c', component_c[selected, column]),
		('variant_count', ~np.isnan(component_a[selected, column])),
	):
		weights = np.nan_to_num(values.astype(np.float64))
		block_accumulators[name][index] += np.bincount(blocks, weights=weights, minlength=block_count)


def accumulate_bin(accumulators, block_accumulators, bin_position, components, chunk_blocks, target_names, block_count):
	for target_position, target_name in enumerate(target_names):
		indices = chunk_blocks[target_name]
		selected = indices >= 0
		if not selected.any():
			continue
		index = (bin_position, target_position)
		add_selection(accumulators, index, components, selected)
		add_blocks(block_accumulators, index, components, indices[selected], selected, block_count)


def accumulate_chunk(accumulators, genotypes, context, weights_by_bin, chunk_blocks, chunk_genes, block_count):
	for bin_position, (rows_a, rows_b) in enumerate(context['bin_rows']):
		weights_a, weights_b = weights_by_bin[bin_position]
		count_a, number_a, het_a = allele_statistics(genotypes[rows_a], weights_a)
		count_b, number_b, het_b = allele_statistics(genotypes[rows_b], weights_b)
		components = weir_cockerham_components(number_a, number_b, count_a, count_b, het_a, het_b)
		accumulate_bin(
			accumulators['targets'], accumulators['blocks'], bin_position, components,
			chunk_blocks, context['target_names'], block_count)
		add_genes(
			accumulators['genes'], bin_position, components,
			chunk_genes[0], chunk_genes[1], len(context['genes']))


def build_context(config):
	'''
	Everything a run is defined by: the time bins compared, the focal regions,
	the variants and the samples they are read from.
	'''
	populations = load_populations(config['populations_path'])
	time_bins = overlapping_time_bins(populations, config['polygon_a'], config['polygon_b'], config['genotype_source'])
	gene_spans = load_gene_spans(config['annotation_path'], locus_gene_names(config['loci']))
	regions = build_regions(config['loci'], gene_spans)
	source = open_genotypes(config['bed_prefix'], config['threads'])
	variant_index, chromosome, position = autosomal_variants(source)
	sample_ids, row_by_sample_id = union_sample_ids(time_bins)
	genes = load_all_genes(config['annotation_path'], config['gene_biotypes'])
	return {
		'source': source,
		'time_bins': time_bins,
		'regions': regions,
		'variant_index': variant_index,
		'genes': genes,
		'gene_membership': gene_membership(genes, chromosome, position),
		'block_indices': block_indices(region_masks(regions, chromosome, position), config['snp_block_count']),
		'target_names': [region['locus'] for region in regions] + [GENOME_WIDE_TARGET],
		'region_by_target': {region['locus']: region for region in regions},
		'sample_indices': select_sample_indices(source, sample_ids),
		'bin_rows': [
			(sample_rows(row_by_sample_id, time_bin['samples_a']), sample_rows(row_by_sample_id, time_bin['samples_b']))
			for time_bin in time_bins
		],
	}


def all_accumulators(context, block_count):
	'''
	The three sets of sums a run fills: per target and leave-one-out sample,
	per target and block of variants, and per gene of the annotation.
	'''
	bin_count = len(context['time_bins'])
	target_count = len(context['target_names'])
	return {
		'targets': empty_accumulators(bin_count, target_count, max(bin_columns(context['time_bins']))),
		'blocks': empty_accumulators(bin_count, target_count, block_count),
		'genes': empty_gene_accumulators(bin_count, len(context['genes']), ACCUMULATOR_NAMES),
	}


def run_time_series(config):
	'''
	Stream the genotypes once and return the run context together with the
	accumulated sums, from which the estimate, either jackknife and the
	background of genes can be formed.
	'''
	context = build_context(config)
	block_count = config['snp_block_count']
	weights_by_bin = bin_weights(context['time_bins'])
	accumulators = all_accumulators(context, block_count)
	for start, end in chunk_bounds(len(context['variant_index']), config['chunk_size']):
		genotypes = read_genotypes(
			context['source'], context['sample_indices'], context['variant_index'][start:end])
		chunk_blocks = {name: indices[start:end] for name, indices in context['block_indices'].items()}
		chunk_genes = chunk_pairs(context['gene_membership'], start, end)
		accumulate_chunk(
			accumulators, genotypes, context, weights_by_bin, chunk_blocks, chunk_genes, block_count)
	return context, accumulators
