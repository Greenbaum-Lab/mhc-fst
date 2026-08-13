'''
FST of focal regions and of the genome across a series of time bins, for two
populations, with a bootstrap over individuals.

The genotypes are read twice. The first pass measures the call rate of every
variant in every bin, which is what the SNP set alternatives are built from.
The second pass streams the variants again and accumulates, for each bin, the
sums that both multi-locus estimators need, separately for every region and
every SNP set alternative. Focal regions and the genome wide background share
one set of leave-one-out samples per bin, so their intervals are comparable.
'''

import numpy as np

from fst_core import weir_cockerham_components, per_variant_fst
from jackknife import paired_jackknife_weights, column_count
from gene_regions import load_gene_spans, build_regions
from time_populations import load_populations, overlapping_time_bins
from variant_masks import region_masks, filter_masks, FILTER_MODES, GENOME_WIDE_TARGET
from genotype_source import (
	open_genotypes,
	autosomal_variants,
	select_sample_indices,
	read_block,
	allele_statistics,
	called_fraction,
)

ACCUMULATOR_NAMES = ['sum_a', 'sum_b', 'sum_c', 'sum_fst', 'variant_count']


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


def block_rows(row_by_sample_id, sample_ids):
	return np.array([row_by_sample_id[sample_id] for sample_id in sample_ids])


def bin_call_rates(source, variant_index, sample_indices, bin_rows, chunk_size):
	'''
	Lowest call rate of the two populations of a bin at each variant, one row
	per bin. The call rate threshold is applied to this, so a variant is kept
	only when both populations of the bin carry enough calls.
	'''
	call_rates = np.empty((len(bin_rows), len(variant_index)), dtype=np.float32)
	for start, end in chunk_bounds(len(variant_index), chunk_size):
		block = read_block(source, sample_indices, variant_index[start:end])
		for bin_position, (rows_a, rows_b) in enumerate(bin_rows):
			rate_a = called_fraction(block[rows_a])
			rate_b = called_fraction(block[rows_b])
			call_rates[bin_position, start:end] = np.minimum(rate_a, rate_b)
	return call_rates


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


def empty_accumulators(bin_count, target_count, widest_bin):
	shape = (bin_count, target_count, len(FILTER_MODES), widest_bin)
	return {name: np.zeros(shape) for name in ACCUMULATOR_NAMES}


def add_selection(accumulators, index, components, variant_fst, selected):
	'''
	Add the variants of one region and SNP set alternative into the running
	sums, dropping the variants the estimator marked unusable. A bin fills
	only as many columns as it has individuals, so each sum is written into
	the leading columns of its row.
	'''
	component_a, component_b, component_c = components
	columns = component_a.shape[1]
	accumulators['sum_a'][index][:columns] += np.nansum(component_a[selected], axis=0, dtype=np.float64)
	accumulators['sum_b'][index][:columns] += np.nansum(component_b[selected], axis=0, dtype=np.float64)
	accumulators['sum_c'][index][:columns] += np.nansum(component_c[selected], axis=0, dtype=np.float64)
	accumulators['sum_fst'][index][:columns] += np.nansum(variant_fst[selected], axis=0, dtype=np.float64)
	accumulators['variant_count'][index][:columns] += np.count_nonzero(~np.isnan(variant_fst[selected]), axis=0)


def accumulate_bin(accumulators, bin_position, components, variant_fst, target_masks, target_names, bin_filter_masks):
	for target_position, target_name in enumerate(target_names):
		for mode_position, mode_name in enumerate(FILTER_MODES):
			selected = target_masks[target_name] & bin_filter_masks[mode_name]
			if selected.any():
				index = (bin_position, target_position, mode_position)
				add_selection(accumulators, index, components, variant_fst, selected)


def accumulate_chunk(accumulators, block, bin_rows, weights_by_bin, target_masks, target_names, filter_masks_by_mode):
	for bin_position, (rows_a, rows_b) in enumerate(bin_rows):
		weights_a, weights_b = weights_by_bin[bin_position]
		count_a, number_a, het_a = allele_statistics(block[rows_a], weights_a)
		count_b, number_b, het_b = allele_statistics(block[rows_b], weights_b)
		components = weir_cockerham_components(number_a, number_b, count_a, count_b, het_a, het_b)
		variant_fst = per_variant_fst(*components)
		bin_filter_masks = {mode: masks[bin_position] for mode, masks in filter_masks_by_mode.items()}
		accumulate_bin(accumulators, bin_position, components, variant_fst, target_masks, target_names, bin_filter_masks)


def build_context(config):
	'''
	Everything a run is defined by: the time bins compared, the focal regions,
	the variants and the samples they are read from.
	'''
	populations = load_populations(config['populations_path'])
	time_bins = overlapping_time_bins(populations, config['polygon_a'], config['polygon_b'], config['genotype_source'])
	regions = build_regions(load_gene_spans(config['annotation_path'], config['genes']), config['flank_sizes'])
	source = open_genotypes(config['bed_prefix'], config['threads'])
	variant_index, chromosome, position = autosomal_variants(source)
	sample_ids, row_by_sample_id = union_sample_ids(time_bins)
	return {
		'source': source,
		'time_bins': time_bins,
		'regions': regions,
		'variant_index': variant_index,
		'target_masks': region_masks(regions, chromosome, position),
		'target_names': [region['region_id'] for region in regions] + [GENOME_WIDE_TARGET],
		'gene_by_target': {region['region_id']: region['gene'] for region in regions},
		'sample_indices': select_sample_indices(source, sample_ids),
		'bin_rows': [
			(block_rows(row_by_sample_id, time_bin['samples_a']), block_rows(row_by_sample_id, time_bin['samples_b']))
			for time_bin in time_bins
		],
	}


def run_time_series(config):
	'''
	Run both passes and return the run context together with the accumulated
	sums, from which either estimator can be formed.
	'''
	context = build_context(config)
	call_rates = bin_call_rates(
		context['source'], context['variant_index'], context['sample_indices'],
		context['bin_rows'], config['chunk_size'])
	filter_masks_by_mode = filter_masks(call_rates, config['call_rate_threshold'])
	weights_by_bin = bin_weights(context['time_bins'])
	accumulators = empty_accumulators(
		len(context['time_bins']), len(context['target_names']), max(bin_columns(context['time_bins'])))
	for start, end in chunk_bounds(len(context['variant_index']), config['chunk_size']):
		block = read_block(context['source'], context['sample_indices'], context['variant_index'][start:end])
		chunk_targets = {name: mask[start:end] for name, mask in context['target_masks'].items()}
		chunk_filters = {mode: masks[:, start:end] for mode, masks in filter_masks_by_mode.items()}
		accumulate_chunk(
			accumulators, block, context['bin_rows'], weights_by_bin,
			chunk_targets, context['target_names'], chunk_filters)
	return context, accumulators
