'''
Turning the accumulated sums into the results table and the values behind it.

Both jackknives are formed from the same sums: dropping an individual is a
column of the accumulators, dropping a block of variants is the total minus
that block.
'''

import itertools
import numpy as np
import pandas as pd

from fst_core import ratio_of_averages
from jackknife import individual_standard_error, block_standard_error, interval, POINT_ESTIMATE_COLUMN

TABLE_COLUMNS = [
	'polygon_a', 'polygon_b', 'time_start', 'time_end',
	'target', 'phenotype', 'trend', 'time_bp', 'fst',
	'ci_low_individuals', 'ci_high_individuals', 'standard_error_individuals',
	'ci_low_snp_blocks', 'ci_high_snp_blocks', 'standard_error_snp_blocks',
	'n_samples_a', 'n_samples_b', 'n_variants', 'n_blocks',
]
EMPTY_REGION = {'phenotype': '', 'trend': '', 'time_bp': None}


def estimator_values(sums):
	'''
	FST over whatever axis the sums carry.
	'''
	return ratio_of_averages(sums['sum_a'], sums['sum_b'], sums['sum_c'])


def delete_block_sums(accumulators, block_accumulators):
	'''
	Sums with each block of variants removed, as the total of the whole sample
	minus that block.
	'''
	return {
		name: accumulators[name][..., POINT_ESTIMATE_COLUMN][..., None] - block_accumulators[name]
		for name in block_accumulators
	}


def result_row(config, context, accumulators, block_accumulators, values, block_values, index):
	'''
	One row of the results table: the estimate, an interval over individuals,
	an interval over blocks of variants, and the counts they rest on.
	'''
	time_bin = context['time_bins'][index[0]]
	target = context['target_names'][index[1]]
	region = context['region_by_target'].get(target, EMPTY_REGION)
	count_a, count_b = len(time_bin['samples_a']), len(time_bin['samples_b'])
	point_estimate = float(values[index][POINT_ESTIMATE_COLUMN])
	block_counts = block_accumulators['variant_count'][index]
	error_individuals = individual_standard_error(values[index], count_a, count_b)
	error_blocks = block_standard_error(block_values[index], block_counts)
	individuals = interval(point_estimate, error_individuals, config['confidence_level'])
	blocks = interval(point_estimate, error_blocks, config['confidence_level'])
	return {
		'polygon_a': config['polygon_a'],
		'polygon_b': config['polygon_b'],
		'time_start': time_bin['time_start'],
		'time_end': time_bin['time_end'],
		'target': target,
		'phenotype': region['phenotype'],
		'trend': region['trend'],
		'time_bp': region['time_bp'],
		'fst': point_estimate,
		'ci_low_individuals': individuals[0],
		'ci_high_individuals': individuals[1],
		'standard_error_individuals': error_individuals,
		'ci_low_snp_blocks': blocks[0],
		'ci_high_snp_blocks': blocks[1],
		'standard_error_snp_blocks': error_blocks,
		'n_samples_a': count_a,
		'n_samples_b': count_b,
		'n_variants': int(accumulators['variant_count'][index][POINT_ESTIMATE_COLUMN]),
		'n_blocks': int(np.count_nonzero(block_counts)),
	}


def build_table(config, context, all_sums):
	'''
	The results table, one row per time bin and target.
	'''
	accumulators, block_accumulators = all_sums['targets'], all_sums['blocks']
	values = estimator_values(accumulators)
	block_values = estimator_values(delete_block_sums(accumulators, block_accumulators))
	positions = itertools.product(range(len(context['time_bins'])), range(len(context['target_names'])))
	rows = [
		result_row(config, context, accumulators, block_accumulators, values, block_values, index)
		for index in positions
	]
	return pd.DataFrame(rows, columns=TABLE_COLUMNS)


def save_jackknife_values(output_path, context, all_sums):
	'''
	The value of every leave-one-out sample, kept so later comparisons between
	a region and the background do not need the genotypes again. Column zero
	of the individual array is the estimate, then one column per individual of
	the first population and one per individual of the second.
	'''
	accumulators, block_accumulators = all_sums['targets'], all_sums['blocks']
	np.savez(
		output_path,
		time_start=np.array([time_bin['time_start'] for time_bin in context['time_bins']]),
		time_end=np.array([time_bin['time_end'] for time_bin in context['time_bins']]),
		samples_a=np.array([len(time_bin['samples_a']) for time_bin in context['time_bins']]),
		samples_b=np.array([len(time_bin['samples_b']) for time_bin in context['time_bins']]),
		targets=np.array(context['target_names']),
		variant_count=accumulators['variant_count'],
		block_variant_count=block_accumulators['variant_count'],
		individuals=estimator_values(accumulators),
		snp_blocks=estimator_values(delete_block_sums(accumulators, block_accumulators)))


def save_regions(output_path, regions):
	pd.DataFrame(regions).to_csv(output_path, index=False)
