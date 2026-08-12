'''
Turning the accumulated sums into the results table and the raw replicates.

Both multi-locus estimators are formed from the same sums, so they differ only
in how the variants of a region are combined, never in the underlying counts.
'''

import itertools
import numpy as np
import pandas as pd

from fst_core import ratio_of_averages, average_of_ratios
from bootstrap import percentile_interval, basic_interval, POINT_ESTIMATE_COLUMN
from variant_masks import FILTER_MODES

TABLE_COLUMNS = [
	'polygon_a', 'polygon_b', 'time_start', 'time_end', 'target', 'filter_mode',
	'estimator', 'fst', 'ci_low_basic', 'ci_high_basic',
	'ci_low_percentile', 'ci_high_percentile', 'bootstrap_mean',
	'n_samples_a', 'n_samples_b', 'n_variants',
]


def estimator_values(accumulators):
	'''
	FST of every bin, target, SNP set alternative and weighting column, under
	each estimator.
	'''
	return {
		'ratio_of_averages': ratio_of_averages(
			accumulators['sum_a'], accumulators['sum_b'], accumulators['sum_c']),
		'average_of_ratios': average_of_ratios(
			accumulators['sum_fst'], accumulators['variant_count']),
	}


def result_row(config, context, accumulators, estimator_name, values, index):
	'''
	One row of the results table: the point estimate, both forms of its
	bootstrap interval and the counts the estimate rests on.
	'''
	bin_position, target_position, mode_position = index
	time_bin = context['time_bins'][bin_position]
	series = values[index]
	point_estimate = float(series[POINT_ESTIMATE_COLUMN])
	replicates = series[POINT_ESTIMATE_COLUMN + 1:]
	ci_low_basic, ci_high_basic = basic_interval(point_estimate, replicates, config['confidence_level'])
	ci_low_percentile, ci_high_percentile = percentile_interval(replicates, config['confidence_level'])
	return {
		'polygon_a': config['polygon_a'],
		'polygon_b': config['polygon_b'],
		'time_start': time_bin['time_start'],
		'time_end': time_bin['time_end'],
		'target': context['target_names'][target_position],
		'filter_mode': FILTER_MODES[mode_position],
		'estimator': estimator_name,
		'fst': point_estimate,
		'ci_low_basic': ci_low_basic,
		'ci_high_basic': ci_high_basic,
		'ci_low_percentile': ci_low_percentile,
		'ci_high_percentile': ci_high_percentile,
		'bootstrap_mean': float(np.nanmean(replicates)) if not np.all(np.isnan(replicates)) else np.nan,
		'n_samples_a': len(time_bin['samples_a']),
		'n_samples_b': len(time_bin['samples_b']),
		'n_variants': int(accumulators['variant_count'][index][POINT_ESTIMATE_COLUMN]),
	}


def build_table(config, context, accumulators):
	'''
	The results table, one row per bin, target, SNP set alternative and
	estimator.
	'''
	values_by_estimator = estimator_values(accumulators)
	positions = itertools.product(
		range(len(context['time_bins'])),
		range(len(context['target_names'])),
		range(len(FILTER_MODES)))
	rows = [
		result_row(config, context, accumulators, estimator_name, values, index)
		for index in positions
		for estimator_name, values in values_by_estimator.items()
	]
	return pd.DataFrame(rows, columns=TABLE_COLUMNS)


def save_replicates(output_path, context, accumulators):
	'''
	The value of every bootstrap replicate, kept so later comparisons between
	a region and the background do not need the genotypes again. Column zero
	of the value arrays is the point estimate, the rest are replicates.
	'''
	values_by_estimator = estimator_values(accumulators)
	np.savez(
		output_path,
		time_start=np.array([time_bin['time_start'] for time_bin in context['time_bins']]),
		time_end=np.array([time_bin['time_end'] for time_bin in context['time_bins']]),
		targets=np.array(context['target_names']),
		filter_modes=np.array(FILTER_MODES),
		variant_count=accumulators['variant_count'],
		**values_by_estimator)


def save_regions(output_path, regions):
	pd.DataFrame(regions).to_csv(output_path, index=False)
