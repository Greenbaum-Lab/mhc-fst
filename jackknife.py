'''
Jackknife over individuals and over blocks of variants.

Both answer a different question. Dropping individuals asks whether a
different set of people would give a different answer. Dropping a block of
neighbouring variants asks whether a different stretch of genome would. For a
genome wide average the first is small and the second is what matters, and for
a region holding few variants the second dominates by far.

Leaving one out is a reweighting, so a whole set of leave-one-out samples is
one matrix product. A leave-one-out sample never holds the same individual
twice, unlike resampling with replacement, so it stays a genuine sample of its
size and the Weir & Cockerham correction remains valid. Nothing here is
random, so a run is reproducible without a seed.
'''

import numpy as np
from statistics import NormalDist

POINT_ESTIMATE_COLUMN = 0


def column_count(count_a, count_b):
	'''
	Columns a time bin needs: the estimate itself and one per individual.
	'''
	return 1 + count_a + count_b


def paired_jackknife_weights(count_a, count_b):
	'''
	Weights of both populations of a time bin, aligned column by column.
	Column zero holds every individual of both populations. The next count_a
	columns drop one individual of the first population and keep the second
	whole, and the remaining count_b columns do the reverse.
	'''
	columns = column_count(count_a, count_b)
	weights_a = np.ones((count_a, columns), dtype=np.float32)
	weights_b = np.ones((count_b, columns), dtype=np.float32)
	weights_a[:, 1:1 + count_a] = 1.0 - np.eye(count_a, dtype=np.float32)
	weights_b[:, 1 + count_a:] = 1.0 - np.eye(count_b, dtype=np.float32)
	return weights_a, weights_b


def group_variance(leave_one_out_values):
	'''
	Variance a group contributes, from the spread of the estimates that drop
	each of its members in turn.
	'''
	count = len(leave_one_out_values)
	if count < 2 or np.all(np.isnan(leave_one_out_values)):
		return np.nan
	deviations = leave_one_out_values - np.nanmean(leave_one_out_values)
	return (count - 1) / count * np.nansum(deviations ** 2)


def individual_standard_error(values, count_a, count_b):
	'''
	Standard error over individuals, summing what each population contributes.
	'''
	values_a = values[1:1 + count_a]
	values_b = values[1 + count_a:1 + count_a + count_b]
	return float(np.sqrt(group_variance(values_a) + group_variance(values_b)))


def block_standard_error(delete_block_values, block_variant_count):
	'''
	Standard error over blocks of variants, counting only the blocks that hold
	variants. Empty blocks would return the estimate unchanged and shrink the
	spread.
	'''
	used = delete_block_values[block_variant_count > 0]
	return float(np.sqrt(group_variance(used)))


def interval(point_estimate, standard_error, confidence_level):
	'''
	Interval around the estimate. The bounds are not clipped at zero, so an
	interval reaching below zero shows a value that cannot be told apart from
	no differentiation.
	'''
	quantile = NormalDist().inv_cdf(1.0 - (1.0 - confidence_level) / 2.0)
	return point_estimate - quantile * standard_error, point_estimate + quantile * standard_error
