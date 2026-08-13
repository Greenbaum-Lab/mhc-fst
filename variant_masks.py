'''
Which variants enter which estimate.

Two independent selections are applied together: the target selects a genomic
region, and the SNP set alternative selects variants by call rate. The three
alternatives are computed side by side so the effect of filtering is visible
rather than assumed.
'''

import numpy as np

GENOME_WIDE_TARGET = 'genome_wide'
FILTER_MODES = ['none', 'per_bin', 'intersection']


def region_masks(regions, chromosome, position):
	'''
	Variants inside each focal region, plus the genome wide target holding
	every autosomal variant.
	'''
	masks = {
		region['region_id']: (chromosome == region['chromosome'])
		& (position >= region['start'])
		& (position <= region['end'])
		for region in regions
	}
	masks[GENOME_WIDE_TARGET] = np.ones(len(position), dtype=bool)
	return masks


def filter_masks(minimum_call_rates, call_rate_threshold):
	'''
	Variants kept per time bin under each SNP set alternative, one row per bin.

	none         every variant, no call rate threshold.
	per_bin      variants above the threshold in that bin, so the variant set
	             changes between bins with coverage.
	intersection variants above the threshold in every bin, so all bins share
	             one fixed variant set.
	'''
	passing = minimum_call_rates >= call_rate_threshold
	shared = passing.all(axis=0)
	return {
		'none': np.ones(passing.shape, dtype=bool),
		'per_bin': passing,
		'intersection': np.broadcast_to(shared, passing.shape),
	}
