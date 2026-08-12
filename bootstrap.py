'''
Bootstrap over individuals.

Resampling individuals with replacement is a reweighting of the individuals,
so a replicate is one column of weights and every replicate of a population is
computed in a single matrix product. Column zero weights each individual once
and is therefore the point estimate, not a replicate.
'''

import numpy as np

POINT_ESTIMATE_COLUMN = 0


def resampling_weights(sample_count, replicate_count, random_generator):
	'''
	Weights of the point estimate and of each bootstrap replicate, one column
	each. Replicate columns are multinomial counts of the individuals drawn
	with replacement, keeping the resampled size equal to the original.
	'''
	probabilities = np.full(sample_count, 1.0 / sample_count)
	replicates = random_generator.multinomial(sample_count, probabilities, size=replicate_count)
	point_estimate = np.ones((1, sample_count))
	return np.vstack([point_estimate, replicates]).T.astype(np.float32)


def percentile_interval(replicate_values, confidence_level):
	'''
	Percentile bootstrap interval of the replicate values, ignoring replicates
	that produced no usable variant.
	'''
	if np.all(np.isnan(replicate_values)):
		return np.nan, np.nan
	tail = (1.0 - confidence_level) / 2.0 * 100.0
	low, high = np.nanpercentile(replicate_values, [tail, 100.0 - tail])
	return float(low), float(high)


def basic_interval(point_estimate, replicate_values, confidence_level):
	'''
	Basic bootstrap interval, the percentile interval reflected about the
	point estimate.

	Resampling individuals with replacement duplicates individuals, which
	inflates the allele frequency variance of a replicate beyond that of a
	genuine sample of the same size. The Weir & Cockerham sample size
	correction is computed from the resampled size and so does not remove it,
	leaving every replicate biased upwards by roughly the reciprocal of the
	smaller population size. Reflecting about the point estimate cancels that
	shift. The bounds are not clipped at zero, so an interval reaching below
	zero shows a value that cannot be told apart from no differentiation.
	'''
	low, high = percentile_interval(replicate_values, confidence_level)
	return 2.0 * point_estimate - high, 2.0 * point_estimate - low
