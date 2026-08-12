'''
Weir & Cockerham 1984 FST components and the two multi-locus estimators.

The component computation is the body of delphi/analyses/fst.py::_compute_fst,
changed only to return a, b and c per variant instead of the ratio a / (a+b+c),
because a ratio of averages needs the components. check_fst_parity.py asserts
that the per-variant ratio rebuilt from these components equals the delphi
implementation exactly.
'''

import numpy as np


def weir_cockerham_components(an1, an2, ac1, ac2, het1, het2):
	'''
	Per-variant a, b and c of the Weir & Cockerham 1984 FST decomposition.
	Input arrays hold one row per variant and one column per weighting of the
	individuals. Variants with two or fewer called alleles in either
	population, or with a zero denominator, are nan in all three components.
	'''
	with np.errstate(invalid='ignore', divide='ignore'):
		valid = (an1 > 2.0) & (an2 > 2.0)
		p1      = ac1 / an1
		p2      = ac2 / an2
		n1      = an1 / 2.0
		n2      = an2 / 2.0
		n_total = n1 + n2
		n_bar   = n_total / 2.0
		n_c     = n_total - (n1 * n1 + n2 * n2) / n_total
		ac_total = ac1 + ac2
		an_total = an1 + an2
		p_bar    = ac_total / an_total
		s2       = (n1 * (p1 - p_bar) ** 2 + n2 * (p2 - p_bar) ** 2) / n_bar
		h_bar    = (het1 + het2) / n_total
		a = (n_bar / n_c) * (s2 - h_bar / (4.0 * n_bar))
		b = (n_bar / (n_bar - 1.0)) * (
			p_bar * (1.0 - p_bar) - 0.5 * s2 - (2.0 * n_bar - 1.0) * h_bar / (4.0 * n_bar)
		)
		c = h_bar / 2.0
		keep = valid & ((a + b + c) != 0.0)
	return np.where(keep, a, np.nan), np.where(keep, b, np.nan), np.where(keep, c, np.nan)


def per_variant_fst(a, b, c):
	'''
	Per-variant FST clipped at zero, the value delphi/analyses/fst.py returns.
	'''
	with np.errstate(invalid='ignore', divide='ignore'):
		return np.maximum(0.0, a / (a + b + c))


def ratio_of_averages(sum_a, sum_b, sum_c):
	'''
	Multi-locus FST as sum(a) / sum(a + b + c) over the variants of a region,
	the standard Weir & Cockerham combination, clipped at zero once at the end.
	'''
	denominator = sum_a + sum_b + sum_c
	with np.errstate(invalid='ignore', divide='ignore'):
		return np.where(denominator != 0.0, np.maximum(0.0, sum_a / denominator), np.nan)


def average_of_ratios(sum_per_variant_fst, variant_count):
	'''
	Multi-locus FST as the mean of the per-variant values, the combination
	delphi/analyses/fst.py uses, where each variant is clipped at zero first.
	'''
	with np.errstate(invalid='ignore', divide='ignore'):
		return np.where(variant_count > 0, sum_per_variant_fst / variant_count, np.nan)
