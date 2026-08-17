'''
Which variants enter which estimate, and how they are blocked.

Every variant is used: there is no call rate threshold, and the only variants
dropped are those the Weir & Cockerham code itself rejects for holding two or
fewer called alleles.
'''

import numpy as np

GENOME_WIDE_TARGET = 'genome_wide'


def region_masks(regions, chromosome, position):
	'''
	Variants inside each focal region, plus the genome wide target holding
	every autosomal variant.
	'''
	masks = {
		region['locus']: (chromosome == region['chromosome'])
		& (position >= region['start'])
		& (position <= region['end'])
		for region in regions
	}
	masks[GENOME_WIDE_TARGET] = np.ones(len(position), dtype=bool)
	return masks


def block_indices(target_masks, block_count):
	'''
	The block each variant belongs to within a target, and -1 for variants the
	target does not hold. Blocks are runs of neighbouring variants holding
	equal counts, so deleting one in turn needs no block weighting. A target
	with fewer variants than blocks leaves some blocks empty.
	'''
	indices = {}
	for target_name, mask in target_masks.items():
		total = int(mask.sum())
		assigned = np.full(len(mask), -1, dtype=np.int16)
		if total > 0:
			ranks = np.arange(total)
			assigned[mask] = np.minimum(ranks * block_count // total, block_count - 1)
		indices[target_name] = assigned
	return indices
