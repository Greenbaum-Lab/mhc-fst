'''
Every gene of the annotation, measured so the genes themselves make a
background.

The genome wide number pools every variant into one value, which is not the
same kind of thing as a gene: it holds a million variants where a gene holds a
few dozen, and it averages over the whole genome rather than over genes. This
module measures each annotated gene on its own, so a focal locus is compared
against a background made of genes like it.

Only the whole sample is carried here, not the leave-one-out samples. Tens of
thousands of genes times a column per individual would not fit in memory, and
the spread across genes is what this background is for.
'''

import gzip
import numpy as np

from gene_regions import parse_attributes


def load_all_genes(annotation_path, gene_biotypes):
	'''
	Name, chromosome, start and end of every gene feature of the annotation
	whose biotype was asked for. The biotype is the annotation's own gene_type,
	so pseudogenes and non-coding genes are left out by naming only the
	biotypes wanted.
	'''
	wanted = set(gene_biotypes)
	genes = []
	opener = gzip.open if annotation_path.endswith('.gz') else open
	with opener(annotation_path, 'rt') as annotation_file:
		for line in annotation_file:
			if line.startswith('#'):
				continue
			fields = line.rstrip('\n').split('\t')
			if len(fields) < 9 or fields[2] != 'gene':
				continue
			attributes = parse_attributes(fields[8])
			if attributes.get('gene_type') not in wanted:
				continue
			genes.append((
				attributes.get('gene_name'),
				fields[0].replace('chr', ''),
				int(fields[3]),
				int(fields[4])))
	if not genes:
		raise ValueError(f'No genes of biotype {sorted(wanted)} in {annotation_path}')
	return genes


def chromosome_slices(chromosome, position):
	'''
	The variants of each chromosome, ordered by position, so a gene resolves
	to a run of them by binary search.
	'''
	slices = {}
	for value in np.unique(chromosome):
		indices = np.where(chromosome == value)[0]
		ordered = indices[np.argsort(position[indices])]
		slices[value] = (ordered, position[ordered])
	return slices


def gene_membership(genes, chromosome, position):
	'''
	Which variants each gene holds, as pairs of variant and gene ordered by
	variant, so a chunk of variants is a slice of the pairs. A variant inside
	several overlapping genes appears once for each of them.
	'''
	slices = chromosome_slices(chromosome, position)
	variant_pairs, gene_pairs = [], []
	for gene_position, (_, gene_chromosome, start, end) in enumerate(genes):
		if gene_chromosome not in slices:
			continue
		ordered, ordered_positions = slices[gene_chromosome]
		low = np.searchsorted(ordered_positions, start, 'left')
		high = np.searchsorted(ordered_positions, end, 'right')
		if high > low:
			variant_pairs.append(ordered[low:high])
			gene_pairs.append(np.full(high - low, gene_position, dtype=np.int32))
	return order_pairs(np.concatenate(variant_pairs), np.concatenate(gene_pairs))


def order_pairs(variant_pairs, gene_pairs):
	order = np.argsort(variant_pairs, kind='stable')
	return variant_pairs[order], gene_pairs[order]


def chunk_pairs(membership, start, end):
	'''
	The pairs of one chunk of variants, with the variant given as its position
	within the chunk.
	'''
	variant_pairs, gene_pairs = membership
	low = np.searchsorted(variant_pairs, start, 'left')
	high = np.searchsorted(variant_pairs, end, 'left')
	return variant_pairs[low:high] - start, gene_pairs[low:high]


def empty_gene_accumulators(bin_count, gene_count, accumulator_names):
	return {name: np.zeros((bin_count, gene_count)) for name in accumulator_names}


def add_genes(gene_accumulators, bin_position, components, chunk_variants, chunk_genes, gene_count):
	'''
	Add a chunk of variants into the per gene sums of the whole sample.
	'''
	component_a, component_b, component_c = components
	for name, values in (
		('sum_a', component_a[chunk_variants, 0]),
		('sum_b', component_b[chunk_variants, 0]),
		('sum_c', component_c[chunk_variants, 0]),
		('variant_count', ~np.isnan(component_a[chunk_variants, 0])),
	):
		weights = np.nan_to_num(values.astype(np.float64))
		gene_accumulators[name][bin_position] += np.bincount(
			chunk_genes, weights=weights, minlength=gene_count)
