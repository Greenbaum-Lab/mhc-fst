import numpy as np
from bed_reader import open_bed, to_bed


def select_sample_indices(source_iid, sample_populations):
	'''
	Return the source row indices for the requested samples and the
	population label aligned to each index, preserving source order.
	Raises KeyError if any requested sample is absent from the source.
	'''
	population_by_id = dict(zip(sample_populations['Poseidon_ID'], sample_populations['population']))
	requested = set(population_by_id)
	present_mask = np.isin(source_iid, list(requested))
	sample_indices = np.where(present_mask)[0]
	missing = requested - set(source_iid[sample_indices])
	if missing:
		raise KeyError(f'Samples not found in source BED: {sorted(missing)}')
	populations = np.array([population_by_id[iid] for iid in source_iid[sample_indices]])
	return sample_indices, populations


def write_subset_bed(source_bed_prefix, sample_populations, output_bed_prefix):
	'''
	Write a PLINK BED subset containing only the requested samples, with
	the family ID set to each sample's population label and all variants
	preserved. Genotypes and variant metadata are copied unchanged.
	'''
	source = open_bed(f'{source_bed_prefix}.bed', num_threads=8)
	sample_indices, populations = select_sample_indices(source.iid, sample_populations)
	genotypes = source.read(index=np.s_[sample_indices, :], dtype='int8')
	properties = {
		'fid': populations,
		'iid': source.iid[sample_indices],
		'chromosome': source.chromosome,
		'sid': source.sid,
		'bp_position': source.bp_position,
		'allele_1': source.allele_1,
		'allele_2': source.allele_2,
	}
	to_bed(f'{output_bed_prefix}.bed', genotypes, properties=properties, num_threads=8)
