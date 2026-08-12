'''
Time bins of the polygon populations defined in populations.json.

A polygon population holds one entry per time bin, each carrying the sample
ids inside the polygon ring dated to that bin. Two polygons are compared only
over the bins both of them define.
'''

import json


def load_populations(populations_path):
	with open(populations_path) as populations_file:
		return json.load(populations_file)


def polygon_time_bins(populations, polygon_name, genotype_source):
	'''
	Sample ids of one polygon population per time bin, keyed by the bin edges
	in years before present.
	'''
	time_bins = {}
	for population in populations:
		definition = population['definition']
		if population.get('source') != 'polygon' or definition.get('polygon') != polygon_name:
			continue
		edges = (definition['time_start'], definition['time_end'])
		time_bins[edges] = population['samples'][genotype_source]
	if not time_bins:
		raise KeyError(f'No polygon population named {polygon_name}')
	return time_bins


def overlapping_time_bins(populations, polygon_a, polygon_b, genotype_source):
	'''
	The time bins both polygons define, ordered from recent to ancient, each
	carrying the sample ids of the two populations compared in that bin.
	'''
	bins_a = polygon_time_bins(populations, polygon_a, genotype_source)
	bins_b = polygon_time_bins(populations, polygon_b, genotype_source)
	shared_edges = sorted(set(bins_a) & set(bins_b))
	return [
		{
			'time_start': time_start,
			'time_end': time_end,
			'samples_a': bins_a[(time_start, time_end)],
			'samples_b': bins_b[(time_start, time_end)],
		}
		for time_start, time_end in shared_edges
	]
