import json
import argparse
from assign_populations import assign_populations
from subset_bed import write_subset_bed


def load_config(config_path):
	'''
	Load the region bounding boxes and time bins that define the
	spatiotemporal populations.
	'''
	with open(config_path) as config_file:
		config = json.load(config_file)
	return config['region_boxes'], config['time_bins']


def build_subset_for_bin(janno_path, source_bed_prefix, region_boxes, time_bin, subset_prefix):
	'''
	Assign samples for a single time bin and write a subset BED whose family
	ID is the population label, yielding one population per region in the bin.
	'''
	sample_populations = assign_populations(janno_path, region_boxes, [time_bin])
	write_subset_bed(source_bed_prefix, sample_populations, subset_prefix)


def build_all_subsets(janno_path, source_bed_prefix, config_path, subset_prefix_template):
	'''
	Write one subset BED per time bin, naming each output by its bin edges
	via the subset prefix template.
	'''
	region_boxes, time_bins = load_config(config_path)
	for time_bin in time_bins:
		low, high = time_bin
		subset_prefix = subset_prefix_template.format(low=low, high=high)
		build_subset_for_bin(janno_path, source_bed_prefix, region_boxes, time_bin, subset_prefix)


def main():
	parser = argparse.ArgumentParser(description='Build one subset BED per time bin')
	parser.add_argument('--janno', required=True)
	parser.add_argument('--source-bed-prefix', required=True)
	parser.add_argument('--config', required=True)
	parser.add_argument('--subset-prefix-template', required=True)
	args = parser.parse_args()
	build_all_subsets(args.janno, args.source_bed_prefix, args.config, args.subset_prefix_template)


if __name__ == '__main__':
	main()
