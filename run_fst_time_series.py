import json
import pathlib
import argparse

from fst_time_series import run_time_series
from fst_results import build_table, save_jackknife_values, save_regions


def load_config(config_path):
	with open(config_path) as config_file:
		return json.load(config_file)


def write_outputs(config, context, accumulators):
	'''
	Write the results table, the leave-one-out values and the focal region
	coordinates the run resolved.
	'''
	output_dir = pathlib.Path(config['output_dir'])
	output_dir.mkdir(parents=True, exist_ok=True)
	table = build_table(config, context, accumulators)
	table.to_csv(output_dir / 'fst_time_series.csv', index=False)
	save_jackknife_values(output_dir / 'fst_jackknife.npz', context, accumulators)
	save_regions(output_dir / 'focal_regions.csv', context['regions'])
	return table


def main():
	parser = argparse.ArgumentParser(description='FST across time bins for two polygon populations')
	parser.add_argument('--config', required=True)
	args = parser.parse_args()
	config = load_config(args.config)
	context, accumulators = run_time_series(config)
	table = write_outputs(config, context, accumulators)
	print(f'{len(context["time_bins"])} time bins, {len(context["target_names"])} targets, {len(table)} rows')


if __name__ == '__main__':
	main()
