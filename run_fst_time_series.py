import json
import pathlib
import argparse

from fst_time_series import run_time_series
from fst_results import build_table, save_jackknife_values, save_regions
from gene_background import gene_table, background_table, save_gene_sums


def load_config(config_path):
	with open(config_path) as config_file:
		return json.load(config_file)


def write_outputs(config, context, all_sums):
	'''
	Write the results table, the leave-one-out values, the focal region
	coordinates the run resolved, and every gene of the annotation with the
	background it makes.
	'''
	output_dir = pathlib.Path(config['output_dir'])
	output_dir.mkdir(parents=True, exist_ok=True)
	table = build_table(config, context, all_sums)
	table.to_csv(output_dir / 'fst_time_series.csv', index=False)
	save_jackknife_values(output_dir / 'fst_jackknife.npz', context, all_sums)
	save_regions(output_dir / 'focal_regions.csv', context['regions'])
	gene_table(context, all_sums['genes']).to_csv(output_dir / 'fst_per_gene.csv', index=False)
	save_gene_sums(output_dir / 'fst_per_gene.npz', context, all_sums['genes'])
	background_table(
		context, all_sums['genes'], config['min_gene_variants'], config['confidence_level']
	).to_csv(output_dir / 'gene_background.csv', index=False)
	return table


def main():
	parser = argparse.ArgumentParser(description='FST across time bins for two polygon populations')
	parser.add_argument('--config', required=True)
	args = parser.parse_args()
	config = load_config(args.config)
	context, all_sums = run_time_series(config)
	table = write_outputs(config, context, all_sums)
	print(f'{len(context["time_bins"])} time bins, {len(context["target_names"])} targets, '
	      f'{len(context["genes"])} annotation genes, {len(table)} rows')


if __name__ == '__main__':
	main()
