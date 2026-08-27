'''
Draws the Manhattan figures of a temporal FST scan table.
'''

import argparse
import pathlib

from temporal_scan.scan_results import prepare, significant_in_every_group
from temporal_scan.manhattan import plot_group, plot_shared

MIN_VARIANTS = 5


def parse_arguments():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('results', help='TSV written by the temporal FST scan')
	parser.add_argument('--out', default='.', help='Directory the figures are written to')
	parser.add_argument('--min-variants', type=int, default=MIN_VARIANTS, help='Keep genes holding more than this many variants')
	return parser.parse_args()


def file_name(group):
	return group.lower().replace(' ', '_')


def main():
	arguments = parse_arguments()
	genes, centres = prepare(arguments.results, arguments.min_variants)
	directory = pathlib.Path(arguments.out)
	directory.mkdir(parents=True, exist_ok=True)
	for group in sorted(genes['group'].unique()):
		block = genes[genes['group'] == group]
		plot_group(genes, centres, group, directory / f'manhattan_{file_name(group)}.png')
		print(f'{group}: {len(block)} genes, {int(block["outside_null"].sum())} outside the null')
	shared = significant_in_every_group(genes)
	plot_shared(shared, centres, directory / 'manhattan_shared.png')
	print(f'shared: {shared["gene"].nunique()} genes outside the null in every group')


if __name__ == '__main__':
	main()
