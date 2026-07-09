import argparse
import pandas as pd
import matplotlib.pyplot as plt

from manhattan import (
	load_pair_labels,
	load_chromosome_windows,
	load_genomewide,
	assign_cumulative_positions,
	apply_ratio,
)

CANDIDATE_COLORS = {
	'filtered': '#e0e0e0',
	'background': '#b0b0b0',
	'valley': '#2980b9',
	'peak': '#c0392b',
}


def load_outliers(outliers_path):
	'''
	Load the outlier_scan.py output table for one bin.
	'''
	return pd.read_csv(outliers_path)


def merge_candidates(windows, outliers):
	'''
	Attach each window's outlier-scan candidate label. Windows absent from
	the outlier table were dropped by its min-SNP filter, and are labelled
	'filtered' rather than 'background' so excluded windows stay visually
	distinct from windows that were tested and found non-significant.
	'''
	merged = windows.merge(
		outliers[['chromosome', 'position', 'candidate']],
		on=['chromosome', 'position'], how='left')
	merged['candidate'] = merged['candidate'].fillna('filtered')
	return merged


def plot_manhattan_with_candidates(windows, pair_label, output_path):
	'''
	Draw the FST-ratio Manhattan plot on a symlog y-axis, colored by
	outlier-scan candidate status. Symlog (linear near zero, log further
	out) is used instead of log2 because WC/Hudson FST windows can be
	negative or zero, which a plain log2 would turn into NaNs and silently
	drop from the plot.
	'''
	figure, axis = plt.subplots(figsize=(14, 5))
	for candidate, color in CANDIDATE_COLORS.items():
		subset = windows[windows['candidate'] == candidate]
		axis.scatter(subset['cumulative_position'], subset['value'], s=8, color=color, label=candidate)
	axis.axhline(1.0, color='black', linewidth=0.8, linestyle='--')
	axis.set_yscale('symlog', linthresh=0.1)
	ticks = windows.groupby('chromosome')['cumulative_position'].mean().sort_index()
	axis.set_xticks(ticks.to_numpy())
	axis.set_xticklabels(ticks.index)
	axis.set_xlabel('Chromosome')
	axis.set_ylabel('FST ratio (window / genome-wide), symlog scale')
	axis.set_title(pair_label)
	axis.legend(loc='upper right', markerscale=2, frameon=False)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)
	figure.tight_layout()
	figure.savefig(output_path, dpi=150)


def build_outlier_manhattan(fst_dir, bin_label, chromosomes, window_size, pair_index, stat_index, outliers_path, output_path):
	'''
	Assemble the ratio Manhattan plot for one bin, exactly as manhattan.py
	does, then color each window by its outlier_scan.py candidate status.
	'''
	window_prefix = f'{fst_dir}/{bin_label}/fst_{bin_label}'
	genomewide_path = f'{fst_dir}/gw/gw_{bin_label}.npy'
	pair_labels = load_pair_labels(window_prefix)
	windows_by_chromosome = [
		load_chromosome_windows(window_prefix, chromosome, window_size, pair_index, stat_index)
		for chromosome in chromosomes]
	windows = assign_cumulative_positions(windows_by_chromosome)
	genomewide_value = load_genomewide(genomewide_path, pair_index, stat_index)
	windows = apply_ratio(windows, genomewide_value)
	outliers = load_outliers(outliers_path)
	windows = merge_candidates(windows, outliers)
	plot_manhattan_with_candidates(windows, pair_labels[pair_index], output_path)


def main():
	parser = argparse.ArgumentParser(description='Manhattan plot of FST ratio colored by outlier-scan candidate status')
	parser.add_argument('--fst-dir', required=True)
	parser.add_argument('--bin', required=True)
	parser.add_argument('--chromosomes', nargs='+', type=int, required=True)
	parser.add_argument('--window-size', type=int, required=True)
	parser.add_argument('--pair-index', type=int, default=0)
	parser.add_argument('--stat-index', type=int, default=0)
	parser.add_argument('--outliers-path', required=True)
	parser.add_argument('--output-path', required=True)
	args = parser.parse_args()
	build_outlier_manhattan(args.fst_dir, args.bin, args.chromosomes, args.window_size,
	                        args.pair_index, args.stat_index, args.outliers_path, args.output_path)


if __name__ == '__main__':
	main()
