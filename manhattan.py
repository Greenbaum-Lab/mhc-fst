import argparse
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load and display the .npy file
# data = np.load('/home/lab2/cluster/fst/10k_window/1000_2000/fst_1000_2000_fst_1.npy')
# data = np.load('/home/lab2/cluster/fst/10k_window/gw/gw_1000_2000.npy')
# print("Shape:", data.shape)
# print("Data type:", data.dtype)
# print("Size:", data.size)
# print("\nData:")
# print(data)


def load_pair_labels(window_prefix):
	'''
	Read the ordered population-pair labels written by the FST engine.
	'''
	pairs_path = pathlib.Path(f'{window_prefix}_pairs.txt')
	with pairs_path.open() as pairs_file:
		return [line.strip() for line in pairs_file]


def load_chromosome_windows(window_prefix, chromosome, window_size, pair_index, stat_index):
	'''
	Load per-window FST for one chromosome and one population pair, returning
	a DataFrame with the window midpoint position and the statistic value.
	'''
	data = np.load(pathlib.Path(f'{window_prefix}_fst_{chromosome}.npy'))
	values = data[:, pair_index, stat_index]
	midpoints = (np.arange(len(values)) + 0.5) * window_size
	return pd.DataFrame({'chromosome': chromosome, 'position': midpoints, 'value': values})


def load_genomewide(genomewide_path, pair_index, stat_index):
	'''
	Read the genome-wide FST scalar for one population pair and statistic.
	'''
	data = np.load(pathlib.Path(genomewide_path))
	return data[pair_index, stat_index]


def assign_cumulative_positions(windows_by_chromosome):
	'''
	Concatenate per-chromosome windows onto a shared genomic axis, offsetting
	each chromosome by the total span of all preceding chromosomes.
	'''
	offset = 0.0
	frames = []
	for windows in windows_by_chromosome:
		shifted = windows.copy()
		shifted['cumulative_position'] = shifted['position'] + offset
		frames.append(shifted)
		offset += windows['position'].max()
	return pd.concat(frames, ignore_index=True)


def apply_ratio(windows, genomewide_value):
	'''
	Divide each window statistic by the genome-wide value to express window
	differentiation relative to the genome-wide background.
	'''
	windows = windows.copy()
	windows['value'] = windows['value'] / genomewide_value
	return windows


def plot_manhattan(windows, pair_label, output_path):
	'''
	Draw a Manhattan plot of the windowed FST ratio along the genome, with
	alternating colors per chromosome.
	'''
	figure, axis = plt.subplots(figsize=(14, 5))
	chromosomes = windows['chromosome'].unique()
	colors = ['#3b6fb0', '#b0563b']
	ticks = []
	for chromosome_index, chromosome in enumerate(chromosomes):
		subset = windows[windows['chromosome'] == chromosome]
		axis.scatter(subset['cumulative_position'], subset['value'],s=8, color=colors[chromosome_index % 2])
		#axis.scatter(subset['cumulative_position'], np.log2(subset['value']), s=8, color=colors[chromosome_index % 2])             
		ticks.append(subset['cumulative_position'].mean())
	axis.set_xticks(ticks)
	axis.set_xticklabels(chromosomes)
	axis.set_xlabel('Chromosome')
	axis.set_ylabel('FST ratio (window / genome-wide)')
	axis.set_title(pair_label)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)
	figure.tight_layout()
	figure.savefig(output_path, dpi=150)


def build_manhattan(fst_dir, bin_label, chromosomes, window_size, pair_index, stat_index, output_path):
	'''
	Assemble windowed FST across chromosomes for one bin, divide by the bin
	genome-wide FST, and write a Manhattan plot of the ratio to disk.
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
	plot_manhattan(windows, pair_labels[pair_index], output_path)


def main():
	parser = argparse.ArgumentParser(description='Manhattan plot of windowed FST ratio')
	parser.add_argument('--fst-dir', required=True)
	parser.add_argument('--bin', required=True)
	parser.add_argument('--chromosomes', nargs='+', required=True)
	parser.add_argument('--window-size', type=int, default=10000)
	parser.add_argument('--pair-index', type=int, default=0)
	parser.add_argument('--stat-index', type=int, default=2)
	parser.add_argument('--output-path', required=True)
	args = parser.parse_args()
	build_manhattan(args.fst_dir, args.bin, args.chromosomes, args.window_size,
	                args.pair_index, args.stat_index, args.output_path)


if __name__ == '__main__':
	main()
