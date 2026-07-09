import argparse
import pathlib
import numpy as np
import pandas as pd

SNP_COUNT_STAT_INDEX = 6

def load_bin_windows(fst_dir, bin_label, chromosomes, window_size, pair_index, stat_index):
	'''
	Load every window across the given chromosomes for one bin, keeping the
	window position, the ranking statistic, and its SNP count.
	'''
	frames = []
	for chromosome in chromosomes:
		path = pathlib.Path(f'{fst_dir}/{bin_label}/fst_{bin_label}_fst_{chromosome}.npy')
		data = np.load(path)
		values = data[:, pair_index, stat_index]
		counts = data[:, pair_index, SNP_COUNT_STAT_INDEX]
		midpoints = (np.arange(len(values)) + 0.5) * window_size
		frames.append(pd.DataFrame({
			'chromosome': chromosome,
			'position': midpoints,
			'fst': values,
			'snp_count': counts,
		}))
	return pd.concat(frames, ignore_index=True)


def filter_by_min_snps(windows, min_snps):
	'''
	Drop windows below the minimum SNP count. Must run before the null is
	built, since sparse windows produce the extreme-FST outliers that would
	otherwise contaminate the null's tail.
	'''
	return windows[windows['snp_count'] >= min_snps].reset_index(drop=True)


def empirical_one_sided_p(values):
	'''
	Return (p_high, p_low) empirical one-sided p-values for each value
	against the full within-bin distribution of `values` itself (no
	simulation):
		p_high(w) = (#windows with FST >= FST(w)) / N   -- peak candidates
		p_low(w)  = (#windows with FST <= FST(w)) / N   -- valley candidates
	Ties are counted on both sides, so a tied value gets the same p on
	whichever side it's evaluated from.
	'''
	n = len(values)
	sorted_values = np.sort(values)
	count_at_least = n - np.searchsorted(sorted_values, values, side='left')
	count_at_most = np.searchsorted(sorted_values, values, side='right')
	return count_at_least / n, count_at_most / n


def benjamini_hochberg(p_values):
	'''
	Return BH-adjusted q-values for an array of p-values.
	'''
	n = len(p_values)
	order = np.argsort(p_values)
	ranked = p_values[order] * n / (np.arange(n) + 1)
	q_ranked = np.minimum.accumulate(ranked[::-1])[::-1]
	q_ranked = np.clip(q_ranked, 0, 1)
	q_values = np.empty(n)
	q_values[order] = q_ranked
	return q_values


def label_candidates(q_high, q_low, q_threshold):
	'''
	Label each window 'peak' (high-FST outlier), 'valley' (low-FST outlier),
	or 'background', preferring peak on the rare tie so a window is never
	both.
	'''
	is_peak = q_high <= q_threshold
	is_valley = (q_low <= q_threshold) & ~is_peak
	labels = np.full(len(q_high), 'background', dtype=object)
	labels[is_valley] = 'valley'
	labels[is_peak] = 'peak'
	return labels


def scan_bin(fst_dir, bin_label, chromosomes, window_size, pair_index, stat_index, min_snps, q_threshold):
	'''
	Run the full within-bin outlier scan: filter sparse windows, rank the
	remainder by raw FST, compute empirical one-sided p-values against that
	bin's own filtered window distribution, apply BH-FDR separately for the
	high and low scans, and label candidates at q_threshold.
	'''
	windows = load_bin_windows(fst_dir, bin_label, chromosomes, window_size, pair_index, stat_index)
	windows = filter_by_min_snps(windows, min_snps)
	p_high, p_low = empirical_one_sided_p(windows['fst'].to_numpy())
	windows['p_high'] = p_high
	windows['p_low'] = p_low
	windows['q_high'] = benjamini_hochberg(p_high)
	windows['q_low'] = benjamini_hochberg(p_low)
	windows['candidate'] = label_candidates(windows['q_high'].to_numpy(), windows['q_low'].to_numpy(), q_threshold)
	return windows


def main():
	parser = argparse.ArgumentParser(description='Within-bin empirical outlier scan of windowed FST')
	parser.add_argument('--fst-dir', required=True)
	parser.add_argument('--bin', required=True)
	parser.add_argument('--chromosomes', nargs='+', type=int, required=True)
	parser.add_argument('--window-size', type=int, required=True)
	parser.add_argument('--pair-index', type=int, default=0)
	parser.add_argument('--stat-index', type=int, default=0)
	parser.add_argument('--min-snps', type=int, required=True)
	parser.add_argument('--q-threshold', type=float, default=0.05)
	parser.add_argument('--output-path', required=True)
	args = parser.parse_args()
	windows = scan_bin(args.fst_dir, args.bin, args.chromosomes, args.window_size,
	                   args.pair_index, args.stat_index, args.min_snps, args.q_threshold)
	windows.to_csv(args.output_path, index=False)


if __name__ == '__main__':
	main()
