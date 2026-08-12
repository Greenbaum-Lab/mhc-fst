# FST across time for two populations

FST between two polygon populations of `populations.json`, in focal gene
regions and genome wide, over every time bin both populations define, with a
bootstrap over individuals.

## Running it

	python run_fst_time_series.py --config config_fst_time.json
	python plot_fst_time_series.py --results results/fst_time_series.csv --output-dir results
	python check_fst_parity.py --reference ../delphi/analyses/fst.py

Everything is set in `config_fst_time.json`: the two polygons, the gene list,
the flank sizes, the call rate threshold, the number of replicates and the
seed. Changing the genes or swapping a population is a config edit, not a code
edit.

## Output

`results/fst_time_series.csv` has one row per time bin, target, SNP set and
estimator, with the point estimate, both bootstrap intervals, the sample
counts and the number of variants the estimate used.
`results/fst_bootstrap.npz` keeps every replicate value, so comparisons
between a region and the background can be made later without the genotypes.
`results/focal_regions.csv` records the coordinates the gene symbols resolved
to.

## Choices that change the numbers

Four choices are reported side by side rather than fixed, because each one
moves the result and none of them is obviously right.

**SNP set.** Ancient samples are patchy, so a variant is only usable where
enough individuals have a call. Call rate is the fraction of individuals in a
population and bin with a call, and the threshold is applied to the lower of
the two populations.

- `none` keeps every variant.
- `per_bin` keeps variants above the threshold in that bin, so the variant set
  changes from bin to bin with coverage. Part of any trend can then come from
  coverage rather than from allele frequencies.
- `intersection` keeps variants above the threshold in every bin, so all bins
  share one fixed set. Fewer variants, but the trend across time is comparable.

**Estimator.** `ratio_of_averages` sums the Weir & Cockerham components over
the region and divides once. `average_of_ratios` averages the per variant
values, which is what the browser does. They differ most where variants are
few, so the gene body regions differ more than the flanked ones.

**Region span.** Each gene is measured twice, over the gene body and over the
body plus 100 kb on each side. The gene body of SLC24A5 carries very few
variants on capture data, and few variants means a noisy estimate.

**Bootstrap interval.** Resampling individuals with replacement duplicates
individuals, which makes a replicate more variable than a real sample of the
same size. The Weir & Cockerham sample size correction is computed from the
resampled size and does not remove that, so replicates sit above the point
estimate by roughly the reciprocal of the smaller population size, about
+0.006 at n=60 and +0.03 at n=11. `ci_low_basic` and `ci_high_basic` reflect
the interval about the point estimate and remove the shift, and are what the
figures show. `ci_low_percentile` and `ci_high_percentile` are the raw
interval, kept for comparison. A basic interval reaching below zero means the
value cannot be told apart from no differentiation.

## Assumptions worth knowing

**Pseudo-haploid genotypes.** Most ancient AADR individuals are pseudo-haploid,
one allele called at random and written as a homozygote. The estimator is the
diploid one and counts two alleles per individual, so it treats n individuals
as 2n independent alleles when they are really n. On simulated data at a true
FST of 0.02 this reads back as 0.052 at n=30, 0.036 at n=60 and 0.023 at
n=300, an inflation of about 1/n that fades as the sample grows. Sample size
changes from bin to bin, so part of any trend across time is this artefact
rather than population history. Within a bin both lines carry the same
inflation, so the gap between a focal region and the background is far more
trustworthy than the height of either. Passing one allele per pseudo-haploid
individual instead of two would remove it, and is not done here.

**Population membership is inherited, not decided here.** Which individual
belongs to which polygon and bin, and which individuals are excluded as low
coverage or outliers, was settled when `populations.json` was built. This
pipeline reads those rosters and does not re-derive them.

**One date per individual.** Each individual sits in exactly one 1000 year bin
by its median date. Dating uncertainty is ignored, and bins are treated as
independent snapshots.

**Coordinates are hg19.** Gene spans come from the GENCODE release passed in,
which must match the build of the genotypes. GENCODE 19 and the AADR are both
hg19.

**Autosomes only,** and the genome wide background includes the focal regions,
which are a negligible fraction of it.

**A rule inside the estimator cannot be switched off.** Variants with two or
fewer called alleles in either population are dropped by the Weir & Cockerham
code itself. The `none` SNP set therefore means no call rate threshold, not
literally every variant.
