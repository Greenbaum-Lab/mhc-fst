# FST across time for two populations

FST between two polygon populations of `populations.json`, in focal loci and
genome wide, over every time bin both populations define, with a jackknife
over individuals and a jackknife over blocks of variants.

## Running it

	sbatch fst_time_series.sh

Submit from inside the clone. The job takes the directory it was submitted
from, so the scripts are found wherever the repository sits. Step by step,
from the same directory:

	python run_fst_time_series.py --config config_fst_time.json
	python plot_fst_time_series.py --results results/fst_time_series.csv --output-dir results

After editing the locus list, resolve it before submitting a job. This reads
only the annotation, takes seconds, and names a symbol the annotation does not
carry rather than failing later:

	python -c "import json; from gene_regions import load_gene_spans, locus_gene_names; \
		config = json.load(open('config_fst_time.json')); \
		print(load_gene_spans(config['annotation_path'], locus_gene_names(config['loci'])))"

The estimator can be checked against the one the browser uses, without any
cluster data:

	python check_fst_parity.py --reference ../delphi/analyses/fst.py

## Configuring loci

A locus is either a set of gene symbols, spanning from the first start to the
last end of those genes, or a chromosome with coordinates given directly:

	{"label": "LCT", "genes": ["LCT"]}
	{"label": "TLR6_TLR1_TLR10", "genes": ["TLR6", "TLR1", "TLR10"]}
	{"label": "MHC", "chromosome": "6", "start": 29000000, "end": 35000000}

Every locus is measured twice, over its span and over the span plus 100 kb on
each side, and gets its own figures. Genes named by several loci are resolved
once. Everything else in `config_fst_time.json` is the two polygons, the flank
sizes, the number of variant blocks and the confidence level. Nothing in the
run is random, so the same inputs always give the same numbers.

## Output

`results/fst_time_series.csv` has one row per time bin, target and estimator,
with the estimate, both intervals, both standard errors, the sample counts and
the number of variants and blocks behind it. `results/fst_jackknife.npz` keeps
every leave-one-out value. `results/focal_regions.csv` records the coordinates
each locus resolved to. Figures are named by locus, estimator and jackknife.

## Which variants are used

All of them. There is no call rate threshold, so no variant is dropped for
being thinly covered, and a locus is measured on everything the array carries
there. The only variants that do not contribute are those the Weir & Cockerham
code itself rejects, holding two or fewer called alleles in either population,
which cannot be switched off without changing the estimator.

This means coverage differences between bins are carried into the result
rather than filtered out. A bin whose individuals are thinly covered has
noisier frequencies, and that noise is what the jackknives are there to show.

## Choices that change the numbers

**Estimator.** `ratio_of_averages` sums the Weir & Cockerham components over
the region and divides once. `average_of_ratios` averages the per variant
values, which is what the browser does. They differ most where variants are
few, so a gene span differs more than the flanked version of it.

**Region span.** Span against span plus 100 kb. A short gene carries very few
variants on capture data, and few variants means a noisy estimate.

## The two jackknives

**Over individuals.** Each individual of each population is dropped in turn.
This asks whether a different set of people would give a different answer. On
simulated data the standard error matches the spread of estimates over fresh
draws of individuals: 0.00031 against 0.00040 at n=294/60, 0.00048 against
0.00051 at n=166/29, and 0.00086 against 0.00073 at n=51/11.

**Over blocks of variants.** The variants of a target are cut into 20 blocks
of neighbouring variants holding equal counts, and each block is dropped in
turn. This asks whether a different stretch of genome would give a different
answer, which for a small region is the larger question by far, since such a
region holds few variants and they are linked. Blocks hold equal counts by
construction, so no block weighting is needed. A target with fewer variants
than blocks leaves some blocks empty, and only the blocks holding variants
count.

The two are reported side by side and drawn in separate figures because they
are not interchangeable. For the genome wide line the first is very small and
the second is what matters. For a gene span the second is usually several
times the first.

A bootstrap over individuals was tried first and does not work here.
Resampling with replacement holds the same individual twice, which makes a
replicate more variable than a real sample of its size, and the Weir &
Cockerham correction, computed from the resampled size, does not remove that.
Every replicate then lands above the estimate, by more than the replicates
spread among themselves, so the percentile interval sits entirely above the
estimate and reflecting it puts the interval entirely below. Both miss. A
leave-one-out sample holds nobody twice and has neither problem.

## Assumptions worth knowing

**Pseudo-haploid genotypes.** Most ancient AADR individuals are pseudo-haploid,
one allele called at random and written as a homozygote. The estimator is the
diploid one and counts two alleles per individual, so it treats n individuals
as 2n independent alleles when they are really n. On simulated data at a true
FST of 0.02 this reads back as 0.052 at n=30, 0.036 at n=60 and 0.023 at
n=300, an inflation of about 1/n that fades as the sample grows. Sample size
changes from bin to bin, so part of any trend across time is this artefact
rather than population history. Within a bin both lines carry the same
inflation, so the gap between a focal locus and the background is far more
trustworthy than the height of either. Passing one allele per pseudo-haploid
individual instead of two would remove it, and is not done here.

**Population membership is inherited, not decided here.** Which individual
belongs to which polygon and bin, and which individuals are excluded as low
coverage or outliers, was settled when `populations.json` was built. This
pipeline reads those rosters and does not re-derive them.

**One date per individual.** Each individual sits in exactly one 1000 year bin
by its median date. Dating uncertainty is ignored, and bins are treated as
independent snapshots.

**Coordinates are hg19.** Spans come from the GENCODE release passed in, which
must match the build of the genotypes. GENCODE 19 and the AADR are both hg19.
Coordinates written directly into a locus are trusted as hg19 and not checked.

**Gene symbols are the annotation's, not today's.** GENCODE 19 dates from 2013
and carries the symbols of that time. DARC was renamed ACKR1 in 2015, so the
config asks for DARC. A symbol the annotation does not know stops the run
before any genotype is read.

**Autosomes only,** and the genome wide background includes the focal regions,
which are a negligible fraction of it. A locus on X or Y resolves to a region
holding no variants and produces an empty panel, so G6PD cannot be measured
here. Hemizygous males make the diploid estimator wrong on X, which is why the
chromosome is left out rather than quietly included.
