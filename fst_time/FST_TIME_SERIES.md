# FST across time for two populations

FST between two polygon populations of `populations.json`, in focal loci and
genome wide, over every time bin both populations define, with a jackknife
over individuals and a jackknife over blocks of variants. The estimate is the
Weir & Cockerham ratio of averages over the variants a locus spans.

## Layout

The repository holds two projects, each in its own folder with the code, the
configuration and the job scripts it needs. `fst_time` is this one. The older
window scan lives in `window_scan` and shares only two helpers with it, the
genomic axis and the Benjamini Hochberg procedure, which `manhattan_genes.py`
imports from there rather than repeating.

Both are Python packages, so their entry points run as modules from the
repository root and never care where the clone sits.

## Running it

	sbatch fst_time/fst_time_series.sh

Submit from the repository root. The job takes the directory it was submitted
from. Step by step, from the same place:

	python -m fst_time.run_fst_time_series --config fst_time/config_fst_time.json
	python -m fst_time.plot_trend_grid --results results/fst_time_series.csv \
		--periods fst_time/time_periods.json \
		--gene-background results/gene_background.csv --output-dir results

After editing the locus list, resolve it before submitting a job. This reads
only the annotation, takes seconds, and names a symbol the annotation does not
carry rather than failing later:

	python -c "import json; from fst_time.loci import FOCAL_LOCI; \
		from fst_time.gene_regions import load_gene_spans, locus_gene_names; \
		config = json.load(open('fst_time/config_fst_time.json')); \
		print(load_gene_spans(config['annotation_path'], locus_gene_names(FOCAL_LOCI)))"

`plot_trend_grid.py` is the one figure a run draws: every locus on one page, in
a column for the trend expected of it, with the archaeological periods shaded
and the expected time of change marked. Each panel holds that locus, the genome
wide line and the mean over all annotated genes, the last two identical in every
panel. The loci drawn are the loci measured, so the figure is changed by editing
`loci.py` and running again. A trend no locus expects takes no column.

It reads only the results table, `gene_background.csv` and
`fst_time/time_periods.json`, so copy them into one directory and run it there:

	python plot_trend_grid.py

`plot_final_figure.py` is the finished version of that page for the paper. It
holds its own locus list and its own periods, so it needs only the results
table and `gene_background.csv`, and it draws both jackknives in one run:

	python plot_final_figure.py --results-dir . --output-dir .

The loci it draws, their order and their columns are the list at the top of
the script, which is a choice of the figure rather than of the run, so a
measured locus can be left off the page without measuring anything again.
Every panel of a column shares one FST axis, the time axis ends where the data
ends, and the mean over genes is drawn as a line with a band like the others
rather than as error bars.

`manhattan_genes.py` scans every gene of a run, one figure per time bin,
written into a folder of its own:

	python -m fst_time.manhattan_genes --per-gene results/fst_per_gene.csv \
		--output-dir results/manhattan

The estimator can be checked against the one the browser uses, without any
cluster data:

	python -m fst_time.check_fst_parity --reference ../delphi/analyses/fst.py

## Configuring loci

The loci are chosen in `loci.py` and nowhere else. The configuration file holds
only paths and parameters, so changing which loci are studied is an edit to one
python list, and both the computation and the figure follow it.

A locus is either a set of gene symbols, spanning from the first start to the
last end of those genes, or a chromosome with coordinates given directly. Each
one carries what it is studied for and what is expected of it, which is what
the figure is titled and grouped by:

	{'label': 'LCT', 'genes': ['LCT'],
	 'phenotype': 'lactase persistence', 'trend': 'low', 'time_bp': 4500}
	{'label': 'FADS1_FADS2', 'genes': ['FADS1', 'FADS2'],
	 'phenotype': 'lipid metabolism', 'trend': 'low', 'time_bp': 8500}

`trend` is what the locus is expected to do against the genome wide
background, one of high, low or neutral. `time_bp` is the year the expectation
changes, or None where there is none. Neither enters the computation, only the
figure. A locus is measured over its span alone, on its own variants and with
its own jackknives, not assembled from the per gene table. Genes named by
several loci are resolved once. Nothing in the run is random, so the same
inputs always give the same numbers.

## Output

`results/fst_time_series.csv` has one row per time bin and target, with the
estimate, both intervals, both standard errors, the sample counts and the
number of variants and blocks behind it. `results/fst_jackknife.npz` keeps
every leave-one-out value. `results/focal_regions.csv` records the coordinates
each locus resolved to. Figures are named by locus and jackknife.

`results/fst_per_gene.csv` holds every gene of the annotation, one row per
gene and time bin, with its own FST and the variants it rests on.
`results/fst_per_gene.npz` holds the sums those came from, so any set of genes
can be recombined without the genotypes. `results/gene_background.csv` is the
mean over genes at each time bin, which the figures draw.

## The background of genes

Every gene of the annotation is measured, not only the focal loci, and the
mean over those genes is drawn beside each locus. The genome wide line is a
different kind of thing: it pools a million variants into one number, where a
gene holds a few dozen, so a locus sitting above or below it may only be
showing that genes are not the genome. Averaging over genes compares like with
like.

The error bars on that mean are the spread between genes, not the uncertainty
of any one of them. Every gene the annotation carries counts towards it: there
is no threshold on the number of variants a gene holds, since that would be a
filter on the data. A gene the estimator left without a usable variant has no
value and so does not enter the mean.

A gene here is a `gene` feature of the annotation whose `gene_type` is one of
`gene_biotypes`, which is `protein_coding` alone, so pseudogenes and long
non-coding genes are left out of the background and out of the saved table.
Adding a biotype is a config edit. The span runs from the first base of the
gene to its last, introns and untranslated regions included, and a variant
inside two overlapping genes counts for both.

## Which variants are used

All of them. There is no call rate threshold, so no variant is dropped for
being thinly covered, and a locus is measured on everything the array carries
there. The only variants that do not contribute are those the Weir & Cockerham
code itself rejects, holding two or fewer called alleles in either population,
which cannot be switched off without changing the estimator.

This means coverage differences between bins are carried into the result
rather than filtered out. A bin whose individuals are thinly covered has
noisier frequencies, and that noise is what the jackknives are there to show.

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

## Scanning every gene

`manhattan_genes.py` asks which genes stand out, and writes
`manhattan/gene_significance.csv` with a p and q value for every gene,
`manhattan/top_genes.csv` with the named ones, and one figure per time bin.

A gene measured from twenty variants scatters far more than one measured from
six hundred, so genes are ranked within strata of variant count, ten by
default. Within a stratum the null is the bulk of the genes themselves, its
centre and spread taken from the median and the median absolute deviation so
that the outliers being looked for do not widen the null they are measured
against. A gene's p value is the normal tail beyond its distance from that
centre, and Benjamini Hochberg turns those into q values.

Ranking a gene against the empirical distribution of the others cannot serve
here, because a p value read off the ranks is uniform by construction and
every q value comes back at one. That fraction is still written out, as
`percentile`, since it says plainly where a gene sits.

Two assumptions come with the test. Most genes are taken to be ordinary, which
is what makes the bulk a null, and the tail is taken to be normal, which the
skew of an FST distribution only approximately obeys. Read a p value as a way
of ordering candidates rather than as a precise probability.

The axis is capped, at 30 by default, and genes beyond it are drawn as
triangles with their number given in the title, so one gene far out in the
tail does not flatten the rest.

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
and carries the symbols of that time, so a gene renamed since then is asked for
by its name of that time. A symbol the annotation does not know stops the run
before any genotype is read.

**Autosomes only,** and the genome wide background includes the focal regions,
which are a negligible fraction of it. A locus on X or Y holds no variants and
its panel is drawn empty and labelled as such. Hemizygous males make the
diploid estimator wrong on X, which is why the chromosome is left out rather
than quietly included.
