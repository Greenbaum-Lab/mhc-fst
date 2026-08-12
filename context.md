# Figure 0 / Case Study Pipeline — Context

This file defines the goal, scope, and open decisions for the analysis that will
produce the general (non-MHC) schematic figure of the paper (currently Fig 0,
placeholder). It is meant to keep any pipeline work anchored to the paper's
argument and to make explicit what is decided vs. still open.

## 1. Purpose in the paper

The paper's core mechanism (currently illustrated only as a hand-drawn schematic
for the MHC, Fig 1) is general and not specific to immunity:

- When gene flow occurs between two populations, genome-wide FST is expected to
  decline gradually, dominated by drift/admixture.
- A locus under selection can decouple from that neutral trajectory in two
  directions, depending on whether the selective pressure travels with the
  migrants or stays local:
  - **Transmissible pressure** (e.g. a shared pathogen, a cultural practice such
    as dairying) -> selection reinforces convergence -> FST at the locus drops
    earlier / faster than the genome-wide background.
  - **Non-transmissible / place-bound pressure** (e.g. UV exposure, altitude, a
    locally restricted pathogen) -> local adaptation persists despite gene flow
    -> FST at the locus stays elevated relative to the genome-wide background.

Fig 0's job is to demonstrate this general principle empirically, using loci
that are NOT the MHC and are independently established in the literature, before
the paper narrows to the MHC as its focal case study in Sections 2-4. This gives
the reader the mechanism first, then the MHC-specific application.

## 2. Confirmed design decisions

- **Ancient DNA source:** AADR (Allen Ancient DNA Resource), already available
  on the working cluster. Exact AADR version/release to be confirmed before
  any data pull.
- **Case study framings:** BOTH of the following, one case study per framing
  (not necessarily the same locus for both):
  - **Time-axis case study:** FST between a fixed population pair (or a
    population and its ancestral/pre-admixture source) tracked across ancient
    time bins/transects at a candidate locus, compared against the
    genome-wide FST trajectory over the same time bins. This mirrors the
    logic of paper Section 2 (Europe vs. East Asia).
  - **Chromosome-axis case study:** FST at a candidate locus at a fixed
    population pair / timepoint, shown against the genome-wide 1 Mb window
    background (same style as the empirical null used for MHC in Section 3 /
    Fig. S2), highlighting the locus's departure from the background
    distribution.
- **Locus selection:** User-provided (not a fresh genome-wide scan for
  candidates). Genes/loci will be supplied by name; positions must be
  resolved against the **hg19/GRCh37** reference (same build already used
  elsewhere in the project), converting gene symbol -> genomic coordinates
  before any FST computation.

## 3. Open questions (need to be confirmed before implementation)

These are explicitly NOT decided yet and should not be assumed:

1. **Population pairs for the time-axis case studies** — which ancient
   population(s)/time transect(s) define the "before" and "after" gene-flow
   states for each candidate locus? Does this reuse a specific published
   admixture event (e.g. Steppe expansion, Neolithic transition) or is it
   locus-specific?
2. **Time binning** — fixed archaeological period bins (e.g. Mesolithic /
   Neolithic / Bronze Age) vs. continuous/sliding time windows vs.
   admixture-date-anchored bins per population.
3. **Genome-wide background reuse** — should the chromosome-axis case study
   reuse the exact 1 Mb window null-distribution pipeline built for the MHC
   modern-data analysis (Section 3), applied here to ancient genotypes, or
   does ancient data (lower coverage, imputation, ascertainment) require a
   separate background procedure?
4. **FST estimator/correction** — should case studies use the same corrected
   Weir & Cockerham FST implementation (the in-house Python package used for
   the MHC analysis), applied without modification to ancient data, or does
   the correction (polymorphism-dependent bound) need reassessment given
   typically lower SNP density in ancient/capture data?
5. **Locus list** — to be provided by the user; positions to be resolved from
   gene symbol to hg19 coordinates once supplied (need a defined window size
   around each gene, e.g. gene body only vs. +/- flanking region).
6. **Sample size / coverage filtering** for ancient individuals per time bin
   per population (minimum count, coverage/genotype quality thresholds) —
   not yet defined.

## 4. Explicitly out of scope for Fig 0

- No MHC loci in this figure (MHC is reserved for Sections 1-4).
- No de novo genome-wide selection scan to discover new candidate loci.
- No modern-only data — this figure is anchored on ancient DNA (AADR).

## 5. Next steps

- User supplies the candidate gene/locus list.
- Confirm AADR version and cluster path.
- Resolve open questions in Section 3 (population pairs, time bins, background
  method, FST estimator reuse) before writing any pipeline code.
