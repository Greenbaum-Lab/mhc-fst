'''
Focal regions resolved from a locus list against a GENCODE GTF.

A locus is either a set of gene symbols, in which case it spans from the first
start to the last end of those genes, or a chromosome and a pair of
coordinates given directly. Positions come from the annotation file, so the
reference build of the regions is the build of the file passed in. GENCODE
release 19 is hg19, the build the AADR genotypes use.
'''

import gzip


def parse_attributes(attributes_field):
	'''
	Parse a GTF attributes field into a dictionary of key to value.
	'''
	attributes = {}
	for part in attributes_field.strip().split(';'):
		part = part.strip()
		if not part or ' ' not in part:
			continue
		key, value = part.split(' ', 1)
		attributes[key] = value.strip().strip('"')
	return attributes


def collect_gene_features(annotation_path, gene_names):
	'''
	Every gene feature of the requested symbols, as chromosome, start and end.
	'''
	wanted = set(gene_names)
	features = {}
	opener = gzip.open if annotation_path.endswith('.gz') else open
	with opener(annotation_path, 'rt') as annotation_file:
		for line in annotation_file:
			if line.startswith('#'):
				continue
			fields = line.rstrip('\n').split('\t')
			if len(fields) < 9 or fields[2] != 'gene':
				continue
			gene_name = parse_attributes(fields[8]).get('gene_name')
			if gene_name in wanted:
				chromosome = fields[0].replace('chr', '')
				features.setdefault(gene_name, []).append((chromosome, int(fields[3]), int(fields[4])))
	return features


def merge_spans(label, spans):
	'''
	The extent of several spans. Spans on more than one chromosome are
	ambiguous and raise rather than being resolved silently.
	'''
	chromosomes = {chromosome for chromosome, _, _ in spans}
	if len(chromosomes) > 1:
		raise ValueError(f'{label} covers chromosomes {sorted(chromosomes)}')
	starts = [start for _, start, _ in spans]
	ends = [end for _, _, end in spans]
	return chromosomes.pop(), min(starts), max(ends)


def load_gene_spans(annotation_path, gene_names):
	'''
	The hg19 span of each requested gene symbol. Raises if a symbol is absent
	from the annotation, so a symbol the annotation does not carry stops the
	run instead of becoming a region without variants.
	'''
	features = collect_gene_features(annotation_path, gene_names)
	missing = set(gene_names) - set(features)
	if missing:
		raise KeyError(f'Genes not found in annotation: {sorted(missing)}')
	return {gene_name: merge_spans(gene_name, features[gene_name]) for gene_name in gene_names}


def locus_gene_names(loci):
	'''
	Every gene symbol the locus list names, which is what has to be resolved.
	'''
	return [gene_name for locus in loci for gene_name in locus.get('genes', [])]


def locus_span(locus, gene_spans):
	'''
	Chromosome, start and end of a locus, either the extent of the genes it
	names or the coordinates it gives directly.
	'''
	if 'chromosome' in locus:
		return locus['chromosome'], locus['start'], locus['end']
	return merge_spans(locus['label'], [gene_spans[gene_name] for gene_name in locus['genes']])


def region_label(locus_label, flank_size):
	if flank_size == 0:
		return f'{locus_label}_span'
	return f'{locus_label}_flank{flank_size // 1000}kb'


def locus_regions(locus, gene_spans, flank_sizes):
	chromosome, start, end = locus_span(locus, gene_spans)
	return [
		{
			'region_id': region_label(locus['label'], flank_size),
			'locus': locus['label'],
			'chromosome': chromosome,
			'start': max(1, start - flank_size),
			'end': end + flank_size,
		}
		for flank_size in flank_sizes
	]


def build_regions(loci, gene_spans, flank_sizes):
	'''
	One region per locus and flank size, the locus span padded on both sides.
	'''
	return [region for locus in loci for region in locus_regions(locus, gene_spans, flank_sizes)]
