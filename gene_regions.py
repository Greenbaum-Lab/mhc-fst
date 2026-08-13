'''
Focal regions resolved from gene symbols against a GENCODE GTF.

Positions come from the annotation file, so the reference build of the regions
is the build of the file that is passed in. GENCODE release 19 is hg19, the
build the AADR genotypes use.
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


def resolve_span(gene_name, gene_features):
	'''
	The span of one gene, merging several annotated features of the same
	symbol. A symbol annotated on more than one chromosome is ambiguous and
	raises rather than being resolved silently.
	'''
	chromosomes = {chromosome for chromosome, _, _ in gene_features}
	if len(chromosomes) > 1:
		raise ValueError(f'Gene {gene_name} is annotated on chromosomes {sorted(chromosomes)}')
	starts = [start for _, start, _ in gene_features]
	ends = [end for _, _, end in gene_features]
	return chromosomes.pop(), min(starts), max(ends)


def load_gene_spans(annotation_path, gene_names):
	'''
	The hg19 span of each requested gene symbol. Raises if a symbol is absent
	from the annotation, so a typo cannot pass as a region without variants.
	'''
	features = collect_gene_features(annotation_path, gene_names)
	missing = set(gene_names) - set(features)
	if missing:
		raise KeyError(f'Genes not found in annotation: {sorted(missing)}')
	return {gene_name: resolve_span(gene_name, features[gene_name]) for gene_name in gene_names}


def region_label(gene_name, flank_size):
	if flank_size == 0:
		return f'{gene_name}_body'
	return f'{gene_name}_flank{flank_size // 1000}kb'


def build_regions(gene_spans, flank_sizes):
	'''
	One region per gene and flank size, the gene span padded on both sides.
	'''
	return [
		{
			'region_id': region_label(gene_name, flank_size),
			'gene': gene_name,
			'chromosome': chromosome,
			'start': max(1, start - flank_size),
			'end': end + flank_size,
		}
		for gene_name, (chromosome, start, end) in gene_spans.items()
		for flank_size in flank_sizes
	]
