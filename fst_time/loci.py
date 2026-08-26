'''
The loci a run measures and the figure draws.

This list is the one place the loci are chosen. A locus carries a label, the
gene symbols it spans, what it is studied for, and the trend expected of it
against the genome wide background, one of high, low or neutral. A locus of
several symbols spans from the first start to the last end of them.

`time_bp` is the year the expectation changes, or None where there is none.
Neither the trend nor that year enters the computation, only the figure, which
is grouped into a column per trend and marks the year.

Symbols are resolved against the GENCODE release the config names. Release 19
carries the gene names of 2013, so a symbol renamed since then is asked for by
its name of that time.
'''

FOCAL_LOCI = [
	{'label': 'SLC24A5', 'genes': ['SLC24A5'], 'phenotype': 'pigmentation', 'trend': 'high', 'time_bp': 8500},
	{'label': 'SLC25A2', 'genes': ['SLC25A2'], 'phenotype': 'pigmentation', 'trend': 'high', 'time_bp': 4500},
	{'label': 'OCA2', 'genes': ['OCA2'], 'phenotype': 'pigmentation', 'trend': 'high', 'time_bp': 8500},
	{'label': 'EDAR', 'genes': ['EDAR'], 'phenotype': 'ectodermal morphology', 'trend': 'high', 'time_bp': None},
	{'label': 'EPAS1', 'genes': ['EPAS1'], 'phenotype': 'high altitude hypoxia', 'trend': 'high', 'time_bp': None},
	{'label': 'LCT', 'genes': ['LCT'], 'phenotype': 'lactase persistence', 'trend': 'low', 'time_bp': 4500},
	{'label': 'FADS1_FADS2', 'genes': ['FADS1', 'FADS2'], 'phenotype': 'lipid metabolism', 'trend': 'low', 'time_bp': 8500},
	{'label': 'ENSA', 'genes': ['ENSA'], 'phenotype': 'insulin signalling', 'trend': 'low', 'time_bp': 8500},
	{'label': 'TLR6_TLR1_TLR10', 'genes': ['TLR6', 'TLR1', 'TLR10'], 'phenotype': 'innate immunity', 'trend': 'low', 'time_bp': 8500},
	{'label': 'MARK3', 'genes': ['MARK3'], 'phenotype': 'infectious disease immunity', 'trend': 'low', 'time_bp': 4500},
	{'label': 'IL23R', 'genes': ['IL23R'], 'phenotype': 'immune signalling', 'trend': 'low', 'time_bp': 4500},
	{'label': 'TYK2', 'genes': ['TYK2'], 'phenotype': 'immune signalling', 'trend': 'low', 'time_bp': 4500},
	{'label': 'ILRL1', 'genes': ['IL1RL1'], 'phenotype': 'immune signalling', 'trend': 'low', 'time_bp': 3150},
	{'label': 'ERAP', 'genes': ['ERAP2'], 'phenotype': 'antigen processing, plague', 'trend': 'low', 'time_bp': 700},
	{'label': 'ABO', 'genes': ['ABO'], 'phenotype': 'blood group, pathogen binding', 'trend': 'low', 'time_bp': None},
]
