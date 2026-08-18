'''
Check that the FST rebuilt from the components equals delphi's implementation.

The reference function is compiled straight out of the delphi source instead of
being copied here, so the check fails if the two implementations ever diverge.
'''

import ast
import argparse
import numpy as np

from fst_time.fst_core import weir_cockerham_components, per_variant_fst

DEFAULT_REFERENCE_PATH = '../delphi/analyses/fst.py'
REFERENCE_FUNCTION = '_compute_fst'


def load_reference_function(module_path, function_name):
	'''
	Compile one function out of a module without importing the module, whose
	own imports are not available outside the delphi runtime.
	'''
	with open(module_path) as module_file:
		tree = ast.parse(module_file.read())
	definition = next(
		node for node in tree.body
		if isinstance(node, ast.FunctionDef) and node.name == function_name)
	namespace = {'np': np}
	module = ast.Module(body=[definition], type_ignores=[])
	exec(compile(module, module_path, 'exec'), namespace)
	return namespace[function_name]


def random_inputs(variant_count, random_generator):
	'''
	Allele statistics covering the awkward cases too: uncalled variants,
	single called individuals and monomorphic variants.
	'''
	an1 = random_generator.integers(0, 60, variant_count).astype(float)
	an2 = random_generator.integers(0, 60, variant_count).astype(float)
	ac1 = np.floor(random_generator.random(variant_count) * (an1 + 1))
	ac2 = np.floor(random_generator.random(variant_count) * (an2 + 1))
	het1 = np.minimum(ac1, an1 / 2.0)
	het2 = np.minimum(ac2, an2 / 2.0)
	return an1, an2, ac1, ac2, het1, het2


def main():
	parser = argparse.ArgumentParser(description='Compare the component form of FST with the delphi implementation')
	parser.add_argument('--reference', default=DEFAULT_REFERENCE_PATH)
	parser.add_argument('--variants', type=int, default=200000)
	parser.add_argument('--seed', type=int, default=0)
	args = parser.parse_args()
	reference_function = load_reference_function(args.reference, REFERENCE_FUNCTION)
	inputs = random_inputs(args.variants, np.random.default_rng(args.seed))
	expected = reference_function(*inputs)
	observed = per_variant_fst(*weir_cockerham_components(*inputs))
	np.testing.assert_array_equal(np.isnan(expected), np.isnan(observed))
	np.testing.assert_array_equal(expected[~np.isnan(expected)], observed[~np.isnan(observed)])
	print(f'identical on {args.variants} variants, {np.count_nonzero(~np.isnan(expected))} usable')


if __name__ == '__main__':
	main()
