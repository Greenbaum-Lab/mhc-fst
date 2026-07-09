import numpy as np
import pandas as pd


def load_janno(janno_path):
	'''
	Load a Poseidon .janno table, keeping only the fields needed for
	spatiotemporal population assignment, and derive years before present
	from the median calendar date.
	Samples without a numeric median date are dropped.
	'''
	columns = ['Poseidon_ID', 'Latitude', 'Longitude', 'Date_BC_AD_Median']
	table = pd.read_csv(janno_path, sep='\t', usecols=columns, na_values=[])
	table['Latitude'] = pd.to_numeric(table['Latitude'], errors='coerce')
	table['Longitude'] = pd.to_numeric(table['Longitude'], errors='coerce')
	table['Date_BC_AD_Median'] = pd.to_numeric(table['Date_BC_AD_Median'], errors='coerce')
	table = table.dropna(subset=['Latitude', 'Longitude', 'Date_BC_AD_Median'])
	table['date'] = np.maximum(0, 1950 - table['Date_BC_AD_Median']).astype(int)
	return table


def assign_region(latitude, longitude, region_boxes):
	'''
	Return the label of the first region box containing the coordinate,
	or None if the coordinate falls outside every box.
	'''
	for label, box in region_boxes.items():
		low_lat, high_lat = box['lat']
		low_lon, high_lon = box['lon']
		if low_lat <= latitude <= high_lat and low_lon <= longitude <= high_lon:
			return label
	return None


def assign_bin(date, time_bins):
	'''
	Return the '{low}_{high}' label of the first time bin containing the
	date (low inclusive, high inclusive), or None if outside every bin.
	'''
	for low, high in time_bins:
		if low <= date <= high:
			return f'{low}_{high}'
	return None


def assign_populations(janno_path, region_boxes, time_bins):
	'''
	Assign each sample to a population labelled '{region}_{low}_{high}'
	based on its coordinate and years-before-present date.
	Samples outside all region boxes or all time bins are dropped.
	Returns a DataFrame with columns Poseidon_ID and population.
	'''
	table = load_janno(janno_path)
	table['region'] = table.apply(
		lambda row: assign_region(row['Latitude'], row['Longitude'], region_boxes), axis=1)
	table['bin'] = table['date'].apply(lambda date: assign_bin(date, time_bins))
	table = table.dropna(subset=['region', 'bin'])
	table['population'] = table['region'] + '_' + table['bin']
	return table[['Poseidon_ID', 'population']].reset_index(drop=True)
