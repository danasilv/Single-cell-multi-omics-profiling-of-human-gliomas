#!/usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
from collections import Counter

import numpy as np
import pandas as pd
import pickle
import sys
import os
import glob

DIRECTORY = sys.argv[1]

# GLOBAL VARIABLES: specify for individual datasets
MET_FILE = os.path.join(DIRECTORY, 'site_matrix.csv')
BED_INTERSECT = os.path.join(DIRECTORY, 'site_intersect.bed')

# Function to create site -> gene dictionary from intersect file
CHROM_ID_LENGTH = 10
CHROM_COL = 5
POS_COL = 6
GENE_COL = 3

def create_gene_dict(frame):
    chroms = list(frame[CHROM_COL].values)
    pos = list(frame[POS_COL].values)
    sites = [i + '-' + str(j).zfill(CHROM_ID_LENGTH) for i,j in zip(chroms, pos)]
    genes = list(frame[GENE_COL].values)
    return dict(zip(sites, genes))

# Load intersectBed output and create dictionary
intersect_df = pd.read_table(BED_INTERSECT, delim_whitespace=True, header=None)
gene_dict = create_gene_dict(intersect_df)

met = pd.read_csv(MET_FILE)
print(met.head())

# Remove sites that do not overlap with any of the gene windows and add gene column (takes a while)
PRINT_EVERY = 100000

def remove_sites(frame, site_dict):
    new_frame = frame
    new_frame.loc[:, 'gene_name'] = None
    for idx, row in new_frame.iterrows():
        if idx % PRINT_EVERY == 0:
            print('Processing Site: ' + str(idx))
        try:
            row['gene_name'] = site_dict[row['site']]
        except KeyError:
            continue
    return new_frame[new_frame['gene_name'].notnull()]

met_subset_cut = remove_sites(met, gene_dict)

# Compute fraction of methylated sites for each gene
CPG_THRESHOLD = 5
def met_pct(x):
    total_calls = x.count('1') + x.count('0')
    if total_calls >= CPG_THRESHOLD:
        return x.count('1') / float(total_calls)
    else:
        return None

def met_count(x):
    total_calls = x.count('1') + x.count('0')
    if total_calls >= CPG_THRESHOLD:
        return total_calls
    else:
        return None

met_by_gene_pct = (met_subset_cut.drop('site', axis=1).groupby(['gene_name']).sum()).applymap(met_pct).reset_index()
met_by_gene_raw  = (met_subset_cut.drop('site', axis=1).groupby(['gene_name']).sum()).applymap(met_count).reset_index()

# Write DataFrames to CSV
met_by_gene_pct.to_csv(os.path.join(DIRECTORY, 'region_matrix_pct.csv'), index = False)
met_by_gene_raw.to_csv(os.path.join(DIRECTORY, 'region_matrix_count.csv'), index = False)
