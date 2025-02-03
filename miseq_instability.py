#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
import argparse as ap

__author__ = 'kc199'

parser = ap.ArgumentParser(description="Calculate instability metrics from Miseq data.")
parser.add_argument("cag", help="Miseq-analyzed CAG distribution.")
parser.add_argument("alleles", help="Identified CAG alleles for this sample.")
parser.add_argument("--peak_bias", type=int, default=0, help="Peak-bias for when to begin summing peak-proportions. "
                                                             "Default = 0 (begin at main peak)")
parser.add_argument("--cutoff_override", action="store_true", help="Override 35 CAG cutoff.")
parser.add_argument("--a1", help="Allele 1 manual override.")
parser.add_argument("--a2", help="Allele 2 manual override.")

args = parser.parse_args()
peak_bias = args.peak_bias

cag = pd.read_csv(args.cag, header=None, names=['cag_dist'])
cag = pd.DataFrame(data=cag.groupby('cag_dist').size(), columns=['cag_count'])

# Make dataframe of all cag sizes
cag_all = pd.DataFrame(data=np.arange(0, 200), columns=['cag_size'])

# Merge left on to all CAG sizes
cag = pd.DataFrame(data=pd.merge(cag_all, cag, how='left', left_index=True, right_index=True))

# Set NaN to zero
cag.loc[pd.isnull(cag['cag_count']), 'cag_count'] = 0

# Get allele information
alleles = pd.read_csv(args.alleles, sep='\t')
a1 = alleles.loc[0, "miseq_cag1"]
a2 = alleles.loc[0, "miseq_cag2"]

total_seq = sum(cag.cag_count)

if a1 == 'None':
    alleles['expansion_index'] = 'NA'
    alleles['peak-proportion_sum'] = 'NA'
    alleles['pps-bias'] = peak_bias
    alleles['total_seq'] = total_seq
    alleles['expanded_seq'] = 'NA'
else:

    normal_cag_count = cag.loc[a2, 'cag_count']
    expanded_cag_count = cag.loc[a1, 'cag_count']

    # Analyze allele2 if it is expanded, cut off >= 35
    if (a1 >= 35 or args.cutoff_override) and expanded_cag_count > 0.01 * normal_cag_count:
        # Select expanded data. Take all data from the expanded allele and forward
        expanded = cag.loc[cag.index >= a1]
        expanded = expanded.sort_values(by='cag_size')
        expanded['normalized'] = expanded.cag_count / sum(expanded.cag_count)
        expanded['locations'] = expanded['cag_size'] - a1

        expansion_index = sum(expanded.loc[expanded.index > a1, ]['normalized'] * expanded.loc[expanded.index > a1, ]['locations'])
        peak_fractional_sum = sum(expanded.loc[expanded.index >= a1 + peak_bias]['cag_count'] / expanded.loc[a1, "cag_count"])

        expanded_seq = sum(expanded.cag_count)
        alleles['expansion_index'] = expansion_index
        alleles['peak-proportion_sum'] = peak_fractional_sum
        alleles['pps-bias'] = peak_bias
        alleles['total_seq'] = total_seq
        alleles['expanded_seq'] = expanded_seq

        #alleles.to_csv(sys.stdout, sep='\t', header=True, index=False)
    else:
        alleles['expansion_index'] = 'NA'
        alleles['peak-proportion_sum'] = 'NA'
        alleles['pps-bias'] = peak_bias
        alleles['total_seq'] = total_seq
        alleles['expanded_seq'] = 'NA'

alleles.to_csv(sys.stdout, sep='\t', header=True, index=False)
