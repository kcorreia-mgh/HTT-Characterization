#!/usr/bin/env python

import sys
import argparse as ap
import numpy as np
import pandas as pd

__author__ = 'kc199'


def results_check(res, res_start):
    res.seek(res_start)

    f = res.readlines()

    if len(f) <= 5:
        return False
    else:
        return True


def build_cag_dist(cag_dist):
    # Read cag sizes and build a distribution
    cag = pd.read_csv(cag_dist, header=None, names=['cag_dist'])
    cag = pd.DataFrame(data=cag.groupby('cag_dist').size(), columns=['cag_count'])
    # Make dataframe of all cag sizes
    cag_all = pd.DataFrame(data=np.arange(0, 80), columns=['cag_size'])

    # Merge left on to all CAG sizes
    cag = pd.DataFrame(data=pd.merge(cag_all, cag, how='left', left_index=True, right_index=True))
    cag.loc[pd.isnull(cag['cag_count']), 'cag_count'] = 0

    return cag


def find_peak_regions(cag):
    # Smooth the distribution
    cag['smooth'] = cag.cag_count.rolling(window=5).mean()
    # Replace NaN with 0
    cag.loc[pd.isnull(cag.smooth), 'smooth'] = 0

    # Define peak regions as areas > threshold
    threshold = np.mean(cag.loc[cag.smooth != 0, 'smooth'])
    cag['peak_regions'] = False
    cag.loc[cag.smooth >= threshold, 'peak_regions'] = True

    return cag


def find_peak_height(cag, peaks):
    heights = []
    for p in peaks:
        height = cag[cag.cag_size == p]['cag_count'].values[0]
        heights.append(height)

    return max(heights)


def get_peak_regions(cag):
    peaks = cag[cag.peak_regions]
    peak_regions = []

    prior = peaks.index[0]
    start = 0
    for i, v in enumerate(peaks.index[1:]):
        if i == len(peaks.index[1:]) - 1:
            peak_regions.append(peaks.index.values[start:])
            break
        elif v - prior == 1:
            prior = v
        else:
            peak_regions.append(peaks.index.values[start:i + 1])
            start = i + 1
            prior = v

    if len(peak_regions) == 0:
        peak_regions.append(peaks.index.values)
    else:
        # If there are more than one peak region, take only the two with the tallest peak
        for i, p in enumerate(peak_regions):
            peak_regions[i] = [p, find_peak_height(cag, p)]

        peak_regions = sorted(peak_regions, key=lambda x: x[1], reverse=True)

        peak_regions = peak_regions[:2]
        for i, p in enumerate(peak_regions):
            peak_regions[i] = p[0]

    return peak_regions


def find_peak(cag, region):
    # Add a buffer region of 2 to find peaks
    prior = list(range(region[0] - 2, region[0]))
    post = list(range(region[-1] + 1, region[-1] + 3))
    region = prior + list(region) + post

    peaks = {}
    # Find all peaks in this region
    for r in region:
        if cag.loc[r - 1, 'cag_count'] <= cag.loc[r, 'cag_count'] > cag.loc[r + 1, 'cag_count']:
            peaks[cag.loc[r, 'cag_count']] = r

    # If no peaks found take the highest peak as a peak
    if len(peaks) == 0:
        peak_max = max(cag.loc[region, 'cag_count']['cag_count'])
        peak_index = cag.loc[cag.cag_count == peak_max, :].index[0]

        peaks[peak_max] = peak_index

    return peaks


def identify_cag_alleles(cag, peak_regions):
    def check_homozygote():
        # Function to check homozygosity
        p0 = peaks[max(peaks)]
        if cag.loc[p0 + 1, 'cag_count'] / cag.loc[p0, 'cag_count'] >= 0.5:
            # Not homozygote
            c0 = cag.loc[p0, 'cag_size']
            c1 = cag.loc[p0 + 1, 'cag_size']
            return [c0, c1]
        else:
            # Homozygote
            c0 = cag.loc[p0, 'cag_size']
            return [c0, c0]

    # If more than one peak region, find one peak in each region, the maximum peak:
    if len(peak_regions) == 2:
        cag_alleles = []
        for r in peak_regions:
            peaks = find_peak(cag, r)
            peak_loc = peaks[max(peaks)]
            cag_alleles.append(cag.loc[peak_loc, 'cag_size'])

        return cag_alleles
    # if only one region, look for two peaks (less than 50% difference in height)
    # if not, check the cag size next to the identified peak if it is significant enough to consider an allele
    # else classify as homozygote
    elif len(peak_regions) == 1:
        # Check for two peaqks
        peaks = find_peak(cag, peak_regions[0])
        if len(peaks) >= 2:
            # Get the top two peaks present
            p0, p1 = sorted(peaks.keys(), reverse=True)[:2]
            if p1 / p0 >= 0.5:
                c0 = cag.loc[peaks[p0], 'cag_size']
                c1 = cag.loc[peaks[p1], 'cag_size']
                return [c0, c1]
            else:
                return check_homozygote()

        else:
            return check_homozygote()


parser = ap.ArgumentParser(description="Reads through a MiSeq analyzed summarized results file of each sample and "
                                       "extracts the most frequent structures for the two present alleles.")
parser.add_argument("sample", help="Sample name.")
parser.add_argument("cag_dist", help="CAG Distribution of the MiSeq analyzed sample.")
parser.add_argument("res_f", help="Result file with structure frequencies.")
parser.add_argument("output", help="Output prefix.")


args = parser.parse_args()
sample, cag_dist, res_f = args.sample, args.cag_dist, args.res_f

with open(res_f, 'r') as res, open(f"{args.output}.top_alleles.txt", 'w') as out:
    # Skip header
    res.readline()
    res_start = res.tell()
    # Check if file is empty
    if res.readline() == '':
        cag1 = None
        cag2 = None
        r = [None, None, None, None]
        records = [r, r]
    # Check that there's at least 5 records
    elif results_check(res, res_start) is False:
        cag1 = None
        cag2 = None
        r = [None, None, None, None]
        records = [r, r]


    else:
        # Build cag distribution df
        cag = build_cag_dist(args.cag_dist)
        pd.set_option("display.max_rows", cag.shape[0] + 1)

        # Find the peak regions
        try:
            cag = find_peak_regions(cag)
            # Get the peak regions
            regions = get_peak_regions(cag)

            cag_alleles = identify_cag_alleles(cag, regions)
            # Obtain alleles
            cag1, cag2 = sorted(cag_alleles, reverse=True)

        except:
            cag1 = None
            cag2 = None
            r = [None, None, None, None]
            records = [r, r]

        # Iterate over lines, looking for the first structure with the identified expanded CAG size
        if cag1 is not None:
            res.seek(res_start)
            records = []
            for cag in [cag2, cag1]:
                while True:
                    structure = res.readline().strip().split('\t')
                    if structure[1].split(' ')[0] == "CAG%s" % cag:
                        records.append(structure)
                        res.seek(res_start)
                        break

            # For homozygous samples, check the top 10 structures for a CAG size of greater than 30. If found, replace cag1
            if cag1 == cag2:
                i = 0
                while i < 10:
                    structure = res.readline().strip().split('\t')
                    cag_size = int(structure[1].split(' ')[0][3:])
                    if cag_size >= 30:
                        records[1] = structure
                        cag1 = cag_size
                        break
                    else:
                        i += 1

    header = ("sample", "miseq_cag2", "cag2_structure", "cag2_count", "cag2_freq", "miseq_cag1", "cag1_structure",
              "cag1_count", "cag1_freq", "miseq_run")
    metrics = (sample, cag2, *records[0][1:-1], cag1, *records[1][1:])
    out.write("\t".join(header) + '\n')
    out.write("\t".join(["%s"] * len(metrics)) % metrics + '\n')
