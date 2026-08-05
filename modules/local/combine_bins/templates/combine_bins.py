#!/usr/bin/env python3
"""
Join the output from both runs of multimap_class_binning - run once at the specific
category level (snRNA | TE | tRNA), and then once at the region level (3UTR | CDS | intron).
"""

import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    prog='Merge assigned reads from binning',
    description='Merges the outputs from the two rounds of binning to produce a .csv for multiqc'
)

parser.add_argument('-a',
                    help = 'Dataframe a, original output')

parser.add_argument('-b',
                    help = 'Dataframe b, secondary output')

args = parser.parse_args()


TE_counts = pd.read_csv(args.a)

region_counts = pd.read_csv(args.b)

merged = pd.merge(TE_counts, region_counts, on="sample", how="left")

merged.to_csv("multimap_binning.csv", index=False)
