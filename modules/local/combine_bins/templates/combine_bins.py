#!/usr/bin/env python3
"""
Join the output from both runs of multimap_class_binning - run once at the specific
category level (snRNA | TE | tRNA), and then once at the region level (3UTR | CDS | intron).
"""

import pandas as pd
import argparse
import re, os

parser = argparse.ArgumentParser(
    prog='Merge assigned reads from binning',
    description='Creates bam from unassigned reads - reads which were discarded by previous iteration of binning due to '
    'ambiguity.'
)

parser.add_argument('-a',
                    help = 'Dataframe a, original output')

parser.add_argument('-b',
                    help = 'Dataframe b, secondary output')

parser.add_argument('-premap',
                    nargs = "+",
                    help = "List of bowtie logs from pre-mapping")

args = parser.parse_args()


TE_counts = pd.read_csv(args.a)

region_counts = pd.read_csv(args.b)

## Create df from bowtie logs

linenum = 0
d = []

for logfile in args.premap:
    sample_id = os.path.basename(logfile).replace(".out", "").replace("_ncrna", "")
    with open(logfile, "rt") as myfile:
        for line in myfile:
            linenum += 1
            if line.find("Reported ") != -1 :
                first = line.find("Reported ") + 9
                second = line.find(" alignments")
                val = line[first:second]
                d.append({"sample":sample_id, "pre-mapped":int(val)})
            else:
                continue

df = pd.DataFrame(d)

print(df)

merged = pd.merge(TE_counts, region_counts, on="sample", how="left")

print(merged)

merged_withpremap = pd.merge(merged, df, on = "sample", how = "left")

print(merged_withpremap)

merged_withpremap.to_csv("multimap_binning.csv", index=False)