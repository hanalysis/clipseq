#!/usr/bin/env python3
"""
Include only reads discarded due to ambiguous TE/ncRNA identity.

Produce a sorted and indexed .bam file from this to be run through
multimap_class_binning again.
"""
import pysam
import pandas as pd
import argparse
import os

# Argparse

parser = argparse.ArgumentParser(
    prog='Create bam with unassigned reads',
    description='Creates bam from unassigned reads - reads which were discarded by previous iteration of binning due to '
    'ambiguity.'
)

parser.add_argument('-b', '--bam',
                    help = 'Original bam file, previously run through multimap_class_binning.py')

parser.add_argument('-d', '--discarded',
                    help = 'List of reads discarded by multimap_class_binning.py due to class' \
                    'ambiguity, in .tsv format')

args = parser.parse_args()

bam = args.bam

# Get ID
label = os.path.splitext(os.path.basename(bam))[0]

## Read in discarded read names

discarded_reads =  pd.read_csv(args.discarded, sep='\t')

discarded_read_names = set(discarded_reads['read_name'])

# Open original bam

orig_bam = pysam.AlignmentFile(bam, "rb")

with pysam.AlignmentFile(f"{label}.discarded_reads.bam", "wb", template=orig_bam) as out:
    for read in orig_bam:
        if read.query_name in discarded_read_names:
            out.write(read)
        else:
            continue

# sort and index

pysam.sort("-o", f"{label}.discarded_reads_sorted.bam", f"{label}.discarded_reads.bam")
pysam.index(f"{label}.discarded_reads_sorted.bam")
