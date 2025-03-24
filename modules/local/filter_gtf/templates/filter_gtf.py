#!/usr/bin/env python3

"""Filter GTF file for optimal clipseq execution. Filters GENCODE or ENSEMBL genomic annotation in GTF format."""

import platform
import pandas as pd
import csv
import argparse
import numpy as np
import os


def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")
        out_f.write("    pandas: " + pd.__version__ + "\n")

def ParseIds(annotation_string):
    # Check for NaN or missing values
    if not isinstance(annotation_string, str):
        raise ValueError("Annotation string is missing or NaN.")

    ids = {}
    for id_type in ['gene_id', 'transcript_id']:
        lst = [v for v in annotation_string.split(';') if f'{id_type} ' in v]
        # If the list is empty (gene feature will not have tx id) return np.nan, otherwise clean up the id to remove the id_type and double quotes
        if len(lst) == 0:
            val = np.nan
        elif len(lst) == 1:
            val = lst[0].replace(f'{id_type} ', '').replace('"', '')
        else:
            raise ValueError(f'More than one {id_type} found - corrupted annotation file: {annotation_string}')
        ids[id_type] = val
    return ids['gene_id'], ids['transcript_id']


def read_ann(gtf):
    df_ann = pd.read_csv(gtf,
                            sep='\t',
                            names=['chrom', 'source', 'feature', 'start', 'end', 'name', 'strand', 'name2', 'annotations'],
                            header=None,
                            comment='#',
                            dtype={
                                "chrom": str,
                                "source": str,
                                "feature": str,
                                "start": int,
                                "end": int,
                                "name": str,
                                "strand": str,
                                "name2": str,
                                "annotations": str,
                                })
    ann_cols = df_ann.columns.tolist()
     # Apply ParseIds directly and assign new columns
    df_ann[['gene_id', 'transcript_id']] = df_ann['annotations'].apply(lambda x: pd.Series(ParseIds(x)))
    return df_ann, ann_cols


def Annotate_EnsCan(df_gtf):
    df_gtf['Ensembl_canonical'] = df_gtf['annotations'].apply(lambda x: True if 'tag "Ensembl_canonical"' in x else False)
    # Check if each gene ID has one Ensembl canonical transcript; if not raise error
    all_geneids = set(df_gtf.loc[df_gtf['feature'] == 'gene', 'gene_id'].unique())
    canonicalids = set(df_gtf.loc[(df_gtf['Ensembl_canonical']==True),'gene_id'].unique())
    if all_geneids == canonicalids:
        print("All genes have an Ensembl canonical transcript.")
        return df_gtf
    else:
        print("Some genes do not have one Ensembl canonical transcript.")
        raise ValueError("Some genes do not have one Ensembl canonical transcript.")


def main(process_name, gtf, output=None):

    if output is None:
        output = f"filtered.{os.path.basename(gtf)}"

    # Export versions
    dump_versions(process_name)

    # Parse gtf file into pandas dataframe.
    print("Reading annotation file.")
    input_annotation, ann_cols = read_ann(gtf)

    # Annotate Ensembl canonical transcripts
    print("Annotating Ensembl canonical transcripts.")
    input_annotation = Annotate_EnsCan(input_annotation.copy())

    # Filter
    filt_annot = input_annotation.loc[(input_annotation['Ensembl_canonical'] == True) | (input_annotation['feature'] == 'gene')].copy()

    # Test if all genes have a representative transcript; raise error if not
    minTx = filt_annot.loc[filt_annot.feature.isin(['gene', 'transcript'])].groupby(['gene_id'])['transcript_id'].nunique().min()
    if minTx < 1:
        raise ValueError("Some genes do not have a representative Ensembl_canonical transcript.")
    else:
        print("All genes have a representative Ensembl_canonical transcript")
        # Save filtered annotation
        print("Saving filtered annotation.")
        # filt_annot[ann_cols].to_csv(f"{output}/filtered.{os.path.basename(gtf)}", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)
        filt_annot[ann_cols].to_csv(output, header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)
        return

if __name__ == "__main__":
    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--gtf", default="!{gtf}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    main(args.process_name, args.gtf, args.output)
