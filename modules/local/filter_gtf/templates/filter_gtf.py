#!/usr/bin/env python3

"""Filter GTF file for optimal clipseq execution. Filters GENCODE or ENSEMBL genomic annotation in GTF format."""

import pandas as pd
import csv
import argparse
import numpy as np
import os

def cli():
    parser = argparse.ArgumentParser(description='Filter genomic annotation from GENCODE or ENSEMBL in GTF format by tag \"basic\" and transcript_support_level.'\
    ' These flags are currently (20250312) annotated for Homo sapiens and Mus musculus organisms.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-a', '--annotation', type=str, required=True,
                        help='Annotation file from GENCODE or ENSEMBL in GTF format.')
    required.add_argument('-o', '--outputdir', type=str, required=True,
                        help='Path to output folder.')
    args = parser.parse_args()
    print(args)
    return(args.annotation, args.outputdir)


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


def read_ann(gtf_file):
    df_ann = pd.read_csv(gtf_file,
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
    
        

def filter_gff(gtf_file, outputdir):

    # Parse gtf file into pandas dataframe.
    print("Reading annotation file.")
    input_annotation, ann_cols = read_ann(gtf_file)
    
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
        filt_annot[ann_cols].to_csv(f"{outputdir}/filtered.{os.path.basename(gtf_file)}", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)
        return

def main():
    (gtf_file, outputdir) = cli()
    filter_gff(gtf_file, outputdir)

if __name__ == '__main__':
    main()


# import platform
# import argparse
# from sys import exit
# import csv
# import pandas as pd


# def dump_versions(process_name):
#     with open("versions.yml", "w") as out_f:
#         out_f.write(process_name + ":\n")
#         out_f.write("    python: " + platform.python_version() + "\n")
#         out_f.write("    pandas: " + pd.__version__ + "\n")


# def main(process_name, gtf, output):
#     # Dump version file
#     dump_versions(process_name)

#     # Parse gtf file into pandas dataframe.
#     print("Reading annotation file.")
#     annotation = pd.read_csv(
#         gtf,
#         sep="\t",
#         names=["chrom", "source", "feature", "start", "end", "name", "strand", "name2", "annotations"],
#         header=None,
#         comment="#",
#         dtype={
#             "chrom": str,
#             "source": str,
#             "feature": str,
#             "start": int,
#             "end": int,
#             "name": str,
#             "strand": str,
#             "name2": str,
#             "annotations": str,
#         },
#     )

#     # Filter transcripts by basic tag
#     print("Number of entries in input annotation:", len(annotation))

#     # Check if annotation contains tag "basic"
#     print("Checking for basic flag...")

#     basic = annotation["annotations"].str.contains("basic", regex=True)
#     if basic.any():
#         print("Basic flag available.")

#         nbasic = basic.value_counts()[True]
#         print(f"{nbasic} entries flagged as basic.")

#         annotation = annotation.loc[
#             annotation["annotations"].str.contains('tag "basic"') | (annotation["feature"] == "gene"), :
#         ]
#         print('Number of entries after filtering for tag "basic":', len(annotation))

#         # Filter annotation gene-by-gene by transcript level support (TSL), to keep higher confidence transcripts where possible (TSL1 and 2)
#         df_TSL = annotation.loc[
#             annotation["annotations"].str.contains(
#                 'transcript_support_level "1|transcript_support_level "2', regex=True
#             ),
#             :,
#         ]
#         gene_ids = df_TSL["annotations"].str.split(";", n=1, expand=True)[0].unique().tolist()
#         print("Number of genes that contain TSL1 or TSL2 transcripts:", len(gene_ids))

#         # Keeping only TSL1 and TSL2 entries for genes that contain them, discardig other entries (no TSL information or TSL3-5)
#         print("Filtering out low-confidence transcripts.")
#         df_t = annotation.loc[
#             (annotation["feature"] != "gene") & (annotation["annotations"].str.contains("|".join(gene_ids))), :
#         ]
#         df_t = df_t.loc[
#             ~df_t["annotations"].str.contains('transcript_support_level "1"|transcript_support_level "2"', regex=True)
#         ]
#         annotation.drop(index=df_t.index, inplace=True)

#         print("Number of entries in filtered annotation.", len(annotation))
#         print("Saving filtered gtf file.")
#     else:
#         print('No tag "basic". Returning input annotation as output. Exiting.')
#     annotation.to_csv(output, header=None, index=None, sep="\t", quoting=csv.QUOTE_NONE)


# if __name__ == "__main__":
#     # Allows switching between nextflow templating and standalone python running using arguments
#     parser = argparse.ArgumentParser()
#     parser.add_argument("--process_name", default="!{process_name}")
#     parser.add_argument("--gtf", default="!{gtf}")
#     parser.add_argument("--output", default="!{output}")
#     args = parser.parse_args()

#     main(args.process_name, args.gtf, args.output)
