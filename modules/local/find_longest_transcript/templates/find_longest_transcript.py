#!/usr/bin/env python3

"""Calculates a table of CDS lengths for each protein coding transcript and selects the largest one."""

import platform
import argparse
from sys import exit
import pyranges as pr
import pandas as pd
import warnings
import csv

def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")


def main(process_name, gtf, user_transcripts, output):
    # Dump version file
    dump_versions(process_name)
    print(type(user_transcripts))
    # Parse gtf
    df_gtf = pr.read_gtf(gtf, as_df=True)
    df_gtf['length'] = df_gtf['End'] - df_gtf['Start']
    input_genes = set(df_gtf.gene_id.unique())
    print("number of total genes in: ", len(input_genes))
    # if not any(g in df_gtf.columns for g in g_types):
    #     raise ValueError(f"Gene type column not specified with {g_types}")

    # Calculate CDS and exon lengths per transcript ID
    cds_sums = df_gtf.loc[df_gtf['Feature'] == 'CDS'].groupby(['gene_id', 'transcript_id'])['length'].sum()
    exon_sums = df_gtf.loc[df_gtf['Feature'] == 'exon'].groupby(['gene_id', 'transcript_id'])['length'].sum()
    # Join the two dataframes
    df_txlengths = pd.concat([cds_sums, exon_sums], axis=1, keys=['cds_length', 'exon_length']).reset_index()

    # # Check if all transcripts have exon length
    # if df_txlengths.cds_length.isna().any():
    #     warnings.warn("Some transcripts do not have exon length.", RuntimeWarning)
    # Fill in nan with 0
    df_txlengths[['cds_length', 'exon_length']] = df_txlengths[['cds_length', 'exon_length']].fillna(0)

    # If user specified transcripts are provided, use these to filter gtf; else find longest CDS / exon transcript.
    if user_transcripts != "":
        with open(user_transcripts, "r") as file:
            transcript_ids = [line.replace('\n', '') for line in file]

        print(user_transcripts)
        print(transcript_ids[1:5])

        # Check if all user provided ids are in the gtf
        if not set(transcript_ids).issubset(set(df_txlengths.transcript_id.unique())):
            raise ValueError("Some user provided transcript ids are not in the gtf.")

        # Filter based on user provided transcript id
        df_txlengths_filtered = df_txlengths.loc[df_txlengths.transcript_id.isin(transcript_ids)]

        # Check if all genes have exactly one user defined tx id
    else:
        # Sort
        df_txlengths = df_txlengths.sort_values(by=['cds_length', 'exon_length', 'transcript_id'], ascending=[False, False, True])
        # Drop duplicates on gene_id, keep first row
        df_txlengths_filtered = df_txlengths.drop_duplicates(subset='gene_id', keep='first')
        transcript_ids = df_txlengths_filtered.transcript_id.unique().tolist()

    # Check if all of the initial genes have a transcript assigned
    if df_txlengths_filtered.transcript_id.isna().any() or df_txlengths_filtered.transcript_id.isna().any():
        raise ValueError("Some genes do not have a representative transcript.")

    set_filtgenes = set(df_txlengths_filtered.gene_id.unique())
    print(len(set_filtgenes))
    print(len(input_genes))

    if set_filtgenes != input_genes:
        raise ValueError("Some genes are missing after filtering for cds and exon lengths.")

    # Filter
    filt_annot = df_gtf.loc[(df_gtf['transcript_id'].isin(transcript_ids)) | (df_gtf['Feature'] == 'gene')].copy()
    print(filt_annot.dtypes)
    # Convert the filtered annotation to pyranges
    pr_gtf = pr.PyRanges(filt_annot[[c for c in filt_annot.columns if c!='length']])
    # Save the filtered annotation
    pr_gtf.to_gtf(f"{output}.filtered.gtf")

    # Get a list of transcript ids for saving
    pr_gtf = pr_gtf.as_df()

    transcript_ids_sorted = []
    for g in pr_gtf.loc[pr_gtf.Feature=='transcript', 'transcript_id'].values.tolist():
        if g not in transcript_ids_sorted:
            transcript_ids_sorted.append(g)

    # Save to file
    print("Saving representative transcript per gene...")
    with open(output + ".txt", "w") as f:
        f.write("\n".join(map(str, transcript_ids_sorted)) + "\n")

    # Make a transcript fai file
    # Set transcript_id as index
    df_txlengths_filtered.set_index('transcript_id', inplace=True)
    # Sort in order of transcript_ids
    df_txlengths_filtered = df_txlengths_filtered.loc[transcript_ids_sorted]
    # Reset index
    df_txlengths_filtered.reset_index(inplace=True)

    # Save transcript id and exon length to .fai file, without header
    df_txlengths_filtered[['transcript_id', 'exon_length']].to_csv(f"{output}.fai", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)

    # Make the transcript gtf file
    df_txlengths_filtered['src'] = 'src'
    df_txlengths_filtered['gene'] = 'gene'
    df_txlengths_filtered['start'] = 1
    df_txlengths_filtered['score'] = '.'
    df_txlengths_filtered['strand'] = '+'
    df_txlengths_filtered['frame'] = '.'
    df_txlengths_filtered['attributes'] = df_txlengths_filtered.apply(lambda row: f'gene_id:{row["gene_id"]}; transcript_id:{row["transcript_id"]}', axis=1)
    # Order into gtf format
    df_txlengths_filtered[['transcript_id', 'src', 'gene', 'start', 'exon_length', 'score', 'strand', 'frame', 'attributes']].to_csv(f"{output}.gtf", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)
    return

        # else:
        #     # Get the gene type column that is in the gtf
        #     gene_type_col = [v for v in df_gtf.columns if v in g_types][0]
        #     print(gene_type_col)

        #     # transcript_type_col = [v for v in df_gtf.columns if v in t_types][0]

        #     # FEATURETODO include genes encoding the antibody chains as protein coding (ie. IG_C_gene)

        #     gtf_prot = df_gtf.loc[df_gtf[gene_type_col] == "protein_coding"].copy()
        #     input_prot = set(gtf_prot.gene_id.unique())
        #     print("number of protein coding genes in: ", len(input_prot))
        #     gtf_noncod = df_gtf.loc[df_gtf[gene_type_col] != "protein_coding"].copy()
        #     input_noncod = set(gtf_noncod.gene_id.unique())
        #     print("number of protein non-coding genes in: ", len(input_noncod))

    # Filtering protein-coding gtf; only keep protein_coding transcript types for protein_coding genes
    # gtf_prot = gtf_prot.loc[(gtf_prot['Feature'] == 'gene') | (gtf_prot[transcript_type_col] == 'protein_coding')]

    # Can come after gtf filter stage
    # check_genes = set(gtf_prot.gene_id.unique())
    # minTx = gtf_prot.loc[gtf_prot['Feature'] == 'transcript'].groupby('gene_id')['transcript_id'].nunique().min()
    # if (check_genes != input_prot) or (minTx < 1):
    #     raise ValueError(f"Protein coding gene IDs do not match after filtering or some genes do not have a representative transcript.")


    # # Concatenate the two dataframes
    # df_txlengths = pd.concat([df_txlengths_prot, df_txlengths_nc], ignore_index=True)

    # # Extract transcript ids
    # transcript_ids = df_txlengths.transcript_id.tolist()


if __name__ == "__main__":
    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--gtf", default="!{gtf}")
    parser.add_argument("--user_transcripts", default="!{user_transcripts}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    main(args.process_name, args.gtf, args.user_transcripts, args.output)