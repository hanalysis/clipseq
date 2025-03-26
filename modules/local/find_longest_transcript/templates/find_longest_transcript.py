#!/usr/bin/env python3

"""
Calculates a table of CDS and exon lengths for each transcript and filters
genome annotation by representative transcripts.

Handles both user-provided transcripts and automatic selection of the longest
transcript per gene based on CDS and exon length.

Outputs:
- A filtered genome annotation GTF containing only entries containing representative transcripts
- A list of selected transcript IDs (.txt)
- A transcript length index file (.fai) with transcript ID and exon length
- A transcript GTF
"""

import platform
import argparse
from sys import exit
import pyranges as pr
import pandas as pd
import warnings
import csv
import logging


def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")

def parse_gtf_and_calculate_lengths(gtf_path):
    """
    Reads a GTF file and calculates CDS and exon lengths per transcript.

    Returns:
    - df_gtf: full parsed GTF as a DataFrame (with added 'length' column)
    - df_txlengths: DataFrame with gene_id, transcript_id, cds_length, exon_length
    - input_genes: set of all gene_ids present in the GTF
    """
    df_gtf = pr.read_gtf(gtf_path, as_df=True)
    df_gtf['length'] = df_gtf['End'] - df_gtf['Start']
    input_genes = set(df_gtf.gene_id)

    # Calculate CDS and exon lengths per transcript
    cds_sums = df_gtf[df_gtf['Feature'] == 'CDS'].groupby(['gene_id', 'transcript_id'])['length'].sum()
    exon_sums = df_gtf[df_gtf['Feature'] == 'exon'].groupby(['gene_id', 'transcript_id'])['length'].sum()
    tx_sums = df_gtf[df_gtf['Feature'] == 'transcript'].copy().set_index(['gene_id', 'transcript_id'])['length']
    df_txlengths = pd.concat([cds_sums, exon_sums, tx_sums], axis=1, keys=['cds_length', 'exon_length', 'unspliced_length']).reset_index()

    # Warn if any transcript is missing exon length
    missing_exons = df_txlengths[df_txlengths.exon_length.isna()]
    if not missing_exons.empty:
        warnings.warn("Some transcripts do not have exon length.", RuntimeWarning)
        print("Transcripts missing exon length:")
        print(missing_exons[['gene_id', 'transcript_id']])

    # Fill missing values with 0
    df_txlengths[['cds_length', 'exon_length', 'unspliced_length']] = df_txlengths[['cds_length', 'exon_length', 'unspliced_length']].fillna(0)

    # Check if all genes are present for filtering
    missing_genes = input_genes - set(df_txlengths.gene_id)
    if len(missing_genes) > 0:
        logging.error("Some genes don't have transcripts associated with them")
        raise ValueError("FATAL: Some genes don't have transcripts associated with them. Corruped GTF file.")

    return df_gtf, df_txlengths, input_genes


def main(process_name, gtf, user_transcripts, output):
    # Dump version file
    dump_versions(process_name)

    # Logging
    log_file = f"{output}.log"
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Load the genome annotation file and calculate CDS and exon lengths per transcript
    df_gtf, df_txlengths, input_genes = parse_gtf_and_calculate_lengths(gtf)
    logging.info(f"Total genes in genome GTF: {len(input_genes)}")

    # If user specified transcripts are provided, use these to filter gtf; else find longest CDS / exon transcript.
    if user_transcripts != "":
        logging.info("Using the transcript IDs in: ", user_transcripts)
        with open(user_transcripts, "r") as file:
            transcript_ids = [line.strip() for line in file]
            # Remove potential empty strings caused by trailing newline
            transcript_ids = [x for x in transcript_ids if x]

        # Check if all user-provided IDs are in the GTF
        missing_ids = set(transcript_ids) - set(df_txlengths.transcript_id)
        if len(missing_ids) > 0:
            logging.warning(f"The following transcript IDs are not in the GTF: {sorted(missing_ids)}")
            raise ValueError("Some user-provided transcript IDs are not in the GTF.")


        # Filter based on user provided transcript ID
        df_txlengths_filtered = df_txlengths.loc[df_txlengths.transcript_id.isin(transcript_ids)]

        # Throw an error if none of the user transcripts are in the GTF
        if df_txlengths_filtered.empty:
            logging.error("FATAL: None of the provided transcript IDs are found in the GTF.")
            raise ValueError("FATAL: None of the provided transcript IDs are found in the GTF.")

    else:
        # Select longest transcript per gene: hierarchy: CDS length > exon length; remaining ties are resolved alphabetically on transcript ID
        # Sort
        df_txlengths = df_txlengths.sort_values(by=['cds_length', 'exon_length', 'unspliced_length', 'transcript_id'], ascending=[False, False, False, True])
        # Drop duplicates on gene_id, keep first row
        df_txlengths_filtered = df_txlengths.drop_duplicates(subset='gene_id', keep='first')
        transcript_ids = df_txlengths_filtered.transcript_id.unique().tolist()

    # Check if all genes have exactly one matching transcript ID

    # Count number of transcripts per gene
    transcripts_per_gene = df_txlengths_filtered.groupby("gene_id").transcript_id.nunique()
    # Raise warning if any gene has not exactly one transcript and print problematic genes
    multi_assigned_genes = transcripts_per_gene[transcripts_per_gene > 1].index.tolist()
    if not (transcripts_per_gene == 1).all():
        logging.error("Some genes do not have exactly one representative transcript. Offending genes (more than 1 transcript assigned): ", multi_assigned_genes)
        raise ValueError("FATAL: Some genes do not have exactly one representative transcript.")

    set_filtgenes = set(df_txlengths_filtered.gene_id)
    logging.info(f"Remaining genes after filtering: {len(set_filtgenes)}")

    missing_genes = input_genes - set_filtgenes
    if len(missing_ids) > 0:
        missing_genes_list = sorted(missing_genes)
        top = missing_genes_list[:10]
        more = f" and {len(missing_genes_list) - 10} more" if len(missing_genes_list) > 10 else ""

        logging.error(f"Missing genes after filtering: {top}{more}")
        raise ValueError("FATAL: Some genes are missing after filtering. Ensure each has one representative transcript.")

    # Filter the genome annotation
    filt_annot = df_gtf.loc[(df_gtf['transcript_id'].isin(transcript_ids)) | (df_gtf['Feature'] == 'gene')].copy()
    # Convert the filtered annotation to pyranges
    pr_gtf = pr.PyRanges(filt_annot[[c for c in filt_annot.columns if c!='length']])
    # Save the filtered annotation
    logging.info(f"Saving genome GTF filtered by representative transcripts to {output}_filtered.gtf")
    pr_gtf.to_gtf(f"{output}_filtered.gtf")

    # Get a list of transcript IDs for saving
    pr_gtf = pr_gtf.as_df()
    transcript_ids_sorted = []
    for g in pr_gtf.loc[pr_gtf.Feature=='transcript', 'transcript_id'].values.tolist():
        if g not in transcript_ids_sorted:
            transcript_ids_sorted.append(g)

    # Save transcript IDs to file
    logging.info(f"Saving representative transcript IDs to {output}.txt")
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
    logging.info(f"Saving representative transcript fai to {output}.fai")
    df_txlengths_filtered[['transcript_id', 'exon_length']].to_csv(f"{output}.fai", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)

    # Create the transcript GTF file
    df_txlengths_filtered['src'] = 'src'
    df_txlengths_filtered['gene'] = 'gene'
    df_txlengths_filtered['start'] = 1
    df_txlengths_filtered['score'] = '.'
    df_txlengths_filtered['strand'] = '+'
    df_txlengths_filtered['frame'] = '.'
    df_txlengths_filtered['attributes'] = df_txlengths_filtered.apply(lambda row: f'gene_id:{row["gene_id"]}; transcript_id:{row["transcript_id"]}', axis=1)

    logging.info(f"Saving representative transcript GTF to {output}.gtf")
    # Order into gtf format
    df_txlengths_filtered[['transcript_id', 'src', 'gene', 'start', 'exon_length', 'score', 'strand', 'frame', 'attributes']].to_csv(f"{output}.gtf", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)
    return


if __name__ == "__main__":
    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--gtf", default="!{gtf}")
    parser.add_argument("--user_transcripts", default="!{user_transcripts}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    main(args.process_name, args.gtf, args.user_transcripts, args.output)