#!/usr/bin/env python3
"""
Count unique reads per BAM that overlap features in the Telescope GTF,
summarised by class_id (repeat elements) or gene_type (ncRNA).

Rules:
 - Each read is counted at most once.
 - If a read's alignments overlap features in MORE THAN ONE category
   (class_id / gene_type), the read is discarded and logged.
 - Reads mapping to exactly one category are counted once for that category.

Uses pybedtools for C-speed interval intersections, then pandas/matplotlib.
"""

import pybedtools
import argparse
import pysam
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import re
import os
from collections import defaultdict

# -- Arguments ----------------------------------------------------------------

parser = argparse.ArgumentParser(
    prog='Count repeat reads',
    description='Counts multi-mapped reads which have been assigned multiple times' \
    'to a single class of repeat element. Discards other reads.'
)

parser.add_argument('--gtf',
                    help = "Annotation for re-assigning multi-mapped reads, gtf format")
parser.add_argument('--bamdir',
                    help = "Directory where input coord-sorted bam files are stored")
parser.add_argument('-o', '--outdir',
                    help = "Directory for program to store output",
                    default = ".")
parser.add_argument('-i', '--input',
                    help = "Comma separated list of coord-sorted bam files",
                    nargs="+")

args = parser.parse_args()

# ── Configuration ────────────────────────────────────────────────────────────
GTF_PATH = args.gtf
BAM_DIR  = args.bamdir
BAM_FILES = args.input
OUTPUT_DIR = args.outdir

# ── Helpers ──────────────────────────────────────────────────────────────────
def extract_sample_name(bam_filename):
    """Sample name - everything before ."""
    return bam_filename.split(".")[0]

def parse_attr(attr_str, key):
    m = re.search(rf'{key}\s+"([^"]+)"', attr_str)
    return m.group(1) if m else None


# ── Ensure BAM index exists ─────────────────────────────────────────────────
def ensure_index(bam_path):
    """Create .bai if missing."""
    bai = bam_path + ".bai"
    if not os.path.exists(bai):
        print(f"  Indexing {os.path.basename(bam_path)} …")
        pysam.index(bam_path)

# ── Ensure BAM is coord-sorted ──────────────────────────────────────────────

def check_sorting(bam):
    sort_order = bam.header.get("HD", {}).get("SO", "unknown")
    if sort_order != "coordinate":
        raise ValueError(
            f"'{bam.filename.decode()}' is '{sort_order}'-sorted, "
            "bam files must be coordinate-sorted."
        )

# ── GTF → sorted BED4 with category ─────────────────────────────────────────
def gtf_to_category_bed(gtf_path):
    """Return a sorted pybedtools BedTool of BED4: chrom start end category.

    category = class_id (TEs) or gene_type (ncRNAs).
    Writes to a temp file to avoid holding 4.8 M intervals in memory.
    """
    print("Building feature BED from GTF …")
    bed_path = os.path.join(OUTPUT_DIR, "_features_unsorted.bed")
    n = 0
    skipped = 0
    with open(gtf_path) as fin, open(bed_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            chrom = cols[0].strip()
            strand = cols[6].strip()
            if not chrom:               # skip lines with empty chromosome
                skipped += 1
                continue
            start = int(cols[3]) - 1    # GTF 1-based → BED 0-based
            end   = int(cols[4])
            attrs = cols[8]

            cat = parse_attr(attrs, "class_id")
            if not cat:
                cat = parse_attr(attrs, "gene_type")
            if not cat:
                continue
            # Replace any whitespace in category with underscore
            cat = cat.replace(" ", "_")
            fout.write(f"{chrom}\t{start}\t{end}\t{cat}\t0\t{strand}\n")
            n += 1

    if skipped:
        print(f"  Skipped {skipped:,} lines with empty chromosome")

    bt = pybedtools.BedTool(bed_path).sort()
    print(f"  {n:,} features loaded")
    return bt


# ── Process one BAM ─────────────────────────────────────────────────────────
def process_bam(bam_path, feature_bed_path):
    """Intersect BAM with feature BED entirely through shell pipeline.

    Phase 1: bedtools bamtobed | sort | bedtools intersect -sorted
             → produces overlap lines in C, fast
    Phase 2: awk collects unique (read, category) pairs per read,
             prints reads with exactly 1 category (counted),
             and dumps multi-category reads to stderr (logged).

    Returns (counts_dict, n_ambiguous).
    """
    import subprocess

    label = os.path.splitext(os.path.basename(bam_path))[0]
    print(f"Processing: {label}")
    ensure_index(bam_path)

    # Build genome file from BAM header for -sorted
    genome_tmp = os.path.join(OUTPUT_DIR, "_genome_tmp.txt")
    import pysam as _pysam
    with _pysam.AlignmentFile(bam_path, "rb") as bam:
        check_sorting(bam)

        chroms = [(sq["SN"], sq["LN"]) for sq in bam.header["SQ"]]
    chroms.sort(key=lambda x: x[0])
    with open(genome_tmp, "w") as f:
        for c, l in chroms:
            f.write(f"{c}\t{l}\n")

    # Entire pipeline in C / awk — no Python iteration
    # intersect output: BED6(read) + BED4(feature) + overlap
    #   $4 = read_name, $10 = category
    #
    # The awk script:
    #   1) For each line, records the read's category set.
    #   2) At END, counts reads with exactly 1 category,
    #      and prints multi-category reads to stderr.
    awk_script = r"""
    {
        rn = $4
        cat = $10
        if (!(rn in cats)) {
            cats[rn] = cat
        } else {
            # append only if new category
            n = split(cats[rn], arr, SUBSEP)
            found = 0
            for (i = 1; i <= n; i++) {
                if (arr[i] == cat) { found = 1; break }
            }
            if (!found) cats[rn] = cats[rn] SUBSEP cat
        }
    }
    END {
        for (rn in cats) {
            n = split(cats[rn], arr, SUBSEP)
            if (n == 1) {
                count[arr[1]]++
            } else {
                # ambiguous → stderr
                printf "%s\t%s\n", rn, cats[rn] > "/dev/stderr"
                ambig++
            }
        }
        for (c in count) print count[c] "\t" c
        printf "AMBIG\t%d\n", ambig+0 > "/dev/stderr"
    }
    """

    with open(os.path.join(OUTPUT_DIR, "_awk_dedup.awk"), "w") as f:
        f.write(awk_script)

    # path for intersected bed
    intersect_path = os.path.join(OUTPUT_DIR, f"{label}_intersected.bed")

    cmd = (
        f"bedtools bamtobed -i '{bam_path}' "
        f"| sort -k1,1 -k2,2n "
        f"| bedtools intersect -a stdin -b '{feature_bed_path}' "
        f"  -wo -sorted -s -g '{genome_tmp}' "
        f"| tee '{intersect_path}' "
        f"| awk -F'\\t' -f '{OUTPUT_DIR}/_awk_dedup.awk'"
    )

    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, executable="/bin/bash",
    )
    if result.returncode != 0 and not result.stdout:
        print(f"  ERROR: {result.stderr.strip()}")
        return {}, 0

    # Parse counts (stdout)
    counts = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t", 1)
        if len(parts) == 2 and parts[0].isdigit():
            counts[parts[1]] = int(parts[0])

    # Parse ambiguous count and reads (stderr)
    ambig_lines = result.stderr.strip().split("\n") if result.stderr.strip() else []
    n_ambiguous = 0
    ambig_reads = []
    for line in ambig_lines:
        if line.startswith("AMBIG\t"):
            n_ambiguous = int(line.split("\t")[1])
        else:
            parts = line.split("\t", 1)
            if len(parts) == 2:
                ambig_reads.append((parts[0], parts[1]))

    # Make bam files of assigned and discarded reads ---------------------------------------------

    # make set for faster lookup
    ambig_reads_set = {rn for rn, _ in ambig_reads}

    # output files for beds
    ambig_bed_path = os.path.join(OUTPUT_DIR, f"{label}_ambiguous.bed")
    assigned_bed_path = os.path.join(OUTPUT_DIR, f"{label}_assigned.bed")

    assigned_reads = []

    # write output files if they don't exist
    with open(intersect_path) as intersected_bed, \
        open(ambig_bed_path, "w") as ambig_out, \
        open(assigned_bed_path, "w") as assigned_out:
        for line in intersected_bed:
            cols = line.rstrip("\n").split("\t")
            intersect_rn = cols[3]
            if intersect_rn in ambig_reads_set:
                ambig_out.write(line)
            else:
                assigned_out.write(line)
                assigned_reads.append(intersect_rn)

    # Make "collapsed" bed file of assigned - one line per read----------------------------------

    # For summarising categories which the reads align to using icount-mini, creating a bed file
    # where each read exists once, keeping only the first occurence of each read.
    # As bucketing assigns reads which only map within the category, it doesn't matter which of the
    # assigned reads icount-mini counts.

    # extract unique read names from list

    assign_reads_dict = dict.fromkeys(assigned_reads, 0)
    collapsed_bed_path = os.path.join(OUTPUT_DIR, f"{label}_collapsed.bed")

    with open(assigned_bed_path) as assigned_out, \
        open(collapsed_bed_path, "w") as collapsed_out:
            for line in assigned_out:
                cols = line.rstrip("\n").split("\t")
                assigned_rn = cols[3]
                for k in assign_reads_dict:
                    if assigned_rn == k:
                        assign_reads_dict[k] = assign_reads_dict[k] + 1
                        if assign_reads_dict[k] == 1:
                            collapsed_out.write(line)
                    else:
                        continue

    # Also export collapsed as a bam for resolve_groups_and_crosslinks workflow

    collapsed_bed = pybedtools.BedTool(collapsed_bed_path)
    collapsed_bam = collapsed_bed.bedtobam(g=genome_tmp)
    #collapsed_bam.saveas(os.path.join(OUTPUT_DIR, f'{label}_collapsed.bam'))

    # and sort and index

    collapsed_sorted_bam_path = os.path.join(OUTPUT_DIR, f'{label}_collapsed.bam')
    pysam.sort("-o", collapsed_sorted_bam_path, collapsed_bam.fn)
    pysam.index(collapsed_sorted_bam_path)

    # -------------------------------------------------------------------------------------------

    total = sum(counts.values())
    print(f"  {total:,} reads kept, {n_ambiguous:,} discarded (multi-category)")

    return counts, ambig_reads


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    feature_bt = gtf_to_category_bed(GTF_PATH)
    # Get path of the sorted BED file for shell pipeline
    feature_bed_path = feature_bt.fn

    all_rows = []
    all_ambiguous = {}

    for bam_file in BAM_FILES:
        sample_name  = extract_sample_name(bam_file)
        bam_path = os.path.join(BAM_DIR, bam_file)
        counts, ambig = process_bam(bam_path, feature_bed_path)

        for cat, n in counts.items():
            all_rows.append({"sample": sample_name, "category": cat, "read_count": n})
        all_ambiguous[sample_name] = ambig

    df = pd.DataFrame(all_rows)

    # ── Save CSV summary ─────────────────────────────────────────────────
    csv_path = os.path.join(OUTPUT_DIR, "repeat_ncRNA_read_counts.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nSummary saved to: {csv_path}")

    pivot = df.pivot_table(
        index="category", columns="sample",
        values="read_count", fill_value=0,
    )
    print(pivot.to_string())

    # ── Log discarded ambiguous reads ────────────────────────────────────
    log_path = os.path.join(OUTPUT_DIR, "discarded_ambiguous_reads.tsv")
    with open(log_path, "w") as f:
        f.write("sample\tread_name\tcategories\n")
        for sample, ambig_list in all_ambiguous.items():
            for rname, cats_str in ambig_list:
                # cats_str uses SUBSEP (0x1c) as separator from awk
                cats = cats_str.replace("\x1c", ",")
                f.write(f"{sample}\t{rname}\t{cats}\n")
    n_total_ambig = sum(len(v) for v in all_ambiguous.values())
    print(f"\n{n_total_ambig:,} ambiguous reads logged to: {log_path}")

    # ── Plot: one facet per sample ──────────────────────────────────────
    samples = sorted(df["sample"].unique())
    n_prot = len(samples)

    fig, axes = plt.subplots(
        1, n_prot,
        figsize=(7 * n_prot, max(6, 0.35 * df["category"].nunique())),
        squeeze=False,
    )

    for i, prot in enumerate(samples):
        ax = axes[0, i]
        sub = df[df["sample"] == prot].sort_values("read_count", ascending=True)
        ax.barh(sub["category"], sub["read_count"], color="steelblue")
        ax.set_xlabel("Unique read count")
        ax.set_title(prot, fontweight="bold")
        ax.tick_params(axis="y", labelsize=8)

    plt.tight_layout()
    for ext in ("pdf", "png"):
        p = os.path.join(OUTPUT_DIR, f"repeat_ncRNA_read_counts.{ext}")
        plt.savefig(p, bbox_inches="tight", dpi=150)
        print(f"Plot saved: {p}")
    plt.close()


if __name__ == "__main__":
    main()
