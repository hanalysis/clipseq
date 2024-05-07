#!/usr/bin/env Rscript

# args <- commandArgs(trailingOnly = TRUE)

# Check if we have exactly three arguments
# if (length(args) != 3) {
#   cat("Usage: script.R <list_of_tsvs> <gtf_file> <output_name>\n")
#   quit(status = 1)
# }

# Load necessary libraries
if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
    !requireNamespace("rtracklayer", quietly = TRUE)) {
  stop("Please install the 'GenomicRanges' and 'rtracklayer' packages.")
}

# Assigning arguments to variables for easier access
# list_of_tsvs <- args[1]
# gtf_file <- args[2]
# output_name <- args[3]

list_of_tsvs <- "!{bedtools_map_outputs}"
gtf_file <- "!{gtf}"
output_name <- "!{output_name}"

# Read the list of TSV files
tsv_files <- strsplit(list_of_tsvs, " ")[[1]]

# Function to read and preprocess each TSV
process_tsv <- function(file) {
  df <- read.table(file, header = TRUE, sep = "\t", quote = "")
  df$PeakID <- apply(df[, 1:4], 1, paste, collapse = "_")
  count_column <- df[, 7, drop = FALSE]
  colnames(count_column) <- sub("\\.[^.]*$", "", basename(file))  # Rename column to the base name of the file
  df <- data.frame(PeakID = df$PeakID, count_column, chr = df[, 1], start = df[, 2], end = df[, 3], strand = df[, 4])
  
  return(df)
}


# Apply the function to each file and merge them
# Process each TSV file into a list of data frames
processed_dfs <- lapply(tsv_files, process_tsv)

# Extract fixed genomic columns only from the first data frame
genomic_data <- processed_dfs[[1]][c("PeakID", "chr", "start", "end", "strand")]

# Function to remove unwanted genomic columns from subsequent data frames
strip_genomic_info <- function(df) {
  df[, setdiff(names(df), c("chr", "start", "end", "strand"))]
}

# Extract only the count columns from all data frames
count_data_frames <- lapply(processed_dfs, strip_genomic_info)

# Merge all count data frames by 'PeakID', which is now the only common column in each
merged_counts <- Reduce(function(x, y) merge(x, y, by = "PeakID", all = TRUE), count_data_frames)

# Combine the genomic data with the merged count data
merged_data <- merge(genomic_data, merged_counts, by = "PeakID", all = TRUE)


# Read the GTF file as a TSV, extract the gene_id from column 9
gtf_data <- read.table(gtf_file, comment.char = "#", stringsAsFactors = FALSE, sep = "\t", header = FALSE)
# Extract gene_id from the attributes field in column 9
gtf_data$gene_id <- sapply(gtf_data$V9, function(x) {
  # Find matches for 'gene_id' followed by any characters up to the first semicolon
  matches <- regmatches(x, regexpr('gene_id ([^;]+);', x))
  if (length(matches) > 0) {
    sub('gene_id ([^;]+);.*', '\\1', matches)
  } else {
    NA  # Return NA if no gene_id is found
  }
})


#print(head(gtf_data$gene_id))

# Create a GRanges object from GTF data
gene_ranges <- GenomicRanges::GRanges(seqnames = gtf_data$V1,
                       ranges = IRanges::IRanges(start = gtf_data$V4, end = gtf_data$V5),
                       strand = gtf_data$V7,
                       gene_id = gtf_data$gene_id)

# Create a GRanges object from the merged TSV data
tsv_ranges <- GenomicRanges::GRanges(seqnames = merged_data$chr,
                      ranges = IRanges::IRanges(start = merged_data$start, end = merged_data$end),
                      strand = merged_data$strand,
                      PeakID = merged_data$PeakID)

# Find overlaps and merge GeneIDs
if (length(tsv_ranges) > 0 && length(gene_ranges) > 0) {
  overlaps <- GenomicRanges::findOverlaps(tsv_ranges, gene_ranges, ignore.strand = FALSE)
  if (length(overlaps) > 0) {
    merged_data$GeneID <- NA  # Initialize GeneID column
    merged_data$GeneID[S4Vectors::queryHits(overlaps)] <- gene_ranges$gene_id[S4Vectors::subjectHits(overlaps)]
  } else {
    cat("No overlaps found between TSV data and gene annotations.\n")
  }
} else {
  cat("One of the range objects is empty. Check your input data.\n")
}

merged_data <- merged_data[c("PeakID", "GeneID", setdiff(names(merged_data), c("PeakID", "GeneID", "chr", "start", "end", "strand")))]
# Write the output file
write.table(merged_data, file=output_name, sep="\t", quote=FALSE, row.names=FALSE)
cat("Merging complete, output saved to", output_name, "\n")

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
rtracklayer.version <- as.character(packageVersion('rtracklayer'))
genomicranges.version <- as.character(packageVersion('GenomicRanges'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    rtracklayer:', rtracklayer.version),
        paste('    GenomicRanges:', genomicranges.version)
    ),
'versions.yml')