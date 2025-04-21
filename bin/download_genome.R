#!/usr/bin/env Rscript

# download_genome.R: fetch one GenBank assembly and use seqkit to rewrite headers

library(optparse)
library(futile.logger)
library(data.table)

option_list <- list(
  make_option(c("--id"),    type="character", help="assembly accession"),
  make_option(c("--url"),   type="character", help="ftp path to assembly without filename"),
  make_option(c("--taxid"), type="integer",   help="matched taxid"),
  make_option(c("--out_dir"), type="character", default="sequences", help="output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

# parameters
id      <- opt$id
base_url<- opt$url
taxid   <- opt$taxid
out_dir <- opt$out_dir

# Configure logger for verbose output
flog.threshold(INFO)

# ensure output folder
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# rsync function with verbose debugging
ncbi_rsync <- function(url, out) {
    flog.info("Starting ncbi_rsync with URL: %s", url)
    # Ensure URL has trailing slash for directory listing
    if (!endsWith(url, "/")) {
        base_url <- paste0(url, "/")
    } else {
        base_url <- url
    }
    rsync_base_url <- sub("^https?://", "rsync://", base_url)
    flog.debug("rsync base URL: %s", rsync_base_url)
    
    # Use --list-only to get directory contents
    flog.debug("Running rsync --list-only to inspect directory contents")
    list_cmd <- paste("rsync --no-motd --list-only", rsync_base_url)
    flog.debug("Command: %s", list_cmd)
    
    ret <- system(list_cmd, intern = TRUE)
    flog.debug("Directory listing results (%s lines):", as.character(length(ret)))
    
    # Find the FASTA file by pattern matching
    fna_files <- grep("_genomic\\.fna\\.gz$", ret, value = TRUE)
    
    flog.debug("Found %s genomic FASTA files:", as.character(length(fna_files)))
    
    if (length(fna_files) == 0) {
        flog.error("No _genomic.fna.gz files found in the directory!")
        return(1)
    }
    
    # Extract the filename from the list output
    line_parts <- strsplit(fna_files[1], "\\s+")[[1]]
    fna_file <- tail(line_parts, 1)
    flog.debug("Extracted filename: %s", fna_file)
    
    # Build complete URL with correct filename
    full_url <- paste0(base_url, fna_file)
    rsync_full_url <- sub("^https?://", "rsync://", full_url)
    
    flog.info("Downloading file from: %s", rsync_full_url)
    flog.info("Saving to: %s", out)
    
    # Download the file
    download_cmd <- paste("rsync --no-motd", rsync_full_url, out)
    flog.debug("Download command: %s", download_cmd)
    
    ret_code <- system(download_cmd)
    flog.debug("rsync download return code: %s", as.character(ret_code))
    
    if (ret_code == 0) {
        flog.info("Download successful!")
    } else {
        flog.error("Download failed with code %s", as.character(ret_code))
    }
    
    return(ret_code)
}

# file paths
work_dir <- getwd()
temp_gz_file <- file.path(work_dir, paste0(id, "_temp.fna.gz"))
final_gz_file <- file.path(out_dir, paste0(id, ".fna.gz"))

flog.info("Starting download process for assembly %s", id)
flog.info("Base URL: %s", base_url)
flog.info("Work directory: %s", work_dir)
flog.info("Temporary file: %s", temp_gz_file)
flog.info("Final output file: %s", final_gz_file)

# download with retry (matching original logic)
for (i in 0:7) {
    flog.info("Download attempt %s for %s", as.character(i+1), id)
    
    if (file.exists(temp_gz_file)) {
        flog.debug("Removing existing temporary file: %s", temp_gz_file)
        unlink(temp_gz_file)
    }
    
    # Try to download the file
    ret <- tryCatch(
        ncbi_rsync(base_url, temp_gz_file),
        error = function(e) { 
            flog.error("Error in rsync: %s", e$message)
            return(1)
        },
        warning = function(e) {
            flog.warn("Warning in rsync: %s", e$message)
            return(1)
        }
    )
    
    if (ret == 0) {
        flog.info("Download successful on attempt %s", as.character(i+1))
        break
    }
    
    flog.warn("Attempt %s failed for %s, retrying after delay...", as.character(i+1), id)
    delay <- 2^i
    flog.debug("Waiting %s seconds before retry", as.character(delay))
    Sys.sleep(delay)  # exponential backoff like original
}

if (ret != 0) {
    flog.error("All download attempts failed for %s from %s", id, base_url)
    stop("Failed downloading ", id)
}

# Verify the file was downloaded
flog.info("Verifying downloaded file: %s", temp_gz_file)
if (!file.exists(temp_gz_file)) {
    flog.error("Expected file %s does not exist after download", temp_gz_file)
    stop("File not found after download: ", temp_gz_file)
}

file_info <- file.info(temp_gz_file)
flog.info("Downloaded file size: %s bytes", as.character(file_info$size))

if (file_info$size == 0) {
    flog.error("Downloaded file is empty!")
    stop("Downloaded file is empty: ", temp_gz_file)
}

# Make sure output directory exists
dir.create(dirname(final_gz_file), recursive = TRUE, showWarnings = FALSE)

# Use seqkit to modify headers directly on gzipped files
flog.info("Using seqkit to rewrite FASTA headers")
seqkit_cmd <- sprintf("seqkit replace -p '^>([^ ]+)' -r '>\\${1}_%s|kraken:taxid|%s \\${1}' %s -o %s",
                    "\\i", as.character(taxid), temp_gz_file, final_gz_file)
flog.debug("Seqkit command: %s", seqkit_cmd)
ret_code <- system(seqkit_cmd)

if (ret_code != 0) {
    flog.error("Failed to process headers with seqkit: %s", as.character(ret_code))
    stop("Header processing failed")
}

# Get sequence count for metadata using seqkit stats
flog.info("Counting sequences with seqkit")
count_cmd <- sprintf("seqkit stats -T %s", final_gz_file)
stats <- system(count_cmd, intern = TRUE)
stats_parts <- strsplit(stats[2], "\t")[[1]]
num_records <- as.numeric(stats_parts[4])
seq_length <- as.numeric(stats_parts[5])

# cleanup temporary files
flog.debug("Removing temporary files")
unlink(temp_gz_file)

# write metadata
flog.info("Writing metadata to manifest file")
gen_meta <- data.table(
  id          = id,
  db          = "genbank",
  filename    = final_gz_file,
  num_records = num_records,
  seqlength   = seq_length
)

manifest_file <- file.path(out_dir, paste0("manifest_genome_", id, ".csv"))
flog.debug("Manifest file: %s", manifest_file)
fwrite(gen_meta, manifest_file)

flog.info("Successfully processed genome for assembly %s", id)