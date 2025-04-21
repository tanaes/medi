#!/usr/bin/env Rscript

library(optparse)
library(data.table)
library(futile.logger)

option_list <- list(
  make_option(c("--contigs_manifest"), type="character", help="path to contigs manifest file"),
  make_option(c("--manifests_dir"), type="character", help="directory containing genome manifests")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Get parameters
contigs_manifest <- opt$contigs_manifest
manifests_dir <- opt$manifests_dir

# Configure logger
flog.threshold(INFO)

flog.info("Loading contigs manifest: %s", contigs_manifest)
contigs <- fread(contigs_manifest)

flog.info("Looking for genome manifests in: %s", manifests_dir)
genome_files <- list.files(manifests_dir, pattern='^manifest_genome_.*\\.csv$', full.names=TRUE)
flog.info("Found %d genome manifest files", length(genome_files))

# Load and combine all manifests
if (length(genome_files) > 0) {
  # Read all manifests and combine them
  genomes_list <- lapply(genome_files, function(f) {
    tryCatch({
      dt <- fread(f)
      return(dt)
    }, error = function(e) {
      flog.error("Error reading file: %s", f)
      flog.error("Error message: %s", as.character(e))
      return(NULL)
    })
  })
  
  # Filter out NULL entries (failed reads)
  genomes_list <- genomes_list[!sapply(genomes_list, is.null)]
  
  if (length(genomes_list) > 0) {
    genomes <- rbindlist(genomes_list, fill=TRUE)
    flog.info("Loaded %d genome records", nrow(genomes))
  } else {
    genomes <- data.table()
    flog.warn("No genome data could be loaded")
  }
} else {
  genomes <- data.table()
  flog.warn("No genome manifest files found")
}

# Combine contigs and genomes
full <- rbind(contigs, genomes, fill=TRUE)
flog.info("Combined manifest contains %d records", nrow(full))

# Write the combined manifest
output_file <- "manifest.csv"
fwrite(full, output_file)
flog.info("Combined manifest written to %s", output_file)