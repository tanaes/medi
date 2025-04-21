#!/usr/bin/env Rscript

# download.R: fetch nucleotide contigs only, compress with pigz

library(optparse)
library(data.table)
library(reutils)
library(futile.logger)
library(Biostrings)

option_list <- list(
  make_option(c("--matches"), type="character", help="path to matches.csv"),
  make_option(c("--threads"), type="integer",   default=1, help="number of cores for pigz"),
  make_option(c("--out_dir"), type="character", default="sequences", help="output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

# params
matches_file <- opt$matches
threads      <- opt$threads
out_dir      <- opt$out_dir
rate         <- if (is.null(getOption("reutils.api.key"))) 0.9 else 9

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# function to download and process nucleotide contigs
download_contigs <- function(hits, taxid) {
  fna <- file.path(out_dir, paste0(taxid, ".fna"))
  flog.info("Downloading contigs for taxon %s...", taxid)

  for (i in 0:5) {
    Sys.sleep(1/rate + min(30, 2^i))
    if (file.exists(fna)) unlink(fna)
    post  <- epost(unique(hits$id), db = "nuccore")
    Sys.sleep(1/rate)
    fetch <- suppressMessages(
      efetch(post, db = "nuccore", rettype = "fasta", retmode = "text")
    )
    if (length(getError(fetch)) == 1) {
      write(content(fetch), fna)
      if (file.exists(fna) && grepl(">", content(fetch))) break
    }
  }
  if (!file.exists(fna) || !grepl(">", content(fetch))) {
    stop("Failed downloading contigs for taxid ", taxid)
  }

  # rewrite headers
  fa <- readDNAStringSet(fna)
  short <- tstrsplit(names(fa), "\\s+")[[1]]
  names(fa) <- paste0(short, "_1|kraken:taxid|", taxid, " ", names(fa))
  writeXStringSet(fa, fna, compress = FALSE)

  # compress with pigz
  system2("pigz", c("-p", threads, fna))

  # prepare metadata
  meta <- data.table(
    id = taxid,
    db = "nucleotide",
    filename = paste0(fna, ".gz"),
    num_records = length(fa),
    seqlength = sum(width(fa))
  )
  return(meta)
}

# main
matches <- fread(matches_file)
contigs_out <- NULL
if (any(matches$db == "nucleotide")) {
  contigs_out <- matches[
    db == "nucleotide",
    download_contigs(.SD, matched_taxid[1]),
    by = "matched_taxid"
  ]
  flog.info("Downloaded contigs for %d taxa", nrow(contigs_out))
}

# write metadata for all contigs
if (!is.null(contigs_out)) {
  fwrite(contigs_out, file.path(out_dir, "manifest_contigs.csv"))
} else {
  # write empty metadata file if no contigs
  fwrite(data.table(), file.path(out_dir, "manifest_contigs.csv"))
}