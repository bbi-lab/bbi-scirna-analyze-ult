#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(argparse)
    library(monocle3)
})

parser = argparse::ArgumentParser(description='Script to run emptyDrops.')
parser$add_argument('sample_name', help='Sample name.')
parser$add_argument('mobs', help='Input monocle objects directory.')
parser$add_argument('scrublet_csv', help='Input scrublet CSV file.')
args = parser$parse_args()

sample_name  <- args$sample_name
mobs_dir     <- args$mobs
scrublet_csv <- args$scrublet_csv

directory_out <- paste0(sample_name, '_cds.raw.mobs')
archive_control <- list(archive_type='none')

cds <- load_monocle_objects(mobs_dir)
tryCatch({scrublet_table <- read.csv(scrublet_csv, header=FALSE)
          if(nrow(scrublet_table) != ncol(cds)) {
            stop('mismatch of cell counts in scrublet table and cds')
          }
          pData(cds)$scrublet_score <- as.numeric(scrublet_table$V1)
          save_monocle_objects(cds, directory_path=directory_out, archive_control=archive_control)},
          error = function(e) {
            # Error reading scrublet_csv (most likely) so save cds as it is.
            save_monocle_objects(cds, directory_path=directory_out, archive_control=archive_control)
          })

