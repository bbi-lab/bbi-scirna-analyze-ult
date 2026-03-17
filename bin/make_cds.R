#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(data.table)
  library(tidyverse)
})

parser = argparse::ArgumentParser(description='Script to make final cds per sample.')
parser$add_argument('sample_name', help='Sample name.')
parser$add_argument('matrix_key', help='Is this \'raw\' or \'filtered\'?')
parser$add_argument('matrix', help='File of umi count matrix.')
parser$add_argument('gene_data', help='File of gene data.')
parser$add_argument('cell_data', help='File of cell data.')
parser$add_argument('barcodes_to_wells', help='File of encoded barcode indices and wells.')
parser$add_argument('umi_counts', help='File of mitochondrial UMI counts.')
parser$add_argument('umi_cutoff', help='UMI cutoff to count as a cell.')
parser$add_argument('counts_per_cell', help='Counts per cell from STARsolo CellReads.stats.')
parser$add_argument('gene_bed', help='Bed file of gene info.')
parser$add_argument('empty_drops', help='RDS file from emptyDrops.')
# parser$add_argument('intron_fraction_file', help='Intron fraction of barcode UMIs file.')
# parser$add_argument('key', help='The sample name prefix.')

args = parser$parse_args()

sample_name <- args$sample_name

umi_cutoff = strtoi(args$umi_cutoff, 10L)

cds <- load_mm_data(mat_path=args$matrix,
                    feature_anno_path=args$gene_data,
                    cell_anno_path=args$cell_data,
                    feature_metadata_column_names=c('gene_short_name', 'gene_expression'),
                    umi_cutoff=umi_cutoff,
                    matrix_control=list(matrix_class='BPCells'))

#
# Drop useless feature values.
#
cds@rowRanges@elementMetadata@listData[['gene_expression']] <- NULL

#
# Add additional gene information.
#
gene_info <- read.csv(args$gene_bed, header=FALSE, sep='\t')
rownames(gene_info) <- gene_info$V4
colnames(gene_info) <- c('chromosome', 'bp1', 'bp2', 'id', 'integer', 'gene_strand')
rowData(cds) <- cbind(rowData(cds), gene_info[rownames(rowData(cds)), c('id', 'chromosome', 'bp1', 'bp2', 'gene_strand')])

#
# Assign percent mitochondrial reads to the cds.
#
counts_per_cell <- fread(args$counts_per_cell,
                         header = TRUE, data.table = F,
                         col.names = c("cell",
                                       "read_total_count",
                                       "umi_count_unique",
                                       "umi_count_multi",
                                       "umi_count_total",
                                       "exonic_read_count",
                                       "intronic_read_count",
                                       "mito_read_count"))

#
# Add percent mitochondrial UMIs.
#
umi_counts <- read.csv(args$umi_counts, header=FALSE, sep='\t')
perc_mito_umis <- umi_counts[3] / (umi_counts[2] + umi_counts[3]) * 100.0
rownames(perc_mito_umis) <- umi_counts$V1
colData(cds)['perc_mitochondrial_umis'] <- perc_mito_umis[rownames(colData(cds)),]

#
# Add number of genes expresses by cell.
#
# cds <- detect_genes(cds)

#
# Add emptyDrops information.
#
emptydrops_data <- readRDS(args$empty_drops)

if(is(emptydrops_data, 'DFrame')) {
  pData(cds)[['emptyDrops_FDR']]         <- emptydrops_data[pData(cds)[,'cell'],]@listData[['FDR']]
  pData(cds)[['emptyDrops_Limited']]     <- emptydrops_data[pData(cds)[,'cell'],]@listData[['Limited']]
  metadata(pData(cds))$emptyDrops_lower  <- metadata(emptydrops_data)[['lower']]
  metadata(pData(cds))$emptyDrops_niters <- metadata(emptydrops_data)[['niters']]
  metadata(pData(cds))$emptyDrops_alpha  <- metadata(emptydrops_data)[['alpha']]
  metadata(pData(cds))$emptyDrops_retain <- metadata(emptydrops_data)[['retain']]
  metadata(pData(cds))$emptyDrops_ignore <- metadata(emptydrops_data)[['ignore']]
  metadata(pData(cds))$emptyDrops_round  <- metadata(emptydrops_data)[['round']]
  ed <- as.data.frame(pData(cds))[,c('cell', 'n.umi', 'emptyDrops_FDR')]
} else {
  ed <- as.data.frame(pData(cds))[,c('cell', 'n.umi')]
}

#
# Add barcode well string to colData(cds).
#
wells <- read.table(args$barcodes_to_wells, sep='\t', row.names=1)
colData(cds)['wells'] <- wells[row.names(colData(cds)),1]

# Extract meta info from well name
df <- as.data.frame(colData(cds))
meta_types <-  meta_types <- c("P5_barcode", "P7_barcode", "RT_barcode", "Ligation_barcode")
meta<- separate(df, wells, into=meta_types, sep="_", remove=FALSE)

for (m in meta_types) {
  colData(cds)[,m] <- meta[[m]]
}

cds <- cds[,Matrix::colSums(counts(cds)) != 0]
cds <- estimate_size_factors(cds)

cds <- detect_genes(cds)

if(ncol(counts(cds)) >= 51) {
  cds <- preprocess_cds(cds)
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds)
  ggp_obj <- suppressMessages(plot_cells(cds))
  file_name <- paste0(sample_name, '_umap.', args$matrix_key, '.png')
  ggsave(filename=file_name, ggp_obj, device='png', width=5, height=5, dpi=600, units='in')
#  saveRDS(cds, file=paste0(sample_name, '_cds.', args$matrix_key, '.rds'))
} else {
  ggp_obj <- ggplot() + geom_text(aes(x = 1, y = 1, label = "Insufficient number of cells")) + monocle3:::monocle_theme_opts() + theme(legend.position = "none") + labs(x="X", y = "Y")
  file_name <- paste0(sample_name, '_umap.', args$matrix_key, '.png')
  ggsave(filename=file_name, ggp_obj, device='png', width=5, height=5, dpi=600, units='in')
}

save_monocle_objects(cds, directory_path=paste0(sample_name, '_cds.', args$matrix_key, '.mobs'), archive_control=list(archive_type='none', archive_compression='none'))

