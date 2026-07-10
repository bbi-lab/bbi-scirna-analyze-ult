#!/usr/bin/env python3

#
# Make UMI and cell statistics for exp_dash.
#

import sys
import argparse
import pandas as pd
import statistics
import csv
import json
import math


#
# Program version string.
#
program_version = '0.1.0'


# cd /net/bbi/vol2/home/bge/bbi/tests/nobackup/RNA3-72-a.bbi_scirna.20251117.add_scrublet/work_analyze/86/fa04eeaca18be2df25558b2254e37e
# umi_counts_filename='SeahubZ01-001_umi_counts.tsv'
# empty_drops_fdr_filename='SeahubZ01-001_empty_drops_fdr.tsv'
# umi_cutoff=100


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program that makes umi and cell statistics for the experiment dashboard.')
  parser.add_argument('-s', '--sample_name', required=True, default=None, help='Input sample name (required string(s)).')
  parser.add_argument('-c', '--umi_cutoff', required=True, default=None, help='Minimum UMI count (required integer).')
  parser.add_argument('-f', '--fdr_cutoff', required=True, default=None, help='Maximum empty drops FDR value (required float).')
  parser.add_argument('-u', '--input_umi_counts', required=True, default=None, help='Input umi counts tsv filename (required string(s)).')
  parser.add_argument('-e', '--input_empty_drops_fdr', required=True, default=None, help='Input empty drops RDS filename (required string(s)).')
  parser.add_argument('-o', '--output', required=True, default=None, help='Output JSON filename (required string(s)).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  sample_name              = args.sample_name
  umi_cutoff               = float(args.umi_cutoff)
  fdr_cutoff               = float(args.fdr_cutoff)
  umi_counts_filename      = args.input_umi_counts
  empty_drops_fdr_filename = args.input_empty_drops_fdr
  filename_json            = args.output 

  #
  # Read UMI counts .tsv file into a pandas dataframe.
  # tsv file contents:
  #   <barcode>\t<non-mito_umis>\t<mito_umis>
  #
  barcode_umi_counts_in = pd.read_csv(filepath_or_buffer=umi_counts_filename, sep='\t', header=None, names=['cell', 'umi_non_mito', 'umi_mito'], index_col='cell')

  #
  # Sum non-mito and mito UMIs.
  #
  barcode_umi_counts_sums = barcode_umi_counts_in.sum(axis=1).to_frame(name='umi_all')

  #
  # Remove low UMI barcodes and join with umi_counts_sums.
  # Notes:
  #   o  filter to save memory use
  #   o  filter on umi_cutoff
  #   o  filter on 100 umi cutoff
  #   o  filter on 1000 umi cutoff
  cell_umi_counts_sums = barcode_umi_counts_sums.loc[barcode_umi_counts_sums['umi_all'] >= umi_cutoff]
  cell_umi_counts = cell_umi_counts_sums.join(other=barcode_umi_counts_in, on='cell', how='inner')

  cell_umi_counts_sums_100_umi_cutoff  = barcode_umi_counts_sums.loc[barcode_umi_counts_sums['umi_all'] >= 100]
  cell_umi_counts_sums_1000_umi_cutoff = barcode_umi_counts_sums.loc[barcode_umi_counts_sums['umi_all'] >= 1000]

  #
  # Calculate
  #   total umis (all barcodes)
  #   median umis (filtered barcodes)
  #   median mitochondrial umis (filtered barcodes)
  #   cells (filtered barcodes)
  umis_total                      = barcode_umi_counts_sums.loc[ :, ['umi_all']].sum().iloc[0]
  umis_median                     = cell_umi_counts_sums.median().iloc[0]
  umis_mito_median                = cell_umi_counts.loc[ :, ['umi_mito']].median().iloc[0]
  cell_counts_umi                 = len(cell_umi_counts_sums)
  cell_counts_umi_100_umi_cutoff  = len(cell_umi_counts_sums_100_umi_cutoff)
  cell_counts_umi_1000_umi_cutoff = len(cell_umi_counts_sums_1000_umi_cutoff)


  umis_total = int(umis_total if( not math.isnan(umis_total)) else 0.0)
  umis_median = int(umis_median if( not math.isnan(umis_median)) else 0.0)
  umis_mito_median = int(umis_mito_median if( not math.isnan(umis_mito_median)) else 0.0)

  #
  # Calculate cells with empty cells FDR <= .01
  # Set cell counts for FDR <= .01 to -1 when the
  # the input_empty_drops_fdr file has no rows;
  # that is, emptyDrops did not run.
  #
  empty_drops_fdr = pd.read_csv(filepath_or_buffer=empty_drops_fdr_filename, sep='\t', header=0, index_col='cell')

  if(len(empty_drops_fdr) > 0):
    barcode_umi_counts = barcode_umi_counts_sums.join(other=barcode_umi_counts_in, on='cell', how='inner')
    barcode_umi_counts_joined = barcode_umi_counts.join(other=empty_drops_fdr, on='cell', how='left')
    cell_umi_counts_fdr = barcode_umi_counts_joined.loc[barcode_umi_counts_joined['FDR'] <= fdr_cutoff]
    cell_counts_fdr = len(cell_umi_counts_fdr)
  else:
    cell_counts_fdr = -1

  #
  # Set up output dictionary.
  #
  umi_cell_counts_dict = dict()
  umi_cell_counts_dict['sample_name']                     = sample_name
  umi_cell_counts_dict['umi_cutoff']                      = umi_cutoff
  umi_cell_counts_dict['empty_drops_fdr_cutoff']          = fdr_cutoff
  umi_cell_counts_dict['umis_total']                      = umis_total
  umi_cell_counts_dict['umis_median']                     = umis_median
  umi_cell_counts_dict['umis_mito_median']                = umis_mito_median
  umi_cell_counts_dict['cell_counts_umi']                 = cell_counts_umi
  umi_cell_counts_dict['cell_counts_fdr']                 = cell_counts_fdr
  umi_cell_counts_dict['cell_counts_umi_100_umi_cutoff']  = cell_counts_umi_100_umi_cutoff
  umi_cell_counts_dict['cell_counts_umi_1000_umi_cutoff'] = cell_counts_umi_1000_umi_cutoff

  #
  # Write output JSON file.
  #
  print('write file %s' % (filename_json), file=sys.stderr)
  try:
    fh = open(filename_json, 'w')
    json.dump(umi_cell_counts_dict, fh, indent=2)
    fh.close()
  except:
    print('Error: unable to write output file \"%s\"' % (filename_json), file=sys.stderr)
    sys.exit(1)

