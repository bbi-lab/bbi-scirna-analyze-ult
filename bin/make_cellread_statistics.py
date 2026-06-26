#!/usr/bin/env python3

#
# Make read and UMI statistics from counts in the STARsolo CellRead.stats file.
#
# Notes:
#   o  the UMI counts in the CellRead.stats files consist of reads that fall
#      in genomic features, apparently. So use the total_feature_reads to
#      calculate the duplication rate.
#

import sys
import argparse
import statistics
import csv
import json
import re


#
# Program version string.
#
program_version = '0.1.0'


def process_tsv(filename, sample_name, umi_cutoff):
  pobj = re.compile(r'^CBnotInPasslist')
  irow = 0
  sum_genome_unique  = 0.0
  sum_genome_multi   = 0.0
  sum_feature_unique = 0.0
  sum_feature_multi  = 0.0
  sum_counted_unique = 0.0
  sum_counted_multi  = 0.0
  sum_exonic         = 0.0
  sum_intronic       = 0.0
  sum_umi_unique     = 0.0
  sum_umi_multi      = 0.0
  all_umi_list       = list()
  with open(filename, 'r') as fp:
    tsv_reader = csv.DictReader(fp, delimiter='\t')
    for row in tsv_reader:
      if(pobj.match(row['CB']) != None):
         continue
      sum_genome_unique   += float(row['genomeU'])
      sum_genome_multi    += float(row['genomeM'])
      sum_feature_unique  += float(row['featureU'])
      sum_feature_multi   += float(row['featureM'])
      sum_counted_unique  += float(row['countedU'])
      sum_counted_multi   += float(row['countedM'])
      sum_exonic          += float(row['exonic'])
      sum_intronic        += float(row['intronic'])
      sum_umi_unique      += float(row['nUMIunique'])
      sum_umi_multi       += float(row['nUMImulti'])

      all_umi = float(row['nUMIunique']) + float(row['nUMImulti'])
      if(all_umi >= umi_cutoff):
        all_umi_list.append(all_umi)
  if(len(all_umi_list) > 0):
    umi_median = statistics.median(all_umi_list)
  else:
    umi_median = 0.0

  statistics_dict = dict()
  statistics_dict['sample_name']              = sample_name
  statistics_dict['sum_genome_reads_unique']  = int(sum_genome_unique)
  statistics_dict['sum_genome_reads_multi']   = int(sum_genome_multi)
  statistics_dict['sum_feature_reads_unique'] = int(sum_feature_unique)
  statistics_dict['sum_feature_reads_multi']  = int(sum_feature_multi)
  statistics_dict['sum_counted_reads_unique'] = int(sum_counted_unique)
  statistics_dict['sum_counted_reads_multi']  = int(sum_counted_multi)
  statistics_dict['sum_exonic_reads']         = int(sum_exonic)
  statistics_dict['sum_intronic_reads']       = int(sum_intronic)
  statistics_dict['sum_umi_unique']           = int(sum_umi_unique)
  statistics_dict['sum_umi_multi']            = int(sum_umi_multi)
  statistics_dict['total_genome_reads']       = int(sum_genome_unique + sum_genome_multi)
  statistics_dict['total_feature_reads']      = int(sum_feature_unique + sum_feature_multi)
  statistics_dict['total_umi']                = int(sum_umi_unique + sum_umi_multi)
  statistics_dict['umi_median']               = int(umi_median)

  return(statistics_dict)


def make_umi_statistics_json(statistics_dict, filename_json):
  try:
    fh = open(filename_json, 'w')
    json.dump(statistics_dict, fh, indent=2)
    fh.close()
  except:
    print('Error: unable to write output file \"%s\"' % (filename_json), file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program that makes some umi statistics from a STARsolo CellReads.stats file.')
  parser.add_argument('-s', '--sample_name', required=True, default=None, help='Input sample name (required string(s)).')
  parser.add_argument('-c', '--umi_cutoff', required=True, default=None, help='Minimum UMI count (required integer).')
  parser.add_argument('-i', '--input', required=True, default=None, help='Input STARsolo CellReads.stats file (required string(s)).')
  parser.add_argument('-o', '--output', required=True, default=None, help='Output JSON filename (required string(s)).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  sample_name         = args.sample_name
  cell_reads_filename = args.input
  umi_cutoff          = float(args.umi_cutoff)
  filename_json       = args.output

  statistics_dict = process_tsv(cell_reads_filename, sample_name, umi_cutoff)

  make_umi_statistics_json(statistics_dict, filename_json)

