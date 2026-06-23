#!/usr/bin/env python3

import sys
import os
import json
import argparse
import re

#
# Program version string.
#
program_version = '0.1.0'

# struct SampleMap {
#   sample_id: String,
#   ranges: String,
#   lanes: String,
#   tissue: String,
#   genome: String,
#   hash_file: String,
#   sample_flags: String,
#   external_sample_name: String,
#   wrap_group: String,
#   rt_file: String,
#   ligation_file: String,
#   p7_file: String,
#   p5_file: String,
#   library: String,
#   process_group: String
# }


#
# Read samplesheet json file.
#
def read_json(filename):
  with open(filename, 'r') as fh:
    json_data = json.load(fh)
  return(json_data)


def make_umi_counts_dict(json_data):
  umi_counts_dict = dict()
  for sample_index_dict in json_data['sample_index_list']:
    sample_name = sample_index_dict['sample_id']
    process_group = sample_index_dict['process_group']
    key = '%s-%03d' % (sample_name, int(process_group))
    if(not key in umi_counts_dict):
      tmp_dict = dict()
      tmp_dict['sample_name'] = key
      tmp_dict['in_matrix'] = '%s_counts.raw.matrix.mtx' % key
      tmp_dict['in_features'] = '%s_counts.raw.features.tsv' % key
      tmp_dict['in_barcodes'] = '%s_counts.raw.cells.tsv' % key
      tmp_dict['out_file'] = '%s_umi_counts.tsv' % key
      umi_counts_dict[key] = tmp_dict
  return(umi_counts_dict)


#
# Make JSON file for running umi_counts process.
#
def make_umi_counts_file_json(umi_counts_dict):
  umi_counts_file_list = []
  for sample_name in umi_counts_dict:
    umi_counts_file_list.append(umi_counts_dict[sample_name])

  # print(umi_counts_file_list)

  try:
    filename_json = 'umi_counts.json'
    fh = open(filename_json, 'w+')
    json.dump(umi_counts_file_list, fh, indent=2)
    fh.close()
  except:
    print('Error: unable to write output file \"%s\"' % (filename_json), file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to make JSON file for running umi_counts process.')
  parser.add_argument('-i', '--input', required=True, default=None, help='Input JSON samplesheet filename (required string).')
#  parser.add_argument('-o', '--output', required=False, default=None, help='Output JSON filename (required string).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  #
  # Read input samplesheet file.
  #
  json_data = read_json(args.input)
#  print(json.dumps(json_data['sample_index_list'], indent=2))
  umi_counts_dict = make_umi_counts_dict(json_data)
#  print(umi_counts_dict)
  make_umi_counts_file_json(umi_counts_dict)

