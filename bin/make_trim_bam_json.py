#!/usr/bin/env python3

import sys
import json
import argparse
import re

#
# Program version string.
#
program_version = '0.1.0'

#
# Gather input BAM files and output filenames for running
# the bbduk adaptor trimming program.
#

#
# Read samplesheet json file.
#
def read_json(filename):
  with open(filename, 'r') as fh:
    json_data = json.load(fh)
  return(json_data)


#
# Expand index list.
# Example: 1-4,8-12
# Return a list of distinct indices.
regex_pattern = r'([0-9]+)([-]([0-9]+))?$'
def expand_index_list(index_string):
  index_list = []
  for index_spec in index_string.split(','):
    mobj = re.match(regex_pattern, index_spec)
    if(mobj == None):
      print('Error: expand_index_list: bad index specification: %s' % (index_spec))
      sys.exit(1)
    index1 = int(mobj.group(1))
    index2 = index1
    if(mobj.group(2) != None):
      index2 = int(mobj.group(3))
    for one_index in range(index1, index2+1):
      index_list.append(one_index)
  return(list(set(index_list)))


#
# Get PCR data from samplesheet JSON file contents.
#
# json_data['sample_index_list']['lanes']
# json_data['sample_index_list']['ranges']
# json_data['sample_index_list']['p7_file']
# json_data['sample_index_list']['p5_file']
#
def get_data_file_dict(json_data):
  data_file_dict = {}
  for sample_index_dict in json_data['sample_index_list']:
    index_ranges = sample_index_dict['ranges'].split(':')
    process_group = sample_index_dict['process_group']
    sample_name = sample_index_dict['sample_id']
    p7_index_list = expand_index_list(index_ranges[1])
    p5_index_list = expand_index_list(index_ranges[2])
    lane_index_list = expand_index_list(sample_index_dict['lanes'])

    if(not process_group in data_file_dict):
      data_file_dict[process_group] = {}

    for p7_index in p7_index_list:
      for p5_index in p5_index_list:
        pcr_pair = '%03d_%03d' % (int(p7_index), int(p5_index))
        if(not pcr_pair in data_file_dict[process_group]):
          data_file_dict[process_group][pcr_pair] = {}
        if(not sample_name in data_file_dict[process_group][pcr_pair]):
          data_file_dict[process_group][pcr_pair][sample_name] = {}
        for lane_index in lane_index_list:
          data_file_dict[process_group][pcr_pair][sample_name][lane_index] = 1
  return(data_file_dict)


def make_data_file_json(data_file_dict):
  trim_bam_list = []
  for process_group in data_file_dict.keys():
    for pcr_pair in data_file_dict[process_group].keys():
      for sample_name in data_file_dict[process_group][pcr_pair].keys():
        merge_dict = {}
        root_file = '%s-%03d_%s' % (sample_name, int(process_group), pcr_pair)
        in_file = '%s.merged.bam' % (root_file)
        out_file = '%s.trimmed.bam' % (root_file)
        merge_dict['sample_name'] = '%s-%03d' % (sample_name, int(process_group))
        merge_dict['root_file'] = root_file
        merge_dict['in_file'] = in_file
        merge_dict['out_file'] = out_file
        trim_bam_list.append(merge_dict)

  try:
    filename_json = 'trim_bam.json'
    fh = open(filename_json, 'w')
    json.dump(trim_bam_list, fh, indent=2)
  except:
    print('Error: unable to write output file \"%s\"' % (filename_json), file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to make JSON file for setting trim BAMs runs.')
  parser.add_argument('-i', '--input', required=True, default=None, help='Input JSON samplesheet filename (required string).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  #
  # Read input samplesheet file.
  #
  json_data = read_json(args.input)
  data_file_dict = get_data_file_dict(json_data)

  make_data_file_json(data_file_dict)



