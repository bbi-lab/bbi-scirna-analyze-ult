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
# Prepare a lists of untrimmed BAM files that belong to
# a combination of sample+process group. The list is
# submitted to the process_hashes program.
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
      print('Error: expand_index_list: bad index specification: %s' % (index_spec), file=sys.stderr)
      sys.exit(1)
    index1 = int(mobj.group(1))
    index2 = index1
    if(mobj.group(2) != None):
      index2 = int(mobj.group(3))
    for one_index in range(index1, index2+1):
      index_list.append(one_index)
  return(list(set(index_list)))


#
# Make a dictionary of hash_file values keyed
# by sample_name+process_group. Check for consistency
# across all sample_index_list entries.
#
def make_sample_hash_dict(json_data):
  sample_hash_dict = {}
  for sample_index_dict in json_data['sample_index_list']:
    process_group = sample_index_dict['process_group']
    sample_name = sample_index_dict['sample_id']
    hash_file = sample_index_dict['hash_file']
    if(sample_hash_dict.get(process_group) == None):
      sample_hash_dict[process_group] = {}
    if(sample_hash_dict[process_group].get(sample_name) == None):
      sample_hash_dict[process_group][sample_name] = []
    sample_hash_dict[process_group][sample_name].append(hash_file)

  # Check that the values in each list are the same.
  for process_group in sample_hash_dict.keys():
    for sample_name in sample_hash_dict[process_group].keys():
      hash_file_list = sample_hash_dict[process_group][sample_name]
      if(not all(x==hash_file_list[0] for x in hash_file_list)):
        print('Error: inconsistent hash_file_names for sample \'%s\'' % (sample_name), file=sys.stderr)
        sys.exit(1)

  return(sample_hash_dict)


#
# Get PCR data from samplesheet JSON file contents.
#
# json_data['sample_index_list']['lanes']
# json_data['sample_index_list']['ranges']
# json_data['sample_index_list']['p7_file']
# json_data['sample_index_list']['p5_file']
#
# A sample in a process group can be spread across
# multiple lanes and PCR pairs and so can be in more
# than one entry in the sample_index_list.
# The lanes and pcr pairs are gathered together
# for the sample+process group, regardless of whether
# the lane and pcr pairs appear together in any of
# the sample_index_list entries.
#
def get_data_file_dict(json_data, sample_hash_dict):
  data_file_dict = {}
  for sample_index_dict in json_data['sample_index_list']:
    index_ranges = sample_index_dict['ranges'].split(':')
    process_group = sample_index_dict['process_group']
    sample_name = sample_index_dict['sample_id']
    p7_index_list = expand_index_list(index_ranges[1])
    p5_index_list = expand_index_list(index_ranges[2])

    #
    # Consider sample_name+process_group if there is a
    # hash_file value for it (with non-zero length).
    #
    if(sample_hash_dict.get(process_group) != None
       and sample_hash_dict[process_group].get(sample_name) != None
       and len(sample_hash_dict[process_group][sample_name][0]) > 0):

      if(not process_group in data_file_dict):
        data_file_dict[process_group] = {}

      if(not sample_name in data_file_dict[process_group]):
        data_file_dict[process_group][sample_name] = {}

      for p7_index in p7_index_list:
        for p5_index in p5_index_list:
          pcr_pair = '%03d_%03d' % (int(p7_index), int(p5_index))
          if(not pcr_pair in data_file_dict[process_group][sample_name]):
            data_file_dict[process_group][sample_name][pcr_pair] = {}
          data_file_dict[process_group][sample_name][pcr_pair] = 1
  return(data_file_dict)


def make_data_file_json(data_file_dict, sample_hash_dict):
  hash_tsv_list = []
  for process_group in data_file_dict.keys():
    for sample_name in data_file_dict[process_group].keys():
      if(sample_hash_dict.get(process_group) != None
        and sample_hash_dict[process_group].get(sample_name) != None
        and len(sample_hash_dict[process_group][sample_name][0]) > 0):

        for pcr_pair in data_file_dict[process_group][sample_name].keys():
          in_file = '%s-%03d_%s.hash_reads.merged.tsv' % (sample_name, int(process_group), pcr_pair)
          hash_file = sample_hash_dict[process_group][sample_name][0]
          out_root = '%s-%03d_%s' % (sample_name, int(process_group), pcr_pair)
          merge_dict = {}
          merge_dict['sample_name'] = '%s-%03d' % (sample_name, int(process_group))
          merge_dict['in_file'] = in_file
          merge_dict['hash_file'] = hash_file
          merge_dict['out_root'] = out_root
          hash_tsv_list.append(merge_dict)

  try:
    filename_json = 'process_hashes.json'
    fh = open(filename_json, 'w')
    json.dump(hash_tsv_list, fh, indent=2)
    fh.close()
  except:
    print('Error: unable to write output file \"%s\"' % (filename_json), file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to make JSON file for setting process_hash runs.')
  parser.add_argument('-i', '--input', required=True, default=None, help='Input JSON samplesheet filename (required string).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  #
  # Read input samplesheet file.
  #
  json_data = read_json(args.input)
  sample_hash_dict = make_sample_hash_dict(json_data)
  data_file_dict = get_data_file_dict(json_data, sample_hash_dict)

  make_data_file_json(data_file_dict, sample_hash_dict)



