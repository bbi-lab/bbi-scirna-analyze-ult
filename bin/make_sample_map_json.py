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
# Merge the lane-specific BAM files that belong to a combination
# of sample+process_group+pcr_pair. As a result, the output file
# contains the reads for sample_name, process_group, and pcr_pair.
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
# A sample in a process group can be spread across
# multiple lanes and PCR pairs and so can be in more
# than one entry in the sample_index_list. Therefore
# the entries for the sample+process group must be
# collected into a single entry in the output JSON
# file.
#
# I want the <sample>-<process_group> map to include
#   o  genome
#   o  hash_file
#   o  sample_flags
#
# Make a dictionary as follows
#   data_file_dict[<process_group>][<sample_name>][<pcr_pair_string>][<lane_index_int>] = 1 (a flag)
#
def get_data_file_dict(samplesheet_data, genomes_data):
  sample_names_processed = list()
  data_dict_list = list()
  data_dict_dict = dict()
  for sample_index_dict in samplesheet_data['sample_index_list']:
    # Gather relevant values.
    process_group = sample_index_dict['process_group']
    sample_name = sample_index_dict['sample_id']
    genome = sample_index_dict['genome']
    hash_file = sample_index_dict['hash_file']
    sample_flags = sample_index_dict['sample_flags']

    star_index = genomes_data[genome]['star_index']
    star_memory = genomes_data[genome]['star_memory']
    genes_bed = genomes_data[genome]['genes_bed']

    sample_name_full = '%s-%03d' % (sample_name, int(process_group))
    if(not sample_name_full in sample_names_processed):
      data_dict = {}
      data_dict['sample_name'] = sample_name_full
      data_dict['genome'] = genome
      data_dict['hash_file'] = hash_file
      data_dict['sample_flags'] = sample_flags
      data_dict['star_index'] = star_index
      data_dict['star_memory'] = star_memory
      data_dict['genes_bed'] = genes_bed
      data_dict_list.append(data_dict)
      data_dict_dict[sample_name_full] = data_dict
      sample_names_processed.append(sample_name_full)
    else:
      if(genome != data_dict_dict[sample_name_full]['genome']):
        print('Error: inconsistent genome values for %s' % (sample_name_full), file=sys.stderr)
        sys.exit(1)
      if(hash_file != data_dict_dict[sample_name_full]['hash_file']):
        print('Error: inconsistent hash_file values for %s' % (sample_name_full), file=sys.stderr)
        sys.exit(1)
      if(sample_flags != data_dict_dict[sample_name_full]['sample_flags']):
        print('Error: inconsistent sample_flag values for %s' % (sample_name_full), file=sys.stderr)
        sys.exit(1)
  return(data_dict_list)


#
# Merge the lane-specific BAM files that belong to a combination
# of sample+process_group+pcr_pair. As a result, the output file
# contains the reads for sample_name, process_group, and pcr_pair.
#
def make_data_file_json(data_file_dict_list):
  try:
    filename_json = 'sample_map.json'
    fh = open(filename_json, 'w')
    json.dump(data_file_dict_list, fh, indent=2)
    fh.close()
  except:
    print('Error: unable to write output file \"%s\"' % (filename_json), file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to make JSON file for merging demux bam files.')
  parser.add_argument('-s', '--samplesheet', required=True, default=None, help='Input JSON samplesheet filename (required string).')
  parser.add_argument('-g', '--genomes', required=True, default=None, help='Input genome data JSON filename (required string).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  #
  # Read input samplesheet file.
  #
  samplesheet_data = read_json(args.samplesheet)
  genomes_data = read_json(args.genomes)
  data_file_dict_list = get_data_file_dict(samplesheet_data, genomes_data)
  make_data_file_json(data_file_dict_list)

