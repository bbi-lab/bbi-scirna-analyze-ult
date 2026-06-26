#!/usr/bin/env python3

import sys
import json
import argparse
import re
from collections import defaultdict

#
# Program version string.
#
program_version = '0.3.0'

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
# Get PCR data from samplesheet JSON file contents.
#
# json_data['sample_index_list']['lanes']
# json_data['sample_index_list']['ranges']
# json_data['sample_index_list']['p7_file']
# json_data['sample_index_list']['p5_file']

def iter_pcr_pairs(ranges):
  fields = ranges.split(':')
  if len(fields) != 3:
    raise ValueError('bad ranges field: %s' % ranges)

  p7_list = expand_index_list(fields[1])
  p5_list = expand_index_list(fields[2])

  for p7 in p7_list:
    for p5 in p5_list:
      yield '%03d_%03d' % (p7, p5)


def get_merge_map(json_data):
  merge_map = defaultdict(set)

  for row in json_data['sample_index_list']:
    sample_id = row['sample_id']
    process_group = row['process_group']
    lanes = expand_index_list(row['lanes'])

    for pcr_pair in iter_pcr_pairs(row['ranges']):
      key = (sample_id, process_group, pcr_pair)
      merge_map[key].update(lanes)

  return merge_map

def make_data_file_json(merge_map, bam_path):
  merge_jobs = []

  for (sample_id, process_group, pcr_pair), lanes in sorted(merge_map.items()):
    sample_full = '%s-%03d' % (sample_id, int(process_group))

    merge_jobs.append({
      'sample_name': sample_full,
      'out_file': '%s_%s.merged.bam' % (sample_full, pcr_pair),
      'in_file_list': [
        '%s/%s-%03d_%s-L%03d.bam' %
        (bam_path, sample_id, lane, pcr_pair, lane)
        for lane in sorted(lanes)
      ],
    })

  with open('merge_demux.json', 'w') as fh:
    json.dump(merge_jobs, fh, indent=2)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to make JSON file for merging demux bam files.')
  parser.add_argument('-i', '--input', required=True, default=None, help='Input JSON samplesheet filename (required string).')
  parser.add_argument('-p', '--bam_path', required=True, default=None, help='Input BAM file path (required string).')
#  parser.add_argument('-o', '--output', required=False, default=None, help='Output JSON output filename (required string).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  #
  # Read input samplesheet file.
  #
  json_data = read_json(args.input)
#  print(json.dumps(json_data['sample_index_list'], indent=2))
#  data_file_dict = get_data_file_dict(json_data)

#  print(json.dumps(data_file_dict, indent=2))

#  make_data_file_json(data_file_dict, args.bam_path)

  merge_map = get_merge_map(json_data)
  make_data_file_json(merge_map, args.bam_path)

