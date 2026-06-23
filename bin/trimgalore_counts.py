#!/usr/bin/env python3

import sys
import argparse
import re
import json


#
# Program version string.
#
program_version = '0.1.0'


def initialize_counter_dict(sample_name, file_list):
  counter_dict = dict()
  counter_dict['sample_name'] = sample_name
  counter_dict['count_source'] = 'trimgalore log files' 
  counter_dict['input_files'] = file_list
  counter_dict['total_reads_in'] = 0
  counter_dict['total_reads_short'] = 0
  return(counter_dict)


def process_trimgalore_log(file_name, counter_dict):
  pobj_dict = dict()
  pobj_dict['total_reads_in'] = re.compile('([0-9]+) sequences processed in total')
  pobj_dict['total_reads_short'] = re.compile('Sequences removed because they became shorter than the length cutoff of [0-9]+ bp:[ \t]+([0-9]+) [(][0-9.]+%[)]')
#  print('file: %s' % (file_name))
  with open(file_name, 'r') as fp:
    tmp_counter_dict = dict()
    for line in fp:
      for key in pobj_dict.keys():
        mobj = pobj_dict[key].match(line)
        if(mobj != None):
          counter_dict[key] += int(mobj.group(1))
  counter_dict['total_reads_out'] = counter_dict['total_reads_in'] - counter_dict['total_reads_short']


def make_trimgalore_json(counter_dict, filename_json):
  try:
    fh = open(filename_json, 'w')
    json.dump(counter_dict, fh, indent=2)
  except:
    print('Error: unable to write output file \"%s\"' % (filename_json), file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to gather counts from trimgalore log files.')
  parser.add_argument('-s', '--sample_name', required=True, default=None, help='Input sample name (required string(s)).')
  parser.add_argument('-i', '--input', required=True, default=None, nargs='+', help='Input trimgalore log files (required string(s)).')
  parser.add_argument('-o', '--output', required=True, default=None, help='Output JSON filename (required string(s)).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  sample_name = args.sample_name
  file_list = args.input
  filename_json = args.output
 
  counter_dict = initialize_counter_dict(sample_name, file_list)
  for file in file_list:
    process_trimgalore_log(file, counter_dict)

  make_trimgalore_json(counter_dict, filename_json)

#  print('total reads in: %d' % (counter_dict['total_reads_in']))
#  print('total reads short: %d' % (counter_dict['total_reads_short']))
#  print('total reads out:   %d' % (counter_dict['total_reads_in'] - counter_dict['total_reads_short']))
