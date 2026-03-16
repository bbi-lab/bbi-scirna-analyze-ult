#!/usr/bin/env python3

import sys
import os
import argparse
import csv


#
# Program version string.
#
program_version = '0.1.0'

def read_summary_file(filename):
  with open(filename, 'r', newline='') as ifh:
    summary_dict = dict()
    csv_reader = csv.reader(ifh, dialect='excel')
    for row in csv_reader:
      if(re_int.match(row[1])):
        summary_dict[row[0]] = int(row[1])
      elif(re_float.match(row[1])):
        summary_dict[row[0]] = float(row[1])
      else:
        summary_dict[row[0]] = row[1]
  return(summary_dict)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to gather STARsolo Features.stats files.')
  parser.add_argument('-i', '--input', required=True, default=None, nargs='+', help='Input TSV file paths (required string).')
  parser.add_argument('-o', '--output', required=True, default=None, help='Output TSV filename (required stgring).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  in_file_list = args.input
  out_file = args.output

  #
  # Begin logging.
  #
  lfh = open('merge_hash_reads.log', 'w')
  print('program: merge_hash_reads', file=lfh)
  print('version: %s' % (program_version), file=lfh)
  print('output file: %s' % (out_file), file=lfh)
  print('input file(s):', file=lfh)
  for file in in_file_list:
    print('  %s' % (file), file=lfh)

  #
  # Merge TSV files.
  #
  num_reads = 0
  with open(out_file, 'w', newline='') as ofh:
    csv_writer = csv.writer(ofh, delimiter='\t', escapechar='\\', quoting=csv.QUOTE_NONE, lineterminator='\n')
    header_ref = list()
    header_write_flag = False
    for i in range(len(in_file_list)):
      with open(in_file_list[i], 'r', newline='') as ifh:
        csv_reader = csv.reader(ifh, delimiter='\t', escapechar='\\', quoting=csv.QUOTE_NONE)
        try:
          header = next(csv_reader)
          print('header: ', header)
          print('header len: %d' % (len(header)))
          if(len(header) == 0):
            print('Error: missing header.')
            sys.exit(-1)
          if(len(header_ref) != 0):
            if(header != header_ref):
              print('Error: inconsistent file headers.', file=sys.stderr)
              sys.exit(-1)
          else:
            if(header_write_flag == True):
              print('Error: inconsistent condition.')
              sys.exit(-1)
            header_ref = header
            csv_writer.writerow(header)
            header_write_flag = True
        except StopIteration:
          continue
        for row in csv_reader:
          csv_writer.writerow(row)
          num_reads += 1

  print('number of reads written: %d' % (num_reads), file=lfh)
  lfh.close()
