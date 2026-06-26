#!/usr/bin/env python3

#
# Concatenate CSV (tab-delimited) files.
#
# Notes:
#   o  check for header consistency
#   o  remove headers from file 2+
#   o  files are unsorted
#   o  the CellReads.stats files have a row with CB
#      value 'CBnotInPasslist'. Modify the name to
#      make them distinct. I am inclined to skip
#      this line because it will need to removed
#      for some downstream uses.
#


import sys
import argparse
import csv
import os.path


def process_files(in_filenames, out_filename, delimiter ):
  ofh = open(out_filename, 'w')
  csv_writer = csv.writer(ofh, delimiter=delimiter, lineterminator='\n')

  for ifile, filename in enumerate(in_filenames):
    with open(file=filename, mode='r', newline='') as ifh:
      csv_reader = csv.reader(ifh, delimiter=delimiter)

      header_row = next(csv_reader)
      if(ifile != 0):
       if(not header_row == ref_header_row):
         print('Error: input CSV file headers differ.', file=sys.stderr)
         sys.exit(1)
      else:
        ref_header_row = header_row
        csv_writer.writerow(ref_header_row)

      for row in csv_reader:
        if(row[0] == 'CBnotInPasslist'):
          basename = os.path.basename(filename)
          row[0] = row[0] + '.' + basename
        csv_writer.writerow(row)

  ofh.close()


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='A program to concatenate STARsolo CellReads.stats files.')
  parser.add_argument('-i', '--input_files', required=True, nargs='+', help='CSV files to concatenate.')
  parser.add_argument('-o', '--output_file', required=True, help='Name of output CSV file.')
  args = parser.parse_args()

  process_files(args.input_files, args.output_file, delimiter='\t')


