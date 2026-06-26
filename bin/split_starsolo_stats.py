#!/usr/bin/env python3

#
# Read CellRead.stats file and split out columns required
# by downstream processes.
#


import sys
import argparse
import csv
import os.path
import re


def process_file(in_filename, sample_name):

  pobj = re.compile(r'^CBnotInPasslist')

  out_filename = '%s_counts_per_cell.txt' % (sample_name)
  ofh_umi_unique = open(file=out_filename, mode='w')
  csv_writer = csv.writer(ofh_umi_unique, delimiter='\t', lineterminator='\n')

  cbIndex = None
  cbMatchIndex = None
  nUMIuniqueIndex = None
  nUMImultiIndex = None
  exonicIndex = None
  intronicIndex = None
  mitoIndex = None

  with open(file=in_filename, mode='r', newline='') as ifh:
    csv_reader = csv.reader(ifh, delimiter='\t')

    header_row = next(csv_reader)

    # Find desired columns.
    nUMIuniqueIndex = None
    exonicIndex = None
    intronicIndex = None
    mitoIndex = None
    for i, item in enumerate(header_row):
      if(item == 'CB'):
        cbIndex = i
      elif(item == 'cbMatch'):
        cbMatchIndex = i
      elif(item == 'nUMIunique'):
        nUMIuniqueIndex = i
      elif(item == 'nUMImulti'):
        nUMImultiIndex = i
      elif(item == 'exonic'):
        exonicIndex = i
      elif(item == 'intronic'):
        intronicIndex = i
      elif(item == 'mito'):
        mitoIndex = i

    if(cbIndex == None or
       cbMatchIndex == None or
       nUMIuniqueIndex == None or
       nUMImultiIndex == None or
       exonicIndex == None or
       intronicIndex == None or
       mitoIndex == None):
      print('Error: missing at least one column in the list \'cbIndex\', \'cbMatchIndex\', \'nUMIuniqueIndex\', \'nUMImultiIndex\', \'exonicIndex\', \'intronicIndex\', \'mitoIndex\'', file=sys.stderr)
      sys.exit(1)


    # Write header row.
    #
    # Note: any modifications to the columns here must be matched in assign_hash.R (counts_per_cell = fread(args$umis_per_cell,...)
    #
    csv_writer.writerow(['cell', 'cbMatch', 'nUMIunique', 'nUMImulti', 'nUMItot', 'exonic', 'intronic', 'mito'])

    # Write data rows.
    for row in csv_reader:
      if(pobj.match(row[0]) != None):
        continue
      csv_writer.writerow([row[cbIndex], row[cbMatchIndex], row[nUMIuniqueIndex], row[nUMImultiIndex], int(row[nUMIuniqueIndex]) + int(row[nUMImultiIndex]), row[exonicIndex], row[intronicIndex], row[mitoIndex]])

  ofh_umi_unique.close()


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='A program to split columns out of STARsolo CellReads.stats file.')
  parser.add_argument('-i', '--input_file', required=True, help='CellReads.stats file to processs.')
  parser.add_argument('-s', '--sample_name', required=True, help='Name of sample.')
  args = parser.parse_args()

  process_file(args.input_file, args.sample_name)

