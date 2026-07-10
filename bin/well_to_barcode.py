#!/usr/bin/env python3

import sys
import argparse
import re

#
# Read a file of base 4 well index barcode strings and
# convert them to a well name string.
#


# convert encoded index to decimal index
# convert decimal index to well
# allow for setting maximum number of plates
# split cell 'barcode' into 7 base segments
# identify oligo by segment

# encoded barcode indices from rna_rtlig_demux:
#
#    let cell_barcode: String = format!("{:A>7}{:A>7}{:A>7}{:A>7}",
#                                       rt_index_encoded,
#                                       lig_index_encoded,
#                                       p7_index_encoded,
#                                       p5_index_encoded);


#
# max_index is the maximum well index possible in the
# samplesheet/barcode files. So if there are 16 RT utiter
# plates, the maximum well index is 96 * 16.
#
max_index = 96 * 16


def pad_well_col(well_col, zero_pad, id_length):
  if zero_pad:
      template = '%%0%sd' % id_length
  else:
      template = '%s'
  col_id = template % (well_col)
  return col_id


def index_to_well(well_index, across_row_first):
  if(well_index < 0):
    return((0, 'none'))
  nrow = 8
  ncol = 12
  ipl = int( well_index / 96 )
  i96 = well_index - ipl * 96
  if across_row_first:
      well_row = chr(65 + int(i96 / ncol))
      well_col = (i96 % ncol) + 1
  else:
      well_row = chr(65 + (i96 % nrow))
      well_col = int(i96 / nrow) + 1

  well_id = '%s%s' % ( well_row, pad_well_col( well_col, True, 2 ) )

  return( (ipl, well_id ) )


#
# Note: The indexes begin with 1 so well A01 has index 1.
#
# def make_well_to_index_dict(across_row_first):
#   well_index_dict = {}
#   for index in range(1, 97):
#     well_index_dict[index_to_well(index-1, across_row_first)[1]] = index
#   return(well_index_dict)


def base_convert(i, b):
  if(i == 0):
    return([0])
  result = []
  while i > 0:
    result.insert(0, i % b)
    i = i // b
  return result


def make_well_encoder(max_index):
  encoder_base = ['A', 'C', 'G', 'T']

  rt_barcode_dict =  {}
  lig_barcode_dict = {}
  p7_barcode_dict =  {}
  p5_barcode_dict =  {}

  # Include index 0, which represents no barcode: e.g., p5.
  for index in range(0, max_index + 1):
    base4_index = base_convert(index, 4)
    encoded4_digit_list = []
    for j in base4_index:
      encoded4_digit_list.append(encoder_base[j])
      encoded4_index = str.join('',encoded4_digit_list).rjust(7, 'A')
    (ipl, swell) = index_to_well(index-1, True)
    rt_barcode_dict['P%02d-%s' % (ipl+1, swell)] = encoded4_index

    (ipl, swell) = index_to_well(index, True)
    lig_barcode_dict['LIG%d' % (index)] = encoded4_index

    (ipl, swell) = index_to_well(index-1, True)
    p7_barcode_dict['P%02d-%s' % (ipl+1, swell)] = encoded4_index
    (ipl, swell) = index_to_well(index-1, False)
    p5_barcode_dict['P%02d-%s' % (ipl+1, swell)] = encoded4_index

    p7_barcode_dict['none'] = 'AAAAAAA'
    p5_barcode_dict['none'] = 'AAAAAAA'

  return((rt_barcode_dict, lig_barcode_dict, p7_barcode_dict, p5_barcode_dict))


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='A program to convert well names to base 4 encoded barcode index string.')
  parser.add_argument('-i', '--input_file', required=True, help='Input filename.')
  parser.add_argument('-o', '--output_file', required=True, help='Output filename.')
  args = parser.parse_args()

  (rt_barcode_dict, lig_barcode_dict, p7_barcode_dict, p5_barcode_dict) = make_well_encoder(max_index)

  ofh = open(args.output_file, 'w')
  with open(args.input_file, 'r') as ifh:
    for line in ifh:
      wells = line.strip()
      wells_parts = wells.split('_')
      rt_well = wells_parts[2]
      lig_well = wells_parts[3]
      p7_well = wells_parts[1]
      p5_well = wells_parts[0]

      #
      # Allow for PCR well strings that don't have the plate #.
      #
      if(re.match(r'[A-H][0-9]+$', p7_well) != None):
        p7_well = 'P01-%s' % (p7_well)
      if(re.match(r'[A-H][0-9]+$', p5_well) != None):
        p5_well = 'P01-%s' % (p5_well)
     
      error_flag = 0 
      if(not rt_well in rt_barcode_dict):
        print('Error: well_to_barcode.py: rt well index not found in dictionary.\nYou may need to increase the value of max_index in this program.')
        error_flag = 1
      if(not lig_well in lig_barcode_dict):
        print('Error: well_to_barcode.py: ligation well index not found in dictionary.\nYou may need to increase the value of max_index in this program.')
        error_flag = 1
      if(not p7_well in p7_barcode_dict):
        print('Error: well_to_barcode.py: p7 well index not found in dictionary.\nYou may need to increase the value of max_index in this program.')
        error_flag = 1
      if(not p5_well in p5_barcode_dict):
        print('Error: well_to_barcode.py: p5 well index not found in dictionary.\nYou may need to increase the value of max_index in this program.')
        error_flag = 1

      if(error_flag == 1):
        sys.exit(1)

      barcode = ''.join([rt_barcode_dict[rt_well], lig_barcode_dict[lig_well], p7_barcode_dict[p7_well], p5_barcode_dict[p5_well]])
      print('%s_%s_%s_%s\t%s' % (p5_well, p7_well, rt_well, lig_well, barcode), file=ofh)

  ofh.close()

