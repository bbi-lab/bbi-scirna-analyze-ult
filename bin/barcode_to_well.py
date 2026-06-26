#!/usr/bin/env python3

import sys
import argparse


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


def make_barcode_decoder(max_index):
  encoder_base = ['A', 'C', 'G', 'T']

  rt_well_dict =  {}
  lig_well_dict = {}
  p7_well_dict =  {}
  p5_well_dict =  {}

  # Include index 0, which represents no barcode: e.g., p5.
  for index in range(0, max_index + 1):
    base4_index = base_convert(index, 4)
    encoded4_digit_list = []
    for j in base4_index:
      encoded4_digit_list.append(encoder_base[j])
      encoded4_index = str.join('',encoded4_digit_list).rjust(7, 'A')
    (ipl, swell) = index_to_well(index-1, True)
    rt_well_dict[encoded4_index] = 'P%02d-%s' % (ipl+1, swell)

    (ipl, swell) = index_to_well(index, True)
    lig_well_dict[encoded4_index] = 'LIG%d' % (index)
 
#     if(index <= 96):
#       (ipl, swell) = index_to_well(index-1, True)
#       p7_well_dict[encoded4_index] = 'P%02d-%s' % (ipl+1, swell)
#       (ipl, swell) = index_to_well(index-1, False)
#       p5_well_dict[encoded4_index] = 'P%02d-%s' % (ipl+1, swell)

    (ipl, swell) = index_to_well(index-1, True)
    if(swell != 'none'):
      p7_well_dict[encoded4_index] = 'P%02d-%s' % (ipl+1, swell)
    else:
      p7_well_dict[encoded4_index] = '%s' % (swell)
    (ipl, swell) = index_to_well(index-1, False)
    if(swell != 'none'):
      p5_well_dict[encoded4_index] = 'P%02d-%s' % (ipl+1, swell)
    else:
      p5_well_dict[encoded4_index] = '%s' % (swell)

  return((rt_well_dict, lig_well_dict, p7_well_dict, p5_well_dict))


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='A program to convert base 4 encoded barcode index string to string of well names.')
  parser.add_argument('-i', '--input_file', required=True, help='Input filename.')
  parser.add_argument('-o', '--output_file', required=True, help='Output filename.')
  args = parser.parse_args()

  (rt_well_dict, lig_well_dict, p7_well_dict, p5_well_dict) = make_barcode_decoder(max_index)

  ofh = open(args.output_file, 'w')
  with open(args.input_file, 'r') as ifh:
    for line in ifh:
      barcode = line.strip()
      rt_code = barcode[0:7]
      lig_code = barcode[7:14]
      p7_code = barcode[14:21]
      p5_code = barcode[21:]

      error_flag = 0
      if(not rt_code in rt_well_dict):
        print('Error: barcode_to_well.py: rt well index not found in dictionary.\nYou may need to increase the value of max_index in this program.')
        error_flag = 1
      if(not lig_code in lig_well_dict):
        print('Error: barcode_to_well.py: lig well index not found in dictionary.\nYou may need to increase the value of max_index in this program.')
        error_flag = 1
      if(not p7_code in p7_well_dict):
        print('Error: barcode_to_well.py: p7 well index not found in dictionary.\nYou may need to increase the value of max_index in this program.')
        error_flag = 1
      if(not p5_code in p5_well_dict):
        print('Error: barcode_to_well.py: p5 well index not found in dictionary.\nYou may need to increase the value of max_index in this program.')
        error_flag = 1

      if(error_flag == 1 ):
        sys.exit(1)

      print('%s\t%s_%s_%s_%s' % (barcode, p5_well_dict[p5_code], p7_well_dict[p7_code], rt_well_dict[rt_code], lig_well_dict[lig_code]), file=ofh)

  ofh.close()

