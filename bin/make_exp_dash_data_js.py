#!/usr/bin/env python3

import sys
import re
import argparse
import json
import math


#
# Program version string.
#
program_version = '0.2.0'


#
# Read json file.
#
def read_json(filename):
  with open(filename, 'r') as fh:
    json_data = json.load(fh)
  return(json_data)


def get_sample_name_list(sample_map):
  sample_name_list = list()
  for sample_dict in sample_map:
    sample_name = sample_dict['sample_name']
    sample_name_list.append(sample_name)
  return(sample_name_list)


def get_barnyard_sample_name(sample_map):
  barnyard_sample_name = None
  for sample_dict in sample_map:
    sample_flags = sample_dict['sample_flags']
    if('B' in sample_flags):
      sample_name = sample_dict['sample_name']
      barnyard_sample_name = sample_name
      break
  return(barnyard_sample_name)


def read_cellread_statistics(sample_name_list):
  cellread_statistics_dict = dict()
  for sample_name in sample_name_list:
    filename = '%s_cellreads_statistics.json' % sample_name
    stats_json = read_json(filename)
    cellread_statistics_dict[sample_name] = stats_json
  return(cellread_statistics_dict)


def read_starsolo_summary(sample_name_list):
  pobj_list    = list()
  for i in range(7):
    new_list = ['', None, int()]
    pobj_list.append(new_list)
  pobj_list[0][0] = 'number_of_reads'
  pobj_list[0][1] = re.compile('Number of Reads[ \t]+([0-9]+)')
  pobj_list[0][2] = int
  pobj_list[1][0] = 'reads_with_valid_barcodes'
  pobj_list[1][1] = re.compile('Reads With Valid Barcodes[ \t]+([0-9.]+)')
  pobj_list[1][2] = float
  pobj_list[2][0] = 'sequencing_saturation'
  pobj_list[2][1] = re.compile('Sequencing Saturation[ \t]+([0-9.]+)')
  pobj_list[2][2] = float
  pobj_list[3][0] = 'reads_mapped_to_genome_unique_and_multiple'
  pobj_list[3][1] = re.compile('Reads Mapped to Genome: Unique[+]Multiple[ \t]+([0-9.]+)')
  pobj_list[3][2] = float
  pobj_list[4][0] = 'reads_mapped_to_genome_unique'
  pobj_list[4][1] = re.compile('Reads Mapped to Genome: Unique[ \t]+([0-9.]+)')
  pobj_list[4][2] = float
  pobj_list[5][0] = 'reads_mapped_to_genefull_ex50pas_unique_and_multiple'
  pobj_list[5][1] = re.compile('Reads Mapped to GeneFull_Ex50pAS: Unique[+]Multiple GeneFull_Ex50pAS[ \t]+([0-9.]+)')
  pobj_list[5][2] = float
  pobj_list[6][0] = 'reads_mapped_to_genefull_ex50pas_unique'
  pobj_list[6][1] = re.compile('Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS[ \t]+([0-9.]+)')
  pobj_list[6][2] = float

  starsolo_summary_dict = dict()
  for sample_name in sample_name_list:
    filename = '%s_Summary.txt' % sample_name
    with open(filename, 'r') as ifh:
      starsolo_summary = dict()
      for line in ifh:
        for pobj in pobj_list:
          mobj = pobj[1].match(line.strip())
          if(mobj != None):
            starsolo_summary[pobj[0]] = pobj[2](mobj.group(1))
            starsolo_summary_dict[sample_name] = starsolo_summary
            break
  return(starsolo_summary_dict)


def read_umi_cell_statistics(sample_name_list):
  umi_cell_statistics_dict = dict()
  for sample_name in sample_name_list:
    filename = '%s_umi_cell.json' % sample_name
    stats_json = read_json(filename)
    umi_cell_statistics_dict[sample_name] = stats_json
  return(umi_cell_statistics_dict)


# [
#   ...
#   {
#     "sample_name": "SeahubZ01-001",
#     "genome": "Zebrafish",
#     "hash_file": "/net/bbi/vol2/home/bge/bbi/tests/nobackup/RNA3-72-a.bbi_scirna_analyze.20260303.exp_dash_hash_read_rate/timecourse_hash2.txt",
#     "sample_flags": "",
#     "star_index": "/net/bbi/vol1/data/genomes_stage/zebrafish/zebrafish_star_2_7",
#     "star_memory": 40,
#     "genes_bed": "/net/bbi/vol1/data/genomes_stage/zebrafish/zebrafish_star_2_7/latest.genes.bed"
#   },
#   ...
# ]
#
#
# Total reads: 8663052
# Hash reads: 8430665
# Hash rate: 0.9732
#
def read_hash_read_rates(sample_map_list):
  pobj_list    = list()
  for i in range(3):
    new_list = ['', None, int()]
    pobj_list.append(new_list)

  pobj_list[0][0] = 'total_reads'
  pobj_list[0][1] = re.compile('Total reads:[ \t]+([0-9]+)')
  pobj_list[0][2] = int
  pobj_list[1][0] = 'hash_reads'
  pobj_list[1][1] = re.compile('Hash reads:[ \t]+([0-9.]+)')
  pobj_list[1][2] = int
  pobj_list[2][0] = 'hash_rate'
  pobj_list[2][1] = re.compile('Hash rate:[ \t]+([0-9.]+)')
  pobj_list[2][2] = float

  hash_read_rate_dict = dict()
  for sample_map in sample_map_list:
    sample_name = sample_map['sample_name']
    if(sample_map['hash_file'] != None):
      hash_file = sample_map['hash_file']
      hash_read_rate = dict()
      if(len(hash_file) > 0):
        filename = '%s_hash_read_rate.txt' % (sample_name)
        with open(filename, 'r') as ifh:
          for line in ifh:
            for pobj in pobj_list:
              mobj = pobj[1].match(line.strip())
              if(mobj != None):
                hash_read_rate[pobj[0]] = pobj[2](mobj.group(1))
                hash_read_rate_dict[sample_name] = hash_read_rate
                break
      hash_read_rate_dict[sample_name] = hash_read_rate
  return(hash_read_rate_dict)


#
# const run_data =
# {
#   "run_name": "RNA3-72-a.canonical_sci_rna.20251117",
#   "cell_counts": ["64882", "-"],
#   "sample_list": [
#     "Sentinel",
#     "24.024",
#     "Fishbowl",
#     "Keyhole",
#     "SeahubZ01"
#   ],
#   "barn_collision": null,
#   "sample_stats": {
#     "Keyhole": {
#       "Sample": "Keyhole",
#       "Total_reads": "    7146",
#       "Total_UMIs": "    4754",
#       "Duplication_rate": "33.5",
#       "Median_UMIs": null,
#       "Median_Mitochondrial_UMIs_Percent": null,
#       "Cells_100_UMIs": "0",
#       "Cells_FDR_p01": "-"
#     },
# 
def make_sample_stats_dict(sample_name_list, cellread_statistics_dict, umi_cell_statistics_dict, starsolo_summary_dict, hash_read_rate_dict):
  sample_stats_dict = dict()
  for sample_name in sample_name_list:
    total_reads              = cellread_statistics_dict[sample_name]['sum_counted_reads_unique'] + cellread_statistics_dict[sample_name]['sum_counted_reads_multi']
    total_umis               = cellread_statistics_dict[sample_name]['total_umi']
    duplication_rate         = (1.0 - (float(total_umis) / float(total_reads))) * 100.0
    median_umis              = umi_cell_statistics_dict[sample_name]['umis_median']
    median_mitochondial_umis = umi_cell_statistics_dict[sample_name]['umis_mito_median']
    cells_100_umis           = umi_cell_statistics_dict[sample_name]['cell_counts_umi_100_umi_cutoff']
    cells_1000_umis          = umi_cell_statistics_dict[sample_name]['cell_counts_umi_1000_umi_cutoff']
    cells_fdr_p01            = umi_cell_statistics_dict[sample_name]['cell_counts_fdr']
#    cells_100_umis           = umi_cell_statistics_dict[sample_name]['cell_counts_umi']
    if(hash_read_rate_dict.get(sample_name)):
      hash_read_rate         = float(hash_read_rate_dict[sample_name]['hash_reads']) / (float(hash_read_rate_dict[sample_name]['total_reads']) + float(starsolo_summary_dict[sample_name]['number_of_reads']))
      hash_match_rate        = float(hash_read_rate_dict[sample_name]['hash_rate'])
    else:
      hash_read_rate         = 'NA'
      hash_match_rate        = 'NA'

    if(median_umis > 0):
      median_mitochondrial_umis_percent = (float(median_mitochondial_umis) / float(median_umis)) * 100.0
    else:
      median_mitochondrial_umis_percent = 0.0

    stats_dict = dict()
    stats_dict['Sample']                            = sample_name
    stats_dict['Total_reads']                       = '%8d' % total_reads
    stats_dict['Total_UMIs']                        = '%8d' % total_umis
    stats_dict['Duplication_rate']                  = '%.1f' % duplication_rate
    stats_dict['Median_UMIs']                       = '%8d' % median_umis
    stats_dict['Median_Mitochondrial_UMIs_Percent'] = '%.1f' % median_mitochondrial_umis_percent
    stats_dict['Cells_100_UMIs']                    = '%d' % cells_100_umis
    stats_dict['Cells_1000_UMIs']                   = '%d' % cells_1000_umis
    stats_dict['Cells_FDR_p01']                     = '%d' % cells_fdr_p01
    stats_dict['Hash_Read_Rate']                    = '%.3f' % float(hash_read_rate) if(hash_read_rate != 'NA') else 'NA'
    stats_dict['Hash_Match_Rate']                   = '%.3f' % float(hash_match_rate) if(hash_match_rate != 'NA') else 'NA'

    sample_stats_dict[sample_name] = stats_dict

  return(sample_stats_dict)


def make_run_data_dict(processing_directory, sample_name_list, sample_stats_dict, barnyard_sample_name):
  #
  # Calculate cell_counts.
  #
  cell_counts = [0] * 2
  for sample_name in sample_name_list:
    cells_100_umis  = int(sample_stats_dict[sample_name]['Cells_100_UMIs'])
    cells_fdr_p01   = int(sample_stats_dict[sample_name]['Cells_FDR_p01'])
    cell_counts[0] += cells_100_umis if(not math.isnan(cells_100_umis)) else 0
    cell_counts[1] += cells_fdr_p01  if((not math.isnan(cells_fdr_p01)) and (cells_fdr_p01 >= 0))  else 0

  cell_counts_str = [''] * 2
  cell_counts_str[0] = '%d' % (cell_counts[0])
  cell_counts_str[1] = '%d' % (cell_counts[1])

  #
  # Read barnyard collision rate.
  #
  barn_collision = None
  if(barnyard_sample_name != None):
    filename = '%s_barnyard_collision.txt' % (barnyard_sample_name)
    with open(filename, 'r') as fp:
      text_line = fp.readline()
      if(len(text_line) > 0):
        barn_collision = text_line

  #
  # Set up run_data_dict, which becomes the data.js file.
  #
  run_data_dict = dict()
  run_data_dict['run_name'] = processing_directory
  run_data_dict['cell_counts'] = cell_counts_str
  run_data_dict['sample_list'] = sample_name_list
  run_data_dict['barn_collision'] = barn_collision
  run_data_dict['barnyard_sample_name'] = barnyard_sample_name

  #
  # Edit stats_dict['Cells_FDR_p01'] there's no emptyDrops data:
  #   o  umi_cell_statistics_dict[sample_name]['cell_counts_fdr'] is
  #      set to -1 for missing emptyDrops data. Set stats_dict['Cells_FDR_p01']
  #      to '-' in these cases.
  #
  for sample_name in sample_name_list:
    if(int(sample_stats_dict[sample_name]['Cells_FDR_p01']) < 0):
      sample_stats_dict[sample_name]['Cells_FDR_p01'] = '-'

  run_data_dict['sample_stats'] = sample_stats_dict

  return(run_data_dict)


def write_data_js(run_data_dict):
  try:
    filename = 'data.js'
    fp = open('data.js', 'w+')
    fp.write('const run_data =\n')
    fp.write(json.dumps(run_data_dict, indent=2))
  except:
    print('Error: unable to write \"data.js\" file.', file=sys.stderr)
    sys.exit(-1)
  return(0)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program .')
  parser.add_argument('-s', '--sample_map_file', required=True, default=None, help='Input sample map file (required string(s)).')
  parser.add_argument('-p', '--processing_directory', required=True, default=None, help='Input processing directory name (required string(s)).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  # sample_name_list = read_sample_name_file(args.sample_name_file)
  sample_map  = read_json(args.sample_map_file) 

  sample_name_list = get_sample_name_list(sample_map)

  barnyard_sample_name = get_barnyard_sample_name(sample_map)

  #
  # Read sample cellreads statistics.
  #
  cellread_statistics_dict = read_cellread_statistics(sample_name_list)
  # print(json.dumps(cellread_statistics, indent=2))

  #
  # Read STARsolo summary reports.
  #
  starsolo_summary_dict = read_starsolo_summary(sample_name_list)

  #
  # Read sample UMI and cell statistics.
  #
  umi_cell_statistics_dict = read_umi_cell_statistics(sample_name_list)
  # print(json.dumps(umi_cell_statistics, indent=2))

  #
  # Read hash read rates.
  #
  hash_read_rate_dict = read_hash_read_rates(sample_map)

  #
  # Make sample data.
  #
  sample_stats_dict = make_sample_stats_dict(sample_name_list, cellread_statistics_dict, umi_cell_statistics_dict, starsolo_summary_dict, hash_read_rate_dict)

  #
  # Make run data.
  #
  run_data_dict = make_run_data_dict(args.processing_directory, sample_name_list, sample_stats_dict, barnyard_sample_name)

  #
  # Write data.js file.
  #
  write_data_js(run_data_dict)


