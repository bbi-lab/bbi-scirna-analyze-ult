#!/usr/bin/env python3

import sys
import os
import json
import argparse
import re

#
# Program version string.
#
program_version = '0.1.0'

# struct SampleMap {
#   sample_id: String,
#   ranges: String,
#   lanes: String,
#   tissue: String,
#   genome: String,
#   hash_file: String,
#   sample_flags: String,
#   external_sample_name: String,
#   wrap_group: String,
#   rt_file: String,
#   ligation_file: String,
#   p7_file: String,
#   p5_file: String,
#   library: String,
#   process_group: String
# }


#
# Read samplesheet json file.
#
def read_json(filename):
  with open(filename, 'r') as fh:
    json_data = json.load(fh)
  return(json_data)


def read_star_genomes_file(file_path):
  star_genomes_dict = {}
  with open(file_path, 'r') as fp:
    for line in fp:
      line = line.strip()
      parts = line.split()
      genome_name = parts[0]
      genome_path = parts[1]
      genome_mem  = parts[2]
      if(not genome_name in star_genomes_dict):
        star_genomes_dict[genome_name] = {}
      star_genomes_dict[genome_name]['genome_path'] = genome_path
      star_genomes_dict[genome_name]['genome_mem']   = genome_mem
  return(star_genomes_dict)


def make_genome_files_dict(json_data, genome_files_dict):
  genome_files_dict = dict()
  for sample_index_dict in json_data['sample_index_list']:
    sample_name = sample_index_dict['sample_id']
    process_group = sample_index_dict['process_group']
    key = '%s-%03d' % (sample_name, int(process_group))
    if(not key in genome_files_dict):
      genome_name = sample_index_dict['genome']
      genome_path = star_genomes_dict[genome_name]['genome_path']
      genome_root = os.path.dirname(genome_path)
      genome_species = os.path.basename(genome_root)
      tmp_genes_bed = os.path.join(genome_root, '%s_rna' % genome_species, 'tmp.genes.bed')
      latest_genes_bed = os.path.join(genome_root, '%s_rna' % genome_species, 'latest.genes.bed')
      hash_file = sample_index_dict['hash_file']

      tmp_dict = dict()
      tmp_dict['sample_name'] = key
      tmp_dict['genome'] = genome_name
      tmp_dict['tmp_genes_bed'] = tmp_genes_bed
      tmp_dict['latest_genes_bed'] = latest_genes_bed
      tmp_dict['hash_file'] = hash_file
      genome_files_dict[key] = tmp_dict
  return(genome_files_dict)


#
# Make JSON file for sourcing genome files.
#
# [
#   {
#     "sample_name": "SeahubZ01-001",
#     "genome": "Zebrafish",
#     "tmp_genes_bed": "/net/bbi/vol1/data/genomes_stage/zebrafish/zebrafish_rna/tmp.genes.bed",
#     "latest_genes_bed": "/net/bbi/vol1/data/genomes_stage/zebrafish/zebrafish_rna/latest.genes.bed"
#   },
#   {
#     "sample_name": "Keyhole-001",
#     "genome": "Fishbowl_seahub",
#     "tmp_genes_bed": "/net/bbi/vol1/data/genomes_stage/fishbowl_seahub/fishbowl_seahub_rna/tmp.genes.bed",
#     "latest_genes_bed": "/net/bbi/vol1/data/genomes_stage/fishbowl_seahub/fishbowl_seahub_rna/latest.genes.bed"
#   },
# ...
# ]
#
def make_genome_files_json(genome_files_dict):
  genome_files_list = []
  for sample_name in genome_files_dict:
    genome_files_list.append(genome_files_dict[sample_name])

  print(genome_files_list)

  try:
    filename_json = 'genome_files.json'
    fh = open(filename_json, 'w+')
    json.dump(genome_files_list, fh, indent=2)
  except:
    print('Error: unable to write output file \"%s\"' % (filename_json), file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to make JSON file for sourcing genome files.')
  parser.add_argument('-i', '--input', required=True, default=None, help='Input JSON samplesheet filename (required string).')
  parser.add_argument('-g', '--genomes_file', required=True, default=None, help='Path to file of STAR aligner genomes data. (required string).')
#  parser.add_argument('-o', '--output', required=False, default=None, help='Output JSON filename (required string).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  #
  # Read input samplesheet file.
  #
  json_data = read_json(args.input)
  star_genomes_dict = read_star_genomes_file(args.genomes_file)
#  print(json.dumps(json_data['sample_index_list'], indent=2))
  genome_files_dict = make_genome_files_dict(json_data, star_genomes_dict)
#  print(genome_files_dict)
  make_genome_files_json(genome_files_dict)

