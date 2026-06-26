#!/usr/bin/env python3

import sys
import os
import argparse
import csv
import re
import math


#
# Program version string.
#
program_version = '0.1.0'

#
# The STARsolo manual overlooks some documentation details. I
# found the following information online.
#
# Features.stats

# URL: https://github.com/alexdobin/STAR/issues/1887
# Marta Benegas:
# 
# Ok, let's see if I got it:
# 
# Barcoded Reads Stats:
# Reads not used for quantification because:
# 
# noNoAdapter: Reads without an Adapter.
# noNoUMI: Reads not Associated to a Valid UMI
# noNoCB: Reads not Associated to a Valid Cell Barcode
# noNinCB: Reads with N's in Cell Barcodes
# noNinUMI: Reads with N's in UMIs
# noUMIhomopolymer: Reads associated with a Homopolymeric UMI
# noNoWLmatch: Reads Without Match to the Whitelist
# noTooManyMM: Reads With Too Many Mismatches to the Whitelist
# noTooManyWLmatches: Reads With Too Many matches to the Whitelist
# Reads used for the quantification:
# 
# yesWLmatchExact: Reads with Exact Match to Whitelist
# yesOneWLmatchWithMM: Reads with One Match to the Whitelist with Mismatches
# yesMultWLmatchWithMM: Reads with Multiple Matches to the Whitelist with Mismatches
# Feature ReadsStats
# Reads not used for quantification because:
# 
# noUnmapped: Unmapped Reads
# noNoFeature: Reads not Mapped to a Feature
# MultiFeature: Reads Aligned to Multiple Features
# subMultiFeatureMultiGenomic ??
# noTooManyWLmatches: Reads Not Counted Because their Barcoded pair has Too Many Matches to the Whitelist
# noMMtoWLwithoutExact: Reads Not Counted Because their Barcoded pair has Mismatches to the Whitelist and there's no more reads supporting that barcode
# Reads used for the quantification:
# 
# yesWLmatch: Reads whose Barcoded Pair has a Match to the Whitelist
# yessubWLmatchExact: ??
# yessubWLmatch_UniqueFeature: ??
# yesCellBarcodes: Reads Assossiated to a Valid Cell Barcode
# yesUMIs: Reads Assossiated to a Valid UMI
# 
# 
# 
# alexdobin:
# subMultiFeatureMultiGenomic: reads mapping to multiple genomic loci and multiple features (a subset of MultiFeature)
# yessubWLmatchExact: reads with cell barcode exactly matched to the whitelist (a subset of yesWLmatch)
# yessubWLmatch_UniqueFeature: reads with matched to the WL and unique feature (a subset of yesWLmatch)
#
#
# URL: https://broadinstitute.github.io/warp/docs/Pipelines/Optimus_Pipeline/starsolo-metrics
# 
# Number of Reads:             Number of reads in the library.
# Reads With Valid Barcodes:   Fraction of reads with valid barcodes.
# Sequencing Saturation:       Proportion of unique molecular identifiers (UMIs) 
#                              that have been sequenced at least once compared 
#                              to the total number of possible UMIs in the sample; 
#                              calculated as: 1-(yesUMIs/yessubWLmatch_UniqueFeature).
# Q30 Bases in CB+UMI:         Fraction of high-quality reads in the cell barcode and UMI read.
# Q30 Bases in RNA read:       Fraction of high-quality reads in the RNA read.
# Reads Mapped to Genome: Unique+Multiple:  Fraction of unique and multimapped reads that mapped to the genome.
# Reads Mapped to Genome: Unique:           Fraction of unique reads that mapped to the genome.
# Reads Mapped to genes: Unique+Multiple:   Fraction of reads that mapped to genes as defined by the �~@~Ssolo-feature parameter.
# Reads Mapped to Genes: Unique:            Fraction of unique reads that mapped to genes.
# Estimated Number of Cells:   Number of barcodes that STARsolo flagged as cells based on UMIs.
# Unique Reads in Cells Mapped to genes:    Total number of unique reads that mapped to genes across all cells
# Fraction of Unique Reads in Cells:        Fraction of unique reads across all cells.
# Mean Reads per Cell:         Mean number of reads per cell.
# Median Reads per Cell:       Median number of reads per cell.
# UMIs in Cells:               Number of UMIs per cell.
# Mean UMI per Cell:           Mean number of UMIs per cell.
# Median UMI per Cell:         Median number of UMI per cell.
# Mean Genes per Cell:         Mean number of genes expressed per cell.
# Median Genes per Cell:       Median number of genes per cell.
# Total Genes Detected:        Total number of genes detected in the overall library.
#
#
# See also:
#   STAR/source/SoloFeature_statsOutput.cpp
#   STAR/source/SoloFeature_processRecords.cpp
#


def read_features_file(filename):
  with open(filename, 'r') as ifh:
    features_dict = dict()
    for line in ifh:
      parts = line.split()
      features_dict[parts[0].strip()] = int(parts[1].strip())
  return(features_dict)


# from https://stackoverflow.com/questions/46647744/checking-to-see-if-a-string-is-an-integer-or-float: Rafe
re_int = re.compile(r"(^[1-9]+\d*$|^0$)")
re_float = re.compile(r"(^\d+\.\d+$|^\.\d+$)")

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


# Features.stats
#  noUnmapped                   13029578
#  noNoFeature                  25056517
#  MultiFeature                 1932371
#  subMultiFeatureMultiGenomic  1445332
#  noTooManyWLmatches           0
#  noMMtoWLwithoutExact         0
#  yesWLmatch                   53641621
#  yessubWLmatchExact           53641621
#  yessubWLmatch_UniqueFeature  51709250
#  yesCellBarcodes              282521
#  yesUMIs                      19539197


# Summary.csv
#  Number of Reads,93402063
#  Reads With Valid Barcodes,1
#  Sequencing Saturation,0.622133
#  Q30 Bases in CB+UMI,-nan
#  Q30 Bases in RNA read,0.880727
#  Reads Mapped to Genome: Unique+Multiple,0.857871
#  Reads Mapped to Genome: Unique,0.742362
#  Reads Mapped to GeneFull_Ex50pAS: Unique+Multiple GeneFull_Ex50pAS,0.574309
#  Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS,0.55362
#  Estimated Number of Cells,7229
#  Unique Reads in Cells Mapped to GeneFull_Ex50pAS,43301562
#  Fraction of Unique Reads in Cells,0.837405
#  Mean Reads per Cell,5989
#  Median Reads per Cell,4531
#  UMIs in Cells,16286345
#  Mean UMI per Cell,2252
#  Median UMI per Cell,1817
#  Mean GeneFull_Ex50pAS per Cell,1033
#  Median GeneFull_Ex50pAS per Cell,948
#  Total GeneFull_Ex50pAS Detected,26766


def float_test(obj):
  if(isinstance(obj, str)):
    return(False)
  elif(math.isnan(obj)):
    return(False)
  return(True)


def make_stats_bundle(starsolo_path, features, summary):
  counts_dict = dict()
  num_read = 0
  for key in summary.keys():
    if(key == 'Number of Reads'):
      counts_dict[key] = int(summary[key])
      num_read = int(summary[key])
    elif(key == 'Reads With Valid Barcodes'):
      counts_dict[key] = int(float(summary[key]) * num_read) if float_test(summary[key]) else summary[key]
    elif(key == 'Sequencing Saturation'):
      counts_dict[key] = -1
    elif(key == 'Q30 Bases in CB+UMI'):
      continue
    elif(key == 'Q30 Bases in RNA read'):
      continue
    elif(key == 'Reads Mapped to Genome: Unique+Multiple'):
      counts_dict[key] = int(float(summary[key]) * num_read) if float_test(summary[key]) else summary[key]
    elif(key == 'Reads Mapped to Genome: Unique'):
      counts_dict[key] = int(float(summary[key]) * num_read)  if float_test(summary[key]) else summary[key]
    elif(key == 'Reads Mapped to GeneFull_Ex50pAS: Unique+Multiple GeneFull_Ex50pAS'):
      counts_dict[key] = int(float(summary[key]) * num_read)  if float_test(summary[key]) else summary[key]
    elif(key == 'Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS'):
       counts_dict[key] = int(float(summary[key]) * num_read)  if float_test(summary[key]) else summary[key]
    elif(key == 'Estimated Number of Cells'):
      counts_dict[key] = int(summary[key])
    elif(key == 'Unique Reads in Cells Mapped to GeneFull_Ex50pAS'):
      counts_dict[key] = int(summary[key])
    elif(key == 'Fraction of Unique Reads in Cells'):
      counts_dict[key] = int(float(summary[key]) * counts_dict['Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS']) if float_test(summary[key]) else summary[key]
    elif(key == 'Mean Reads per Cell'):
      continue
    elif(key == 'Median Reads per Cell'):
      continue
    elif(key == 'UMIs in Cells'):
      counts_dict[key] = int(summary[key])
    elif(key == 'Mean UMI per Cell'):
      continue
    elif(key == 'Median UMI per Cell'):
      continue
    elif(key == 'Mean GeneFull_Ex50pAS per Cell'):
      continue
    elif(key == 'Median GeneFull_Ex50pAS per Cell'):
      continue
    elif(key == 'Total GeneFull_Ex50pAS Detected'):
      counts_dict[key] = int(summary[key])

  if(num_read == 0):
    return(None)

  for key in features.keys():
    if(key == 'noUnmapped'):
      counts_dict[key] = int(features	[key])
    elif(key == 'noNoFeature'):
      counts_dict[key] = int(features	[key])
    elif(key == 'MultiFeature'):
      counts_dict[key] = int(features	[key])
    elif(key == 'subMultiFeatureMultiGenomic'):
      counts_dict[key] = int(features	[key])
    elif(key == 'noTooManyWLmatches'):
      counts_dict[key] = int(features	[key])
    elif(key == 'noMMtoWLwithoutExact'):
      counts_dict[key] = int(features	[key])
    elif(key == 'yesWLmatch'):
      counts_dict[key] = int(features	[key])
    elif(key == 'yessubWLmatchExact'):
      counts_dict[key] = int(features	[key])
    elif(key == 'yessubWLmatch_UniqueFeature'):
      counts_dict[key] = int(features	[key])
    elif(key == 'yesCellBarcodes'):
      counts_dict[key] = int(features	[key])
    elif(key == 'yesUMIs'):
      counts_dict[key] = int(features	[key])

  return(counts_dict)


def merge_stats_bundle_dicts(stats_bundle_dict):
  stats_bundles_merged = dict()
  for dataname in stats_bundle_dict.keys():
    stats_bundle = stats_bundle_dict[dataname]
    for key in stats_bundle:
      count = stats_bundles_merged.setdefault(key, 0)
      val = stats_bundle[key] if (not isinstance(stats_bundle[key], str)) and (not math.isnan(stats_bundle[key])) else 0
      # print('data: %s  val: %f' % (dataname, val), file=sys.stderr)
      stats_bundles_merged[key] = count + val

  return(stats_bundles_merged)


def write_report(stats_bundle_merged, samplename):
  #
  # Report Features.stats.
  #
  filename = '%s_Features.stats' % (samplename)
  with open(filename, 'w') as ofh:
    for key in stats_bundle_merged.keys():
      if(key in ['noUnmapped',
                 'noNoFeature',
                 'MultiFeature', 
                 'subMultiFeatureMultiGenomic',
                 'noTooManyWLmatches',
                 'noMMtoWLwithoutExact',
                 'yesWLmatch',
                 'yessubWLmatchExact',
                 'yessubWLmatch_UniqueFeature',
                 'yesCellBarcodes',
                 'yesUMIs']):
        print('%s\t%d' % (key, stats_bundle_merged[key]), file=ofh)

  #
  # Report summary.csv values.
  #
  num_read = 0
  for key in stats_bundle_merged.keys():
    if(key == 'Number of Reads'):
      num_read = stats_bundle_merged[key]

  filename = '%s_Summary.txt' % (samplename)
  with open(filename, 'w') as ofh:
    print(file=ofh)
    sep = '\t'
    key = 'Number of Reads'
    print('%s%s%d' % (key, sep, int(stats_bundle_merged[key])), file=ofh)
    key = 'Reads With Valid Barcodes'
    print('%s%s%.3f' % (key, sep, float(stats_bundle_merged[key]) / float(num_read)), file=ofh)
    key = 'Sequencing Saturation'
    print('%s%s%.3f' % (key, sep, 1.0 - float(stats_bundle_merged['yesUMIs']) / float(stats_bundle_merged['yessubWLmatch_UniqueFeature'])), file=ofh)
    key = 'Reads Mapped to Genome: Unique+Multiple'
    print('%s%s%.3f' % (key, sep, float(stats_bundle_merged[key]) / float(num_read)), file=ofh)
    key = 'Reads Mapped to Genome: Unique'
    print('%s%s%.3f' % (key, sep, float(stats_bundle_merged[key]) / float(num_read)), file=ofh)
    key = 'Reads Mapped to GeneFull_Ex50pAS: Unique+Multiple GeneFull_Ex50pAS'
    print('%s%s%.3f' % (key, sep, float(stats_bundle_merged[key]) / float(num_read)), file=ofh)
    key = 'Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS'
    print('%s%s%.3f' % (key, sep, float(stats_bundle_merged[key]) / float(num_read)), file=ofh)

#
# The following are calculated using the filtered data.
#
#    key = 'Estimated Number of Cells'
#    print('%s%s%d' % (key, sep, int(stats_bundle_merged[key])), file=ofh)
#    key = 'Unique Reads in Cells Mapped to GeneFull_Ex50pAS'
#    print('%s%s%d' % (key, sep, int(stats_bundle_merged[key])), file=ofh)
#    key = 'Fraction of Unique Reads in Cells'
#    print('%s%s%.3f' % (key, sep, float(stats_bundle_merged[key]) / float(stats_bundle_merged['Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS'])), file=ofh)
#    key = 'UMIs in Cells'
#    print('%s%s%d' % (key, sep, int(stats_bundle_merged[key])), file=ofh)
#    key = 'Total GeneFull_Ex50pAS Detected'
#    print('%s%s%d' % (key, sep, int(stats_bundle_merged[key])), file=ofh)

  return(0)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to gather STARsolo Features.stats files.')
  parser.add_argument('-i', '--input', required=True, default=None, nargs='+', help='STARsolo output paths (required string).')
  parser.add_argument('-s', '--samplename', required=True, default=None, help='Sample name.')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  #
  # path: analyze_out/REF3-001/REF3-001_085_049.trimmed/Solo.out/GeneFull_Ex50pAS/(Features.stats|Summary.csv)
  #
  stats_bundle_dict = dict()
  for i in range(len(args.input)):
    starsolo_path = args.input[i]
    features_filename = starsolo_path + '/Solo.out/GeneFull_Ex50pAS/Features.stats'
    features = read_features_file(features_filename)
    summary_filename = starsolo_path + '/Solo.out/GeneFull_Ex50pAS/Summary.csv'
    summary = read_summary_file(summary_filename)
    stats_bundle = make_stats_bundle(starsolo_path, features, summary)
    if(stats_bundle is None):
      continue
    dataname = os.path.basename(args.input[i])
    if(dataname in stats_bundle_dict):
      print('Error: basename is not distinct: \'%s\'' % (dataname), file=sys.stderr)
      sys.exit(1)
    stats_bundle_dict[dataname] = stats_bundle
    
  stats_bundles_merged = merge_stats_bundle_dicts(stats_bundle_dict)

  write_report(stats_bundles_merged, args.samplename)


