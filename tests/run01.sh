#!/bin/bash

samplesheet_csv="samplesheet${1}.csv"
$HOME/git/bbi-scirna-demux/samplesheet/scirna_samplesheet.py -i $samplesheet_csv -o samplesheet${1}.json -n 4
