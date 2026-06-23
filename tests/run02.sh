#!/bin/bash


samplesheet_json="samplesheet${1}.json"
out_dir='reference_jsons/'

mkdir -p  ${out_dir}

$HOME/git/bbi-scirna-analyze-ult/bin/make_merge_demux_json.py -i ${samplesheet_json} -p '.'
mv merge_demux.json ${out_dir}merge_demux${1}.json

$HOME/git/bbi-scirna-analyze-ult/bin/make_process_hashes_json.py -i ${samplesheet_json}
mv process_hashes.json ${out_dir}process_hashes${1}.json

$HOME/git/bbi-scirna-analyze-ult/bin/make_merge_hash_reads_json.py -i ${samplesheet_json} -p '.'
mv merge_hash_reads.json ${out_dir}merge_hash_reads${1}.json

$HOME/git/bbi-scirna-analyze-ult/bin/make_star_align_json.py -i ${samplesheet_json}
mv star_align.json ${out_dir}star_align${1}.json

$HOME/git/bbi-scirna-analyze-ult/bin/make_merge_align_json.py -i ${samplesheet_json}
mv merge_align.json ${out_dir}merge_align${1}.json

$HOME/git/bbi-scirna-analyze-ult/bin/make_umi_counts_json.py -i ${samplesheet_json}
mv umi_counts.json ${out_dir}umi_counts${1}.json
