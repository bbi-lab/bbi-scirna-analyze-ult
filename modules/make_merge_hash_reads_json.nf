def json_files_out = params.output_dir + '/json_files'

process make_merge_hash_reads_json {
  publishDir path: "${json_files_out}", pattern: "merge_hash_reads.json", mode: 'copy'

  input:
  path(samplesheet_file)
  val(tsv_path)

  output:
  path("merge_hash_reads.json")

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  $workflow.projectDir/bin/make_merge_hash_reads_json.py -i $samplesheet_file -p $tsv_path
  """
}
