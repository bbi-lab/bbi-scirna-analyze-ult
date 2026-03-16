process merge_hash_reads {
  cache 'lenient'
  errorStrategy 'retry'
  maxRetries 2

  input:
  tuple val('sample_name'), val('out_file'), path('files')

  output:
  path("*.hash_reads.merged.tsv")

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  file_list=`ls files*`
  $workflow.projectDir/bin/merge_hash_reads.py -i \${file_list} -o ${out_file}
  """
}

