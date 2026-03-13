process merge_demux {
  cache 'lenient'
  errorStrategy 'retry'
  maxRetries 2

  input:
  tuple val('sample_name'), val('out_file'), path('files')

  output:
  path("*.merged.bam")

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  file_list=`ls files*`
  for file in \$file_list
  do
    samtools sort -@ 4 -m 64G \${file} -o \${file}.sorted
  done
  samtools merge -@ 4 ${out_file} *.sorted
  rm -r *.sorted
  """
}

