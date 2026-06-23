def analyze_out = params.output_dir + '/analyze_out'

process run_scrublet {
  errorStrategy 'retry'
  maxRetries 2

  publishDir path: "${analyze_out}/${sample_name}", pattern: "*_scrublet_out.csv", mode: 'copy'
  publishDir path: "${analyze_out}/${sample_name}", pattern: "run_scrublet.log", mode: 'copy'
  publishDir path: "${analyze_out}/${sample_name}", pattern: "${sample_name}_cds.raw.mobs", mode: 'copy'
  publishDir path: "${analyze_out}/${sample_name}", pattern: "${sample_name}_cds.raw.col_data.tsv", mode: 'copy'

  input:
  tuple val(sample_name), path(mobs), path(umi_counts), val(sample_map)

  output:
  tuple val(sample_name), path("*scrublet_out.csv"), path("*.png"), path('run_scrublet.log'), emit: scrublet_out
  tuple val(sample_name), path("${sample_name}_cds.raw.mobs", includeInputs: true), path(umi_counts), emit: cds
  tuple val(sample_name), path("${sample_name}_cds.raw.col_data.tsv"), emit: col_data

  /*
  ** Don't exit on error. Continue so that
  ** the merge_align is not aborted.
  */
  shell '/bin/bash', '-u'

  script:
  """

  #
  # Move the input mobs directory.
  #

  matrix_filename='scrublet_expression_matrix.mtx'
  if [ "$params.run_scrublet" == "true" ]
  then
    mv ${mobs} tmp.in.mobs
    write_expression_matrix.R $sample_name tmp.in.mobs \$matrix_filename
    run_scrublet.py --sample_name $sample_name --mat \$matrix_filename --run_scrublet > run_scrublet.log
    rm \$matrix_filename
    #
    # Note:
    #   add_scrublet_to_cds.R writes a .mobs directory that has the
    #   name ${sample_name}_cds.raw.mobs. This mobs cds has the
    #   scrublet scores.
    #
    add_scrublet_to_cds.R ${sample_name} tmp.in.mobs ${sample_name}_scrublet_out.csv
  else
    run_scrublet.py --sample_name $sample_name --mat \$matrix_filename > run_scrublet.log 
  fi

  write_col_data.R ${sample_name}_cds.raw.mobs ${sample_name}_cds.raw.col_data.tsv
  """
}

