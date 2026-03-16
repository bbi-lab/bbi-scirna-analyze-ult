
process make_experiment_dashboard {
  errorStrategy 'retry'
  maxRetries 2

  publishDir path: "${params.output_dir}", pattern: "exp_dash", mode: 'copy'

  input:
    path("*") // merge_starsolo_reports.out.cell_reads_stats
    path("*") // merge_starsolo_reports.out.starsolo_summary
    path("*") // make_umi_counts.out.umi_counts_tsv
    path("*") // run_empty_drops.out.empty_drops_fdr
    path("*") // cat_hashes.out.hash_read_rate
    path("*") // make_experiment_dashboard_png_channel_in
    path("*") // make_experiment_dashboard_txt_channel_in
    path(sample_maps_json)
    val(umi_cutoff)
    val(fdr_cutoff)

  output:
    path('exp_dash')


  script:
  """
  # bash watch for errors
  set -ueo pipefail

  #
  # Get sample names from sample_map.json file.
  #
  cat sample_map.json | grep 'sample_name' | sed 's/"sample_name"://' | sed 's/[",]//g' | sed 's/\s*//' > sample_names.txt

  #
  # Extract statistics into sample JSON files.
  #
  for sample_name in `cat sample_names.txt`
  do
    cellread_statistics_tsv="\${sample_name}.starsolo.cell_reads.stats"
    output_cellreads_statistics_json="\${sample_name}_cellreads_statistics.json"
    make_cellread_statistics.py -s \${sample_name} -c ${umi_cutoff} -i \${cellread_statistics_tsv} -o \${output_cellreads_statistics_json}

    umi_counts_tsv="\${sample_name}_umi_counts.tsv"
    empty_drops_fdr_tsv="\${sample_name}_empty_drops_fdr.tsv"
    output_umi_cell_json="\${sample_name}_umi_cell.json"
    make_umi_cell_counts_statistics.py -s \${sample_name} -c ${umi_cutoff} -f ${fdr_cutoff} -u \${umi_counts_tsv} -e \${empty_drops_fdr_tsv} -o \${output_umi_cell_json}
  done

  #
  # Make data.js file for dashboard.
  #
  project_directory=`basename ${workflow.launchDir}`
  make_exp_dash_data_js.py -s ${sample_maps_json} -p \${project_directory}

  #
  # Copy skeleton dashboard to exp_dash directory.
  #
  rm -rf ./exp_dash
  cp -pr ${params.bin_dir}/skeleton_dash ./exp_dash

  #
  # Copy data.js and *.png to exp_dash.
  #
  cp data.js ./exp_dash/js
  cp *.png ./exp_dash/img
  """
}

