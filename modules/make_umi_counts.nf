def make_umi_counts_function(item) {
  def sample_name = item['sample_name']
  def in_matrix = params.object_map.cat_matrices_raw_map[item['in_matrix']]
  def in_features = params.object_map.cat_matrices_raw_map[item['in_features']]
  def in_barcodes = params.object_map.cat_matrices_raw_map[item['in_barcodes']]
  def out_file = item['out_file']
  return([sample_name, out_file, in_matrix, in_features, in_barcodes])
}


def analyze_out = params.output_dir + '/analyze_out'

process make_umi_counts {
  errorStrategy 'retry'

  publishDir path: "${analyze_out}/${sample_name}", pattern: "*_umi_counts.tsv", mode: 'copy'

  input:
  tuple val(sample_name), val(out_file), path(in_matrix), path(in_features), path(in_barcodes), val(sample_map)

  output:
  tuple val(sample_name), path("*_umi_counts.tsv"), emit: umi_counts_tsv


  script:
  """
  # bash watch for errors
  set -ueo pipefail

  mito_umis -m ${in_matrix} -f ${in_features} -b ${in_barcodes} -a ${sample_map['genes_bed']} -o ${out_file}
  """
}
