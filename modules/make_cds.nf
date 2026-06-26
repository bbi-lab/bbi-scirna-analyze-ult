def analyze_out = params.output_dir + '/analyze_out'

process make_cds_raw {
  errorStrategy 'ignore'

  publishDir path: "${analyze_out}/${sample_name}", pattern: "*.png", mode: 'copy'

  input:
  tuple val(sample_name), path(cell_tsv), path(feature_tsv), path(count_matrix), path(barcode_to_wells), path(counts_per_cell), path(empty_drops), path(umi_counts), val(sample_map)
  val(out_file)
  val(umi_cutoff)

  output:
  tuple val(sample_name), path("*.raw.mobs"), path(umi_counts),  emit: cds
  tuple val(sample_name), path("*.png"), emit: png

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  make_cds.R \
  ${sample_name} \
  'raw' \
  ${count_matrix} \
  ${feature_tsv} \
  ${cell_tsv} \
  ${barcode_to_wells} \
  ${umi_counts} \
  ${umi_cutoff} \
  ${counts_per_cell} \
  ${sample_map['genes_bed']} \
  ${empty_drops}
  """
}


/*
** This is unmaintained: there may be errors changes to make_cds_raw that
** were not transferred to make_cds_filtered.
**
process make_cds_filtered {
  errorStrategy 'ignore'

  publishDir path: "${analyze_out}/${sample_name}", pattern: "*.filtered.mobs", mode: 'copy'
  publishDir path: "${analyze_out}/${sample_name}", pattern: "*.png", mode: 'copy'

  input:
//  tuple val(sample_name), path(cell_tsv), path(feature_tsv), path(count_matrix), path(barcode_to_wells)
  tuple val(sample_name), path(cell_tsv), path(feature_tsv), path(count_matrix), path(barcode_to_wells), path(counts_per_cell), path(umi_counts), val(genome), path(latest_genes_bed), val(hash_file), val(sample_map)
  val(out_file)
  val(umi_cutoff)

  output:
  tuple val(sample_name), val(genome), path("*.filtered.mobs"), path(umi_counts), val(hash_file), emit: cds
  tuple val(sample_name), path("*.png"), emit: png

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  make_cds.R \
  ${sample_name} \
  'filtered' \
  ${count_matrix} \
  ${feature_tsv} \
  ${cell_tsv} \
  ${barcode_to_wells} \
  ${umi_counts} \
  ${umi_cutoff} \
  ${counts_per_cell} \
  ${sample_map['genes_bed']} \
  ${empty_drops}
  """
}
*/

