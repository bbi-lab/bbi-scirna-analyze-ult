def merge_starsolo_reports_function(item) {
  def sample_name = item['sample_name']
  def out_file = sample_name + '.starsolo.cell_reads.stats'
  def in_dir_list = []
  def in_root_list = []
  for(def in_dir in item['in_dir_list']) {
    def dir_base_name = in_dir.toString().tokenize('/').last()
    def root_path = params.object_map.merge_align_bam_map[dir_base_name]
    def file_path = params.object_map.merge_align_bam_map[dir_base_name] + '/Solo.out/GeneFull_Ex50pAS/CellReads.stats'
    /*
    ** If the file/value does not exist in
    ** params.object_map.merge_bam_map,
    ** skip this pipeline entry.
    */
    if(file_path == null) {
      continue
    }
    in_dir_list.add(file_path)
    in_root_list.add(root_path)
  }
  return([sample_name, out_file, in_dir_list, in_root_list])
}


def analyze_out = params.output_dir + '/analyze_out' 

process merge_starsolo_reports {
  errorStrategy 'retry'
  maxRetries 2

  publishDir path: "${analyze_out}/${sample_name}", pattern: "*.starsolo.cell_reads.stats", mode: 'copy'
  publishDir path: "${analyze_out}/${sample_name}", pattern: "*Features.stats", mode: 'copy'
  publishDir path: "${analyze_out}/${sample_name}", pattern: "*Summary.txt", mode: 'copy'

  input:
  tuple val('sample_name'), val('out_file'), path('file'), path('root')

  output:
  tuple val(sample_name), path("*.starsolo.cell_reads.stats"), emit: cell_reads_stats
  tuple val(sample_name), path("*Summary.txt"), emit: starsolo_summary
//  tuple path("*Features.stats"), path("*Summary.txt"), emit: dummy

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  cat_starsolo_stats.py -i ${file} -o ${out_file}

  merge_starsolo_reports.py -i ${root} -s ${sample_name}
  """
}

