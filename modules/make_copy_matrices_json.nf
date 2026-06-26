def copy_matrices_function(item) {
  def in_file_list = []
  for(def in_dir in item['in_dir_list']) {
    def dir_base_name = in_dir.toString().tokenize('/').last()
    def file_path = params.object_map.merge_align_bam_map[dir_base_name] + '/Solo.out/GeneFull_Ex50pAS/raw/matrix.mtx'
    in_file_list.add(file_path)
  }
  return(in_file_list)
}


process make_copy_matrices_json {

  input:
  path(samplesheet_file)
  val(dummy)

  output:
  path("copy_matrices.json")

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  $workflow.projectDir/bin/make_copy_matrices_json.py -i $samplesheet_file
  """
}
