def merge_align_function(item) {
  def sample_name = item['sample_name']
  def out_file = item['out_file']
  def in_dir_list = []
  for(def in_dir in item['in_dir_list']) {
    def dir_base_name = in_dir.toString().tokenize('/').last()
    def file_path = params.object_map.merge_align_bam_map[dir_base_name] + '/Aligned.sortedByCoord.out.bam'
    /*
    ** If the file/value does not exist in
    ** params.object_map.merge_bam_map,
    ** skip this pipeline entry.
    */
    if(file_path == null) {
      continue
    }
    in_dir_list.add(file_path)
  }
  return([sample_name, out_file, in_dir_list])
}


def analyze_out = params.output_dir + '/analyze_out' 

process merge_align {
  errorStrategy 'retry'
  maxRetries 2

  publishDir path: "${analyze_out}/${sample_name}", pattern: "*aligned.bam", mode: 'copy'

  input:
  tuple val('sample_name'), val('out_file'), path('files')

  output:
  path("*aligned.bam")

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  nfil=`ls files* | wc -l`

  if [ "\${nfil}" -gt 1 ]
  then
    sambamba merge -t 8 ${out_file} files*
  else
    cp files* ${out_file}
  fi
  """
}

