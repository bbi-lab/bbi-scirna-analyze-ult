def analyze_out = params.output_dir + '/analyze_out' 


/*
** STARsolo 'raw' matrix concatenation.
*/
def cat_matrices_raw_function(item) {
  def sample_name = item['sample_name']
  def out_file = sample_name + '_counts.raw'
  def in_file_list = []
  for(def in_dir in item['in_dir_list']) {
    def dir_base_name = in_dir.toString().tokenize('/').last()
    def file_path = params.object_map.merge_align_bam_map[dir_base_name] + '/Solo.out/GeneFull_Ex50pAS/raw/UniqueAndMult-PropUnique.mtx'
//    def file_path = params.object_map.merge_align_bam_map[dir_base_name] + '/Solo.out/GeneFull_Ex50pAS/raw/matrix.mtx'
    /*
    ** If the file/value does not exist in
    ** params.object_map.merge_bam_map,
    ** skip this pipeline entry.
    */
    if(file_path == null) {
      continue
    }
    in_file_list.add(file_path)
  }
  return([sample_name, out_file, in_file_list])
}


process cat_matrices_raw {
  errorStrategy 'retry'
  maxRetries 2

  publishDir path: "${analyze_out}/${sample_name}", pattern: "*.cells.tsv", mode: 'copy'
  publishDir path: "${analyze_out}/${sample_name}", pattern: "*.features.tsv", mode: 'copy'
  publishDir path: "${analyze_out}/${sample_name}", pattern: "*.matrix.mtx", mode: 'copy'
  publishDir path: "${analyze_out}/${sample_name}", pattern: "*.cells.barcode_to_wells.tsv", mode: 'copy'

  input:
  tuple val('sample_name'), val('out_file'), path('file')

  output:
  tuple val(sample_name), path("*.cells.tsv"), path("*.features.tsv"), path("*.matrix.mtx"), path("*.cells.barcode_to_wells.tsv"), emit: raw_matrix

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  # Use the symbolic link referent because we use the path
  # to find the feature and cell files for cat_sparse_matrix.py.
  #
  file_list=`ls file*`

  in_file_list=''
  for file in \${file_list}
  do
    file_path=`readlink \$file`
    in_file_list="\${in_file_list} \${file_path}"
  done

  if [ -n "\$in_file_list" ]
  then
    cat_sparse_matrix.py -i \$in_file_list -m 'UniqueAndMult-PropUnique.mtx' -f 'features.tsv' -c 'barcodes.tsv' -o ${out_file}
#    cat_sparse_matrix.py -i \$in_file_list -m 'matrix.mtx' -f 'features.tsv' -c 'barcodes.tsv' -o ${out_file}
  fi

  barcode_to_well.py -i ${out_file}.cells.tsv -o ${out_file}.cells.barcode_to_wells.tsv
  """
}


/*
** STARsolo 'filtered' matrix concatenation.
*/
def cat_matrices_filtered_function(item) {
  def sample_name = item['sample_name']
  def out_file = sample_name + '_counts.filtered'
  def in_file_list = []
  for(def in_dir in item['in_dir_list']) {
    def dir_base_name = in_dir.toString().tokenize('/').last()
    def file_path = params.object_map.merge_align_bam_map[dir_base_name] + '/Solo.out/GeneFull_Ex50pAS/filtered/matrix.mtx'
    /*
    ** If the file/value does not exist in
    ** params.object_map.merge_bam_map,
    ** skip this pipeline entry.
    */
    if(file_path == null) {
      continue
    }
    in_file_list.add(file_path)
  }
  return([sample_name, out_file, in_file_list])
}


process cat_matrices_filtered {
  errorStrategy 'retry'
  maxRetries 2

  publishDir path: "${analyze_out}/${sample_name}", pattern: "*.cells.tsv", mode: 'copy'
  publishDir path: "${analyze_out}/${sample_name}", pattern: "*.features.tsv", mode: 'copy'
  publishDir path: "${analyze_out}/${sample_name}", pattern: "*.matrix.mtx", mode: 'copy'

  input:
  tuple val('sample_name'), val('out_file'), path('file')

  output:
  tuple val(sample_name), path("*.cells.tsv"), path("*.features.tsv"), path("*.matrix.mtx"), path("*.cells.barcode_to_wells.tsv"), emit: filtered_matrix

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  #
  #
  file_list=`ls file*`

  #
  # Note: when the STAR input BAM file has few or no reads (after trimming)
  #       the filtered directory does not exist.
  #
  in_file_list=''
  for file in \${file_list}
  do
    file_path=`readlink \$file`
    if [ -f \${file_path} ]
    then
      in_file_list="\${in_file_list} \${file_path}"
    fi
  done

  if [ -n "\$in_file_list" ]
  then
    cat_sparse_matrix.py -i \$in_file_list -m 'matrix.mtx' -f 'features.tsv' -c 'barcodes.tsv' -o ${out_file}
  fi

  barcode_to_well.py -i ${out_file}.cells.tsv -o ${out_file}.cells.barcode_to_wells.tsv
  """
}

