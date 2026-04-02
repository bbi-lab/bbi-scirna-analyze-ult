def align_bam_function(item) {
  def sample_name = item['sample_name']
  def in_file = item['in_file']
  def file_path = params.object_map.merge_bam_map[in_file]
  def out_dir = in_file.take(in_file.lastIndexOf('.'))

  return([sample_name, file_path, out_dir])
}


align_cpus = params.align_cpus < 8 ? params.align_cpus : 8
def analyze_out = params.output_dir + '/analyze_out'

process align_bams {
  errorStrategy 'retry'
  maxRetries 2

  clusterOptions { '-l m_mem_free=' + Math.round(sample_map['star_memory'].toFloat() / align_cpus.toFloat()) + 'G -pe serial ' + align_cpus + " -l cpuid_level=${params.cpuid_level}" }

  publishDir path: "${analyze_out}/${sample_name}", pattern: "CellReads.stats", mode: 'copy'

  input:
  tuple val(sample_name), path(bam_in), val(out_dir), val(sample_map)

  output:
  path(out_dir)
  //path(out_dir/Solo.out/GeneFull_Ex50pAS/CellReads.stats)

/*
**  notes:
**    o  need genome information
**    o  need memory information
**    o  need input bam files
**    o  need output bam files
*/

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  STAR_ALIGNER=${task.ext.star_path}

  #
  # Notes:
  #   o  the cell filter does not work for small numbers of
  #      cells and I wonder about using the multiple testing
  #      correction on subsets of the alignments.
  #   o  so I set --soloCellFilter to None.
  #   o  we run emptyDrops on the concatenated raw matrices
  #      later.
  #
  \${STAR_ALIGNER} \
      --runThreadN ${align_cpus} \
      --genomeDir ${sample_map['star_index']} \
      --soloType CB_UMI_Simple \
      --soloBarcodeReadLength 0 \
      --soloCBwhitelist None \
      --soloCBtype String \
      --soloInputSAMattrBarcodeSeq CB UB \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI nM AS GX GN sM \
      --outSJtype None \
      --outFilterMultimapNmax 6 \
      --soloUMIdedup 1MM_All \
      --soloCellReadStats Standard \
      --soloStrand Forward \
      --soloFeatures GeneFull_Ex50pAS \
      --soloMultiMappers PropUnique \
      --soloCellFilter None \
      --readFilesType SAM SE \
      --readFilesIn \
        ${bam_in} \
      --readFilesCommand samtools view \
      --outFileNamePrefix "${out_dir}/"
  """
}

/*
**
** soloCBtype = Sequence
** This was the standard until today (20260217) when I
** change soloCBtype to String.
**
 \${STAR_ALIGNER} \
      --runThreadN ${align_cpus} \
      --genomeDir ${sample_map['star_index']} \
      --soloCBmatchWLtype Exact \
      --soloType CB_UMI_Simple \
      --soloBarcodeMate 0 \
      --soloCBstart 1 \
      --soloCBlen   28 \
      --soloUMIstart 29 \
      --soloUMIlen 8 \
      --soloCBwhitelist None \
      --soloInputSAMattrBarcodeSeq CB UB \
      --soloInputSAMattrBarcodeQual CY UY \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI nM AS GX GN sM \
      --outSJtype None \
      --outFilterMultimapNmax 6 \
      --soloUMIdedup 1MM_All \
      --soloCellReadStats Standard \
      --soloStrand Forward \
      --soloFeatures GeneFull_Ex50pAS \
      --soloMultiMappers PropUnique \
      --soloCellFilter None \
      --readFilesType SAM SE \
      --readFilesIn \
        ${bam_in} \
      --readFilesCommand samtools view \
      --outFileNamePrefix "${out_dir}/"
*/


/*
** Ss barcode+umi tag
**   emptyDrops_CR parameters: nExpectedCells maxPercentile maxMinRatio indMin indMax umiMin umiMinFracMedian candMaxN FDR simN
  \${STAR_ALIGNER} \
      --runThreadN ${align_cpus} \
      --genomeDir ${genome_dir} \
      --soloCBmatchWLtype Exact \
      --soloType CB_UMI_Simple \
      --soloBarcodeMate 0 \
      --soloCBstart 1 \
      --soloCBlen   28 \
      --soloUMIstart 29 \
      --soloUMIlen 8 \
      --soloCBwhitelist None \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
      --outSJtype None \
      --outSAMmultNmax 1 \
      --outSAMstrandField intronMotif \
      --soloUMIdedup Exact \
      --soloCellReadStats Standard \
      --soloStrand Forward \
      --soloFeatures GeneFull_Ex50pAS \
      --soloMultiMappers PropUnique \
      --soloInputSAMattrBarcodeSeq sS \
      --soloInputSAMattrBarcodeQual sQ \
      --readFilesType SAM SE \
      --readFilesIn \
        ${bam_in} \
      --readFilesCommand samtools view \
      --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 60 0.01 20000 0.01 10000 \
      --outFileNamePrefix "${out_dir}/"
*/
