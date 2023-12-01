/*
========================================================================================
   Annotation Nextflow Workflow
========================================================================================
   Github   :
   Contact  :
----------------------------------------------------------------------------------------
*/


nextflow.enable.dsl=2

// Display pipeline details
println """\
        T R A N S C R I P T - A N N O T A T I O N - N F   P I P E L I N E
        ===================================
        genome      : ${params.genome}
        fastq       : ${params.reads}
        outdir      : ${params.outdir}
        """
        .stripIndent()
        
        
/*
========================================================================================
   Pipeline Modules
========================================================================================
*/

include { read_conf; get_path; get_genome_desc; create_channel_from_path; UNCOMPRESS } from './modules/common.nf'
include { EOULSAN_READ_FILTER_SR } from './modules/filterreads.nf'
include { EOULSAN_INDEX } from './modules/mapping.nf'
include { EOULSAN_MAPPING } from './modules/mapping.nf'
include { EOULSAN_SAM_FILTER } from './modules/filtersam.nf'
include { EOULSAN_EXPRESSION } from './modules/expression.nf'
include { ISOQUANT } from './modules/isoquant.nf'
include { MERGE_FASTQ } from './modules/merge_fastq.nf'
include { RNA_BLOOM } from './modules/rnabloom.nf'
include { RNABLOOM_MINIMAP2 } from './modules/rnabloom_minimap2.nf'
include { RNABLOOM_PATHTOOLS } from './modules/rnabloom_pathtools.nf'
include { RNABLOOM_AGAT } from './modules/agat.nf'
include { SAMTOOLS } from './modules/samtools.nf'
include { SAMTOOLS_MERGE } from './modules/samtools_merge.nf'

/*
========================================================================================
   Create Channels
========================================================================================
*/
genome_ch = Channel.of( params.genome )
annot_ch = Channel.of( params.annotation )
reads_ch = Channel.fromPath( params.reads, checkIfExists:true )
shortread_ch = params.optional_shortread != null ? file(params.optional_shortread, type: "file") : file("no_shortread", type: "file")


// Pipeline Input parameters
params.mapperName = "minimap2"
params.mapperVersion = "2.24"
params.mapperFlavor = ""
params.indexerArguments = "-x splice"
params.mappersArguments = "-x splice --eqx --secondary=no"
params.tmpDir = projectDir + "/tmp"
params.binaryDir = "/tmp"
params.storages = read_conf()
params.readFilteringConf = [ "trimpolynend" : "" ]
params.samFilteringConf = [ "removeunmapped" : "true", "quality.threshold" : "1", "removesupplementary": "true", "removemultimatches" : "true" ]
params.expressionConf = [ "genomic.type" : "exon", "attribute.id" : "gene_id", "stranded" : "no", "overlap.mode" : "union", "remove.ambiguous.cases" : "false" ]


/*
========================================================================================
   WORKFLOW - Transcript Annotation
========================================================================================
*/

workflow {
  // Index creation
  index_ch = EOULSAN_INDEX(genome_ch, params.mapperName, params.mapperVersion, params.mapperFlavor, params.storages, params.tmpDir, params.binaryDir, params.indexerArguments)
  genome_ch = create_channel_from_path(params.genome, params.storages)
  uncompress_ch = UNCOMPRESS(genome_ch)
  
  // Reads filtering
  filterreads_ch = EOULSAN_READ_FILTER_SR(reads_ch, params.readFilteringConf)

  // Mapping
  filterreads_ch.combine(index_ch).set { reads_index_combined_ch }
  mapping_ch = EOULSAN_MAPPING(reads_index_combined_ch, params.mapperName, params.mapperVersion, params.mapperFlavor, params.tmpDir, params.binaryDir, params.mappersArguments)

  // Alignments filtering
  filtersam_ch = EOULSAN_SAM_FILTER(mapping_ch, params.samFilteringConf, params.tmpDir)

  // Expression computation
  filtersam_ch.combine(annot_ch).combine(genome_ch).set { filtersam_annot_combined_ch }
  expression_ch = EOULSAN_EXPRESSION(filtersam_annot_combined_ch, params.expressionConf, "True", params.storages)
  
  // Launch transcript annotation modules
  SAMTOOLS(EOULSAN_SAM_FILTER.out.filtered_sam)
  SAMTOOLS_MERGE(SAMTOOLS.out.samtools_bam.collect())
  ISOQUANT(uncompress_ch, SAMTOOLS_MERGE.out.samtools_mergedbam)
  
  MERGE_FASTQ(EOULSAN_READ_FILTER_SR.out.eoulsan_fasta.collect())
  RNA_BLOOM(MERGE_FASTQ.out.merged_fastq, shortread_ch)
  RNABLOOM_MINIMAP2(uncompress_ch, RNA_BLOOM.out.rnabloom_fasta)
  RNABLOOM_PATHTOOLS(RNABLOOM_MINIMAP2.out.rnabloom_sam)
  RNABLOOM_AGAT(RNABLOOM_PATHTOOLS.out.rnabloom_bed)
}

// Display pipeline execution summary upon completion
workflow.onComplete {
   println (workflow.success ? """
       Pipeline execution summary
       ---------------------------
       Completed at: ${workflow.complete}
       Duration    : ${workflow.duration}
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/
