#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { read_conf; get_path; get_genome_desc } from './modules/common.nf'
include { EOULSAN_READ_FILTER_SR } from './modules/filterreads.nf'
include { EOULSAN_INDEX } from './modules/mapping.nf'
include { EOULSAN_MAPPING } from './modules/mapping.nf'
include { EOULSAN_SAM_FILTER } from './modules/filtersam.nf'
include { EOULSAN_EXPRESSION } from './modules/expression.nf'
include { EOULSAN_TSV2XLSX } from './modules/tsv2xlsx.nf'

params.reads = "$projectDir/data/*.fastq.gz"
params.genome = "genome://hg19ens105"
params.annotation = "gtf://hg19ens105"
params.additionalAnnotation = "additionalannotation://hg19ens105"
params.mapperName = "minimap2"
params.mapperVersion = "2.24"
params.mapperFlavor = ""
params.indexerArguments = "-x splice"
params.mappersArguments = "-x splice --eqx --secondary=no --junc-bed /import/rhodos10/ressources/sequencages/bed12/only_chr_Homo_sapiens_ens105.bed"
params.tmpDir = projectDir + "/tmp"
params.binaryDir = "/tmp"
params.storages = read_conf()
params.readFilteringConf = [ "trimpolynend" : "" ]
params.samFilteringConf = [ "removeunmapped" : "true", "quality.threshold" : "1", "removesupplementary": "true", "removemultimatches" : "true" ]
params.expressionConf = [ "genomic.type" : "exon", "attribute.id" : "gene_id", "stranded" : "no", "overlap.mode" : "union", "remove.ambiguous.cases" : "false" ]

Channel
    .of( params.genome )
    .set { genome_ch }

Channel
    .of( params.annotation )
    .set { annot_ch }

Channel
    .of( params.additionalAnnotation )
    .set { add_annot_ch }

Channel
    .fromPath( params.reads, checkIfExists:true )
    .set { reads_ch }

workflow {

    // Index creation
    index_ch = EOULSAN_INDEX(genome_ch, params.mapperName, params.mapperVersion, params.mapperFlavor, params.storages, params.tmpDir, params.binaryDir, params.indexerArguments)

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

    // TSV to XLSX conversion
    expression_ch.combine(add_annot_ch).set { expression_add_annot_combined_ch }
    xlsx_ch = EOULSAN_TSV2XLSX(expression_add_annot_combined_ch, params.tmpDir, params.storages)
}
