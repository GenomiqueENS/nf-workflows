# nf-workflows
This bioinformatic pipeline is used for transcript annotation of long-read sequencing data using consensus annotated transcriptome from Isoquant and RNA-Bloom tools.
This workflow is written in nextflow and developped by the GenomiqueENS core facility of the Institute of Biology of the Ecole Normale Superieure (IBENS).

<p align="center" width="100%">
    <img width="50%" src="https://github.com/GenomiqueENS/nf-workflows/assets/91611978/b568caf0-7345-47db-9d1c-cfed6cf79a9b">
</p>


Main workflow

1. Define input section of nextflow.config file:
     - path to nanopore data
     - path to genome annotation
     - path to illumina short reads (optional)
2. Read orientation (Eoulsan)
3. RNA-Bloom subworkflow
     - Concatenate fastq files into single fastq (cat)
     - Transcript annotation (RNA-Bloom)
         - optional polishing with Illumina short-reads
     - Mapping reads to genome (minimap2)
     - Convert sam to bed file (pathools)
     - Convert bed to gtf file (agat)
5. Isoquant subworflow
     - Mapping reads to genome (minimap2)
     - Convert sam to bam files (samtools)
     - Concatenate bam files into single bam (samtools)
     - Transcript annotation (Isoquant)
6. Creation of consensus transcript annotation
   
