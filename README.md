# nf-workflows
This bioinformatic pipeline is used for transcript annotation of long-read sequencing data using consensus annotated transcriptome from Isoquant and RNA-Bloom tools.
This workflow is written in nextflow and developped by the GenomiqueENS core facility of the Institute of Biology of the Ecole Normale Superieure (IBENS).

![pipeline_transcript_annot](https://github.com/GenomiqueENS/nf-workflows/assets/91611978/76114ed4-c957-4ce5-96ab-f46720440ce3 |width=100)

Main workflow

1. Define input section of nextflow.config file:
     - Nanopore data
     - Genome annotation
     - Illumina short reads (optional)
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
   
