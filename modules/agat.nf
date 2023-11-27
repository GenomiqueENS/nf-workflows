/*
========================================================================================
   AGAT module
========================================================================================
*/

// Parameter definitions
params.OUTPUT = "result/rnabloom"

/*
* AGAT Conversion bed > gtf
*/

process RNABLOOM_AGAT {
   // where to store the results and in which way
   debug true
   publishDir( params.OUTPUT, mode: 'copy' )

   // show in the log which input file is analysed
   tag( "${bloombed}" )
   
   input:
   path bloombed 
   
   output:
   path("rnabloom.transcripts.gtf"), emit: agat_gtf
   
   script:
   """
   /usr/local/bin/agat_convert_bed2gff.pl \
   --bed ${bloombed} \
   -o rnabloom.transcripts.gtf
   """

}


