/*
========================================================================================
   RNABLOOM_MINIMAP2 module
========================================================================================
*/

// Parameter definitions
params.OUTPUT = "result/rnabloom"

/*
* Transcript genome alignment
*/

process RNABLOOM_MINIMAP2 {
   // where to store the results and in which way
   debug true
   maxForks 1
   publishDir( params.OUTPUT, mode: 'copy' )

   // show in the log which input file is analysed
   tag( "${bloomfasta}" )
   
   input:
   path genome
   path bloomfasta 
   
   output:
   path( "rnabloom_aln.sam" ), emit: rnabloom_sam
   
   script:
   """
   minimap2 -ax splice -uf -k14 \
   ${genome} ${bloomfasta} > rnabloom_aln.sam
   """

}
