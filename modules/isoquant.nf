/*
========================================================================================
   ISOQUANT module
========================================================================================
*/

// Parameter definitions
params.OUTPUT = "result/isoquant"

/*
* Create Isoquant trancript model
*/

process ISOQUANT {

   // where to store the results and in which way
   publishDir( params.OUTPUT, mode: 'copy' )

   // show in the log which input file is analysed
   debug true
   tag( "${bam}" )

   input:
   path genome 
   tuple path(bam), path(bai)

   output:
   path( "result/isoquant/OUT/*.gtf" ), emit: isoquant_gtf

   
   script:
   """
   isoquant.py --reference ${genome}    \
   --bam ${bam} \
   --data_type nanopore    \
   --clean_start    \
   --model_construction_strategy assembly \
   -t 12 \
   -o ${params.OUTPUT}
   """
}  


