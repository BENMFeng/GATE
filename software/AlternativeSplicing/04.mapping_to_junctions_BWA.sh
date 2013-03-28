#!/bin/sh
###################################################################################################
#Description: example shell script used for the mapping of the nomapped reads (fastq format) 
#             to splice junction sequences with BWA      
###################################################################################################

#mapping with BWA, make sure that the bwa have been installed
bwa index -a bwtsw SPresult/spliceJunctions.fa
bwa aln -t 12 -n 2 SPresult/spliceJunctions.fa noMappedToGenome.fq > mappingToJunctions.sai
bwa samse SPresult/spliceJunctions.fa mappingToJunctions.sai noMappedToGenome.fq > mappingToJunctions.sam