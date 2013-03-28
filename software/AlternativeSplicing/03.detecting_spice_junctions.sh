#!/bin/sh
##############################################################################################################
#Description: get candidate splice junctions by the 'sp' step in AlternativeSplicing.jar
##############################################################################################################

#Make sure that you have installed the Samtools and Java with version larger than 1.6
java -jar AlternativeSplicing.jar -sp -i mapping_to_genome.sort.bam -r refgenome_sequence -g gff_refgene
