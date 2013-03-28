#!/bin/sh
##############################################################################################################
#Description: Example shell script to get the unmapped read records from the BAM file which contains reads 
#             mapped to the genome, then covert the unmapped read records to fastq format read sequences.          
##############################################################################################################

#Make sure that you have installed the Samtools and Java with version larger than 1.6
samtools view -f 4 -bh mapping_to_genome.bam > noMappedToGenome.bam
java -jar SamToFastq.jar VALIDATION_STRINGENCY=SILENT I=noMappedToGenome.bam F=noMappedToGenome.fq

