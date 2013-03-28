#!/bin/sh
##########################################################################################################
#Description: example shell script used to convert, sort and index the sam format file to sorted bam format 
#             file with Samtools
##########################################################################################################

samtools view -b -h -S mapping.sam > mapping.bam
samtools sort mapping.bam mapping.sort
samtools index mapping.sort
