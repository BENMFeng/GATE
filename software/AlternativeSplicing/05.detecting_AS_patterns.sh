#!/bin/sh
##############################################################################################################
#Description: get AS event patterns by the 'as' step in AlternatieveSplicing.jar
##############################################################################################################

#Annotate the patterns of alternative splicing patterns
java -jar ./bin/AlternativeSplice.jar -as -i mappingToJunctions.sam -o AS_result.gff -g refgene -s SPresult/spliceSites.txt -t SPresult/transShortRun.txt
#convert and sort the junction read records bam file
samtools view -b -h -S ASresult/JunctionRecord.sam > ASresult/JunctionRecord.bam
samtools sort ASresult/JunctionRecord.bam ASresult/JunctionRecord.sort
samtools index ASresult/JunctionRecord.sort.bam

