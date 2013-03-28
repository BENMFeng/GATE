#!/bin/sh
#################################################################################################################################
#Description: Example shell script used for the mapping of short reads (fastq format) to reference genome sequence (fasta format) 
#             with BWA, The BWA and samtools should be installed in your system and in the system path. 
#################################################################################################################################


#single-end
#Mapping with BWA, modify the number of processors to run the analysis by yourself
bwa aln -t 10 -n 2 refGenome fq > fq.sai
bwa samse refGenome fq.sai fq > mapping_to_genome_se.sam
#Convert SAM to BAM, sort and index BAM
samtools view -b -h -S mapping_to_genome_se.sam > mapping_to_genome_se.bam
samtools sort mapping_to_genome_se.bam mapping_to_genome_se.sort
samtools index mapping_to_genome_se.sort.bam


#pair-end
#Mapping with BWA, modify the number of processors to run the analysis by yourself
bwa aln -t 10 -n 2 refGenome fq1 > fq1.sai
bwa aln -t 10 -n 2 refGenome fq2 > fq2.sai
bwa sampe -s refGenome fq1.sai fq2.sai fq1 fq2 > mapping_to_genome_pe.sam
#Convert SAM to BAM, sort and index BAM
samtools view -b -h -S mapping_to_genome_pe.sam > mapping_to_genome_pe.bam
samtools sort mapping_to_genome_pe.bam mapping_to_genome_pe.sort
samtools index mapping_to_genome_pe.sort.bam

