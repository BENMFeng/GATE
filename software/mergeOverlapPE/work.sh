perl mergeOverlapPE.pl PE1.fq PE2.fq --prefix PE
bwa pemerge PE1.fq PE2.fq > PE_merged_by_bwa.fastq
bwa pemerge -u  PE1.fq PE2.fq  -Q 100 > PE_unmerged.fastq
