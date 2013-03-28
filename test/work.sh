../software/scanAP/scanAP -i test.fastq -a ap.fa -s stat.txt -d detail.txt -k 10 -e 3 -m 10 -l
perl ../bin/trim_seq.pl -edge 35 -len_t 35 -len_p 0.25 --trim_mode both --ap_mode --ap_db ap.fa --verbose -trim_detail detail.txt test.fastq
