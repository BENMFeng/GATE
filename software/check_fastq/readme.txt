1. Compile instruction

Enter into the soure directory.
Type "make", this will generate a runnable program fastqcheck
Type "make install", this will move the executable programs into ../bin/


2. Program instruction

The main C program checkfastq has two main function: (1) caculate distribution of base along reads; (2) caculate distribution of quality along reads.

There is also assist programs to draw figures for the result: distribute_fqcheck.pl draw distribution of base and quality along reads

All the programs here run fast, less than one minute. And only need very small memory, less than 100M.


3. Usage example

(1) run the main C program
  ./checkfastq -i s_1_1.fastq > s_1_1.fqcheck
  ./checkfastq -i s_1_2.fastq > s_1_2.fqcheck

(2) draw distribution figures
  perl ./distribute_fqcheck.pl s_1_1.fqcheck s_1_2.fqcheck -o s_1.fqcheck

