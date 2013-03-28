#!/usr/bin/perl -w
use strict;
die "perl $0 <IN1:sorted.bam> <IN2:gene.gff>\n" if @ARGV!=2;
my ($bam,$gff)=@ARGV;
my $print=0;
open (IN,$gff) || die $!;
while(<IN>){
	my @t=split;
	if ($t[2] eq "mRNA") {
		my $out=`samtools view $bam $t[0]:$t[3]-$t[4]`;
		if (defined $out && $out ne "") {
			$print=1;
		} else {
			$print=0;
		}
	}
	print if ($print==1);
}
close IN;
