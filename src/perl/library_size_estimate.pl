#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my %opts;
my ($Type,$Help);
GetOptions(\%opts,"t:s"=>\$Type,"help"=>\$Help);
die qq(Example: samtools view -F 4 aln.bam | perl $0 -t sam > aln.insert-size.dist.txt
perl $0 reads.fastq -t fq > reads.length.dist.txt\n) if ($Help);
$Type ||= "sam";
if ($Type =~ /fq/i || $Type =~ /fastq/i)
{
	my $cmd=qq(awk 'NR\%4==2{print length\(\$1\)}' $ARGV[0] | sort -n -k 1 | uniq -c | awk '{print \$2"\\t"\$1}');
	system $cmd;
}
elsif ($Type =~ /sam/i || $Type =~ /bam/i)
{
	my %hash;
	if (@ARGV==0)
	{
		while(<>)
		{
			my @t=split;
			if (abs($t[8])>=length($t[9])){
				$hash{abs($t[8])}++;
			}else{
				$hash{length($t[9])}++;
			}
		}
	}
	else
	{
		open (IN,$ARGV[0]) if ($Type =~ /sam/i);
		open (IN,"samtools view $ARGV[0]|") if ($Type =~ /bam/);
		while(<IN>)
		{
			my @t=split;
			if (abs($t[8])>=length($t[9])){
				$hash{abs($t[8])}++;
			}else{
				$hash{length($t[9])}++;
			}
		}
		close IN;
	}
	foreach my $ins(sort{$a<=>$b}keys %hash)
	{
		print "$ins\t$hash{$ins}\n";
	}
}

