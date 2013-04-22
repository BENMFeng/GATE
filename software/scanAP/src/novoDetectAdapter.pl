#!/usr/bin/perl -w
use strict;

use Getopt::Long;
my ($LenRange,$TestNum);
my %opts;
GetOptions(%opts,"len_range:s"=>\$LenRange,"test_num:i"=>\$TestNum);
die "perl $0 <IN:fastq> [-len_range adapter_primers_total_len_range:30-60> [-test_num test_reads_num: 100000]\n" if @ARGV==0;

$LenRange ||= "30-60";
$TestNum  ||= 100000;

my ($kmin,$kmax)=split /\-/,$LenRange;

my @Seq=();
while(<>)
{
	if (/^\@/)
	{
		$_=<>;
		chomp;
		push @Seq,substr($_,length($_)-$kmax,$kmax);
		last if (@Seq>=$TestNum);
		<>;
		<>;
	}
	elsif(/^\>/)
	{
		$_=<>;
		chomp;
		push @Seq,substr($_,length($_)-$kmax,$kmax);
		last if (@Seq>=$TestNum);
	}
	else
	{
		die "it just suitable for FASTQ/FASTA format each sequence in one line\n";
	}
}

my %order=();
my $max=0;
for (my $i=$kmin;$i<=$kmax;$i++)
{
	my %Kmer=();
	foreach my $seq(@Seq)
	{
		if ($i<length($seq))
		{
			for (my $j=0;$j<=length($seq)-$i;$j++)
			{
				my $subseq=substr($seq,$j,$i);
				$Kmer{$subseq}++;
			}
		}
		else
		{
			$Kmer{$seq}++;
		}
	}
	my $n=0;
	my @ary=();
	foreach my $k(sort{$Kmer{$b}<=>$Kmer{$a}}keys %Kmer)
	{
		$n++;
		push @ary,($k,$Kmer{$k});
		last if ($n==2);
	}
	$order{$i}{'drop'}=($ary[1]-$ary[3])/$ary[1];
	$order{$i}{'seq'}=$ary[0];
	$max=($order{$i}{'drop'}>$max)?$order{$i}{'drop'}:$max;
}

my $sec=0;
my @adapter=();
for (my $i=$kmin;$i<=$kmax;$i++)
{
	if ($order{$i}{'drop'} < $max)
	{
		$sec=($order{$i}{'drop'}>$sec)?$order{$i}{'drop'}:$sec;
	}
}

for (my $i=$kmin;$i<=$kmax;$i++)
{
	print qq($order{$i}{'drop'}\t$order{$i}{'seq'});
	if ($order{$i}{'drop'} == $sec)
	{
		push @adapter,$order{$i}{'seq'};
		print "*\n";
	}
	else
	{
		print "\n";
	}
}

print "\nAdapter:\n";
print ((join "\n",@adapter),"\n");

