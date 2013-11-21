#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
my ($Fastq1,$Fastq2,$Help);
GetOptions(\%opts,"fastq1:s"=>\$Fastq1,"fastq2:s"=>\$Fastq2,"help"=>\$Help);

die "perl $0 -fastq1 reads1.fastq -fastq2 reads2.fastq filterlist.files\n" if (!defined $Fastq1 || !defined $Fastq2 || $Help);
my %hash;
if (@ARGV!=0)
{
	while(<>) {
		if (/^\@\S+\:(\d+\:\d+\:\d+\:\d+)/ || /^\@([^\/\s]+)/){
			$hash{$1}=1;
			<>;<>;<>;
			next;
		}
		elsif (/^\>(\S+)/){
			$hash{$1}=1;
		}
	}
}

my $f=0;
my @count=(0,0);
if (keys %hash>0)
{
	foreach my $file($Fastq1,$Fastq2) {
		my $pair=$file;
		$pair=~s/out$/pair\.fq/i;
		$pair=~s/fq$/pair\.fq/i;
		$pair=~s/fastq$/pair\.fastq/i;
		$pair=~s/fq\.gz$/pair\.fq/i;
		$pair=~s/fastq\.gz$/pair\.fastq/i;
		my $single=$file;
		$single=~s/out$/single\.fq/i;
		$single=~s/fq$/single\.fq/i;
		$single=~s/fastq$/single\.fastq/i;
		$single=~s/fq\.gz$/single\.fq/i;
		$single=~s/fastq\.gz$/single\.fastq/i;
		open (IN,$file) if ($file=~/fastq$/i || $file=~/fq$/i);
		open (IN,"zcat $file|") if ($file=~/gz$/);
		open (PE,">$pair") || die $!;
		open (SE,">$single") || die $!;
		my $print=0;
		while(<IN>) {
			if (/^\@\S+\:(\d+\:\d+\:\d+\:\d+)/ || /^\@([^\/\s]+)/){
				if (exists $hash{$1}) {
					$print=0;
				}
				else {
					$print=1;
				}
			}
			if ($print==1) {
				print PE $_;
				$_=<IN>;
				print PE $_;
				$_=<IN>;
				print PE $_;
				$_=<IN>;
				print PE $_;
				$count[$f]++;
			} else {
				print SE $_;
				$_=<IN>;
				print SE $_;
				$_=<IN>;
				print SE $_;
				$_=<IN>;
				print SE $_;
			}
		}
		close IN;
		close PE;
		close SE;
		$f++;
	}
}
if (keys %hash==0 || $count[0] != $count[1] || ($count[0]==0 && $count[1]==0) )
{
	print STDERR "pair reads are not the same or no filterlist.files : $count[0] : $count[1]\n";
	%hash=();
	delete @hash{keys %hash};
	my %pehash;
	open (IN,$Fastq2) || die $!;
	while(<IN>)
	{
		if (/^\@\S+\:(\d+\:\d+\:\d+\:\d+)/ || /^\@([^\/\s]+)/){
			$pehash{$1}++;
		}
		<IN>;
		<IN>;
		<IN>;
	}
	close IN;
	my $i=1;
	foreach my $file($Fastq1,$Fastq2) {
		my $pair=$file;
		$pair=~s/fq$/pair\.fq/i;
		$pair=~s/fastq$/pair\.fastq/i;
		$pair=~s/fq\.gz$/pair\.fq/i;
		$pair=~s/fastq\.gz$/pair\.fastq/i;
		my $single=$file;
		$single=~s/fq$/single\.fq/i;
		$single=~s/fastq$/single\.fastq/i;
		$single=~s/fq\.gz$/single\.fq/i;
		$single=~s/fastq\.gz$/single\.fastq/i;
		open (IN,$file) if ($file=~/fastq$/i || $file=~/fq$/i);
		open (IN,"zcat $file|") if ($file=~/gz$/);
		open (PE,">$pair") || die $!;
		open (SE,">$single") || die $!;
		$i++;
		while(<IN>) {
			my $print=0;
			if (/^\@\S+\:(\d+\:\d+\:\d+\:\d+)/ || /^\@([^\/\s]+)/){
				$pehash{$1}++;
				if (exists $pehash{$1} && $pehash{$1}==$i) {
					$print=1;
				}
				else {
					$print=0;
				}
			}
			if ($print==1) {
				print PE $_;
				$_=<IN>;
				print PE $_;
				$_=<IN>;
				print PE $_;
				$_=<IN>;
				print PE $_;
			} else {
				print SE $_;
				$_=<IN>;
				print SE $_;
				$_=<IN>;
				print SE $_;
				$_=<IN>;
				print SE $_;
			}
		}
		close IN;
		close PE;
		close SE;
	}
}
