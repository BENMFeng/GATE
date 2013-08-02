#!/usr/bin/perl -w
#Author: BENM <binxiaofeng@gmail.com>
#v0.1.1, 2013-03-05
use strict;
use Getopt::Long;
my %opts;
my ($R1,$R2,$R3,$Format,$Proper,$Help);
GetOptions(\%opts,"R1:s"=>\$R1,"R2:s"=>\$R2,"R3:s"=>\$R3,"format:s"=>\$Format,"proper"=>\$Proper,"help"=>\$Help);
die "perl $0 <IN:SAM|pipe_SAM_out> [-R1 R1.fastq] [-R2 R2.fastq] [-R3 R3.fastq] [-f fq|fa]\n" if ($Help);
$Format ||= 'fq';
open (R1,">$R1") if (defined $R1);
open (R2,">$R2") if (defined $R2);
open (R3,">$R3") if (defined $R3);
my $k=0;
while(<>){
	my @t=split /\t+/,$_;
	next if ($Proper && ($t[1] & 0x2)==0);
	if (($t[1] & 0x80) == 128 && defined $R2)
	{
		if ($Format eq "fq")
		{
			print R2 "\@$t[0]\/2\n$t[9]\n+\n$t[10]\n";
		}
		else
		{
			print R2 "\>$t[0]\/2\n$t[9]\n";
		}
	}
	elsif ( ( ($t[1] & 0x40)== 64 ) && defined $R1)
	{
		if ($Format eq "fq")
		{
			print R1 "\@$t[0]\/1\n$t[9]\n+\n$t[10]\n";
		}
		else
		{
			print R1 "\>$t[0]\/1\n$t[9]\n";
		}
	}
	elsif (defined $R3)
	{
		if ($Format eq "fq")
		{
			print R3 "\@$t[0]\n$t[9]\n+\n$t[10]\n";
		}
		else
		{
			print R3 "\>$t[0]\n$t[9]\n";
		}
	}
	else
	{
		if ($Format eq "fq")
		{
			print "\@$t[0]";
			if (($t[1] & 0x80) == 128)
			{
				print "\/2\n";
			}
			elsif (($t[1] & 0x40)== 64)
			{
				print "\/1\n";
			}
			else
			{
				print "\n";
			}
			print "$t[9]\n+\n$t[10]\n";
		}
		else
		{
			print "\>$t[0]";
			if (($t[1] & 0x80) == 128)
			{
				print "\/2\n";
			}
			elsif (($t[1] & 0x40)== 64)
			{
				print "\/1\n";
			}
			else
			{
				print "\n";
			}
			print "$t[9]\n";
		}
	}
}
close R1;
close R2;
close R3;