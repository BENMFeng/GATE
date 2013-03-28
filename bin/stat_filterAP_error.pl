#!/usr/bin/perl -w
use strict;

die "perl $0 <filter_adapter.detail>\n" if (@ARGV==0);

my ($totalAln,$totalMis,$totalGap)=(0,0,0);

foreach my $detail(@ARGV)
{
	open (IN,"<$detail") || die $!;
	my $line=0;
	while(<IN>)
	{
		next if (/^\#/);
		chomp;
		my @t=split /\t+/,$_;
		$totalAln+=$t[8];
		$totalMis+=$t[9];
		$totalGap+=$t[10];
		$line++;
		last if ($line>=1000000);
	}
	close IN;
}

print (int(1000000*($totalMis+$totalGap)/$totalAln+0.4999)/10000,"%\n");