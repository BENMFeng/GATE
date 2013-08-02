#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opts;
my ($Uniq,$Num,$Help);
GetOptions(\%opts,"uniq"=>\$Uniq,"t"=>\$Num,"help"=>\$Help);
die "perl $0 <IN:SAM|pipe_SAM_out> [-t use number instead reads name] [-help]\n" if ($Help);
my $k=0;
my %hash;
while(<>){
	my @t=split;
	next if (($t[1] & 0x4) >0);
	next if ((defined $Uniq) && ($_!~/XT\:A\:U/ || $_=~/XA\:Z\:/));
	if (/MD:Z\:(\S+)/){
		my $tag=$1;
		my $match=0;
		while($tag=~/(\d+)/g){
			$match+=$1;
		}
		$hash{length($t[9])-$match}++;
	}
}
print "MismatchOfReads\tReadsNumber\n";
foreach my $m(sort{$a<=>$b}keys %hash) {
	print "$m\t$hash{$m}\n";
}