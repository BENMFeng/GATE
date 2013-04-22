#!/usr/bin/perl -w
use strict;
die "perl $0 <IN:Seq.fa>\n" if @ARGV==0;

my $Seq="";
while(<>)
{
	s/\s*$//;
	if (/^\>/ || eof)
	{
		$Seq.=$_ if (!/^>/ && $_ ne "");
		trim_polymer($Seq) if ($Seq ne "" && length($Seq)>0);
		$Seq="";
	}
	else
	{
		$Seq.=$_;
	}
}

sub trim_polymer
{
	my $seq=shift;
	my $polymerlen ||= 7;
	while ($seq=~/([^N]{1})(\1{$polymerlen,})/g){
		my ($element,$repeat,$rail,$pos)=($1,$2,$',$-[1]);
		my $primary_element=$element;
		while ($element=~/([ACGT]{1,}?)(\1{1,})/g)
		{
			$primary_element=$1;
			if ($rail=~/^(($primary_element)+)/)
			{
				$repeat.=$1;
			}
		}
		my $len=length("$element$repeat");
		print ("$primary_element\t$len\t$pos-".($pos+$len-1),"\n");
	}
}