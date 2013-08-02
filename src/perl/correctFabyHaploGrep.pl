#!/usr/local/bin/perl -w
use strict;

die "perl $0 <IN1:haplogrep.hsd> <IN2:MT.fa> <SM> [IN3:ref_MT.fa]\n" if @ARGV<3;

my ($Haplo,$Fasta,$SM,$ref)=@ARGV;
my %hash;

open (HP,$Haplo) || die $!;
while(<HP>)
{
	s/\s*$//;
	my @t=split /\s+/,$_;
	next if ($t[0] !~ /$SM/);
	for (my $i=3;$i<@t;$i++)
	{
		$t[$i]="$1$2" if ($t[$i]=~/(\d+)\.\d+([ACGT])/);
		$hash{$1-1}=$2 if ($t[$i]=~/(\d+)([ACGT])/);
		$hash{$1-1}="N" if ($t[$i]=~/(\d+)d/);
	}
}
close HP;

my @REF=();
open (R,$ref) || die $!;
while(<R>)
{
	s/\s*$//;
	if (!/^\>/)
	{
		next if ($_ eq "");
		push @REF,(split //,$_);
	}
}
close R;

my ($name,$seq)=("","");
open (FA,$Fasta) || die $!;
while(<FA>)
{
	s/\s*$//;
	if (/^(>.*)$/ || eof)
	{
		if (eof)
		{
			$seq.=$_;
		}
		if (length($seq)>0)
		{
			print ">$SM\_MT\n";
			correct($seq);
		}
		$name=$1;
	}
	else
	{
		$seq.=$_;
	}
}
close FA;

sub correct
{
	my $str=shift;
	my @s=split //,$str;
	for (my $i=0;$i<@s;$i++)
	{
		if (exists $hash{$i})
		{
			print STDERR ($i+1,"\t$REF[$i]\t$hash{$i}\tSNV\n");
			print STDERR ($i+1,"\t$s[$i]\t=>\t$hash{$i}\tcorrected\n") if (uc($s[$i]) ne $hash{$i});
			$s[$i]=$hash{$i};
		}
		else
		{
			if (defined $ref && $s[$i]!~/N/i && $REF[$i] ne uc($s[$i]))
			{
				print STDERR ($i+1,"\t$s[$i]\t=>\t$REF[$i]\tcorrectedByRef\n");
				$s[$i]=$REF[$i];
			}
		}
	}
	my $out=join '',@s;
	$out=~s/\-//;
	Display_fa($out);
}

sub Display_fa
{
	my $seq_p=shift;
	my $num ||= (@_) ? shift : 100;
	my $disp;
	$seq_p =~ s/\s$//;
	for (my $i=0;$i<length($seq_p);$i+=$num) {
		$disp = substr($seq_p,$i,$num)."\n";
		print $disp;
	}
}
