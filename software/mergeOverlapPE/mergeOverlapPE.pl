#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
my ($Prefix,$Seed,$PE_strand,$QFormat,$Qth,$LowQ,$MaxLoq,$Diffbase,$Help);
GetOptions(\%opts,"prefix:s"=>\$Prefix,"seed:i"=>\$Seed,"pe_strand:s"=>\$PE_strand,"qf:s"=>\$QFormat,
	"qth:i"=>\$Qth,"low:i"=>\$LowQ,"maxloq:i"=>\$MaxLoq,"diff:i"=>\$Diffbase,"help"=>\$Help);

my $usage = qq(
mergeOverlapPE.pl -- merge short libary which contained overlap between paired-ends reads
Author : BENM <BinxiaoFeng\@gmail.com>
Version: 0.1.3 alpha 2013-03-20
Birthday : 2013-03-11 0.1.0
Usage: mergeOverlapPE.pl <IN1:PE1.fa|PE1.fq> <IN2:PE2.fa|PE2.fq> <-prefix output> [option]
	--prefix <str>		output prefix name
	--seed <int>		searching seed size, default: 11
	--pe_strand <str>	PE strand, default: FR
	--qf <str>		Quality threshold format, default: L
	--qth <int>		qualtiy threshold, default: 15
	--low <int>		defined the lowest quality, default: 10
	--maxloq <int>		the maximum number of the lowest quality content, default: 5
	--diff <int>		the maximum number of varation between overlaping, default: 3
	--help			help
);

die ($usage) if (@ARGV<2 || !defined $Prefix || defined $Help);
$Seed ||= 11;
$Qth ||= 15;
$LowQ ||= 10;
$MaxLoq ||= 5;
$Diffbase ||= 3;
$QFormat ||= "L";
$PE_strand ||= "FR";
my %qualth=('S'=>33,'X'=>59,'I'=>64,'J'=>66,'L'=>33);
my ($Reads1,$Reads2)=@ARGV;
my $aq1=check_fafq($Reads1);
my $aq2=check_fafq($Reads2);
my $suffix=($aq1==1)?"fasta":"fastq";
die "$Reads1 is not the same format as $Reads2\n" if ($aq1 != $aq2);
my ($Name1,$Name2,$name1,$name2,$seq1,$seq2,$qual1,$qual2)=();
open (R1,$Reads1) if ($Reads1!~/gz$/);
open (R2,$Reads2) if ($Reads2!~/gz$/);
open (R1,"gzip -cd $Reads1|") if ($Reads1=~/gz$/);
open (R2,"gzip -cd $Reads2|") if ($Reads2=~/gz$/);
open (PE1,">$Prefix\_R1.$suffix") || die $!;
open (PE2,">$Prefix\_R2.$suffix") || die $!;
open (MO,">$Prefix\_merged.$suffix") || die $!;
while(<R1>)
{
	if ($aq1==1)
	{
		$name1=$1 if (/^\>(\S+)/);
		$Name1=$_;
		$_=<R1>;chomp;
		$seq1=$_;
		$_=<R2>;
		$name2=$1 if (/^\>(\S+)/);
		$Name2=$_;
		$_=<R2>;chomp;
		$seq2=$_;
		die "Reads name error:\n$name1\n$name2\n" if ($name1 ne $name2);
		if ($PE_strand =~ /F$/)
		{
			$seq2=reverse($seq2);
			$seq2=~tr/ACGTacgt/TGCAtgca/;
		}
		if ($PE_strand =~ /^R/)
		{
			$seq1=reverse($seq1);
			$seq1=~tr/ACGTacgt/TGCAtgca/;
		}
		my ($mo,$seq)=MergeSeq($seq1,$seq2);
		if (defined $mo && defined $seq)
		{
			chomp $Name1;
			print MO "$Name1 $mo\n$seq\n";
		}
		else
		{
			print PE1 "$Name1$seq1\n";
			print PE2 "$Name2$seq2\n";
		}
	}
	else
	{
		$name1=$1 if (/^\@([^\/\s]+)/);
		$Name1=$_;
		$_=<R1>;chomp;
		$seq1=$_;
		$_=<R1>;$_=<R1>;chomp;
		$qual1=$_;
		$_=<R2>;
		$name2=$1 if (/^\@([^\/\s]+)/);
		$Name2=$_;
		$_=<R2>;chomp;
		$seq2=$_;
		$_=<R2>;$_=<R2>;chomp;
		$qual2=$_;
		die "Reads name error:\n$name1\n$name2\n" if ($name1 ne $name2);
		if ($PE_strand =~ /F$/)
		{
			$seq2=reverse($seq2);
			$seq2=~tr/ACGTacgt/TGCAtgca/;
			$qual2=reverse($qual2);
		}
		if ($PE_strand =~ /^R/)
		{
			$seq1=reverse($seq1);
			$seq1=~tr/ACGTacgt/TGCAtgca/;
			$qual1=reverse($qual1);
		}
		my ($mo,$seq,$qual)=MergeSeq($seq1,$seq2,$qual1,$qual2);
		if (defined $mo && $mo ne "" && defined $seq && $seq ne "" && defined $qual && $qual ne "")
		{
			chomp $Name1;
			print MO "$Name1 $mo\n$seq\n\+\n$qual\n";
		}
		else
		{
			print PE1 "$Name1$seq1\n\+\n$qual1\n";
			print PE2 "$Name2$seq2\n\+\n$qual2\n";
		}
	}
}
close R1;
close R2;
close PE1;
close PE2;
close MO;

sub check_fafq
{
	my $file=shift;
	my $answer=0;
	open (IN,$file) if ($file!~/gz$/);
	open (IN,"gzip -cd $file|") if ($file=~/gz$/);
	$_=<IN>;
	die "failed to read file: $file" if (!defined $_ || $_ eq "");
	if ($_=~/^\>/)
	{
		$_=<IN>;
		if ($_=~/^[ACGTN]+$/)
		{
			$answer=1;
		}
		else
		{
			die "\nunidentified format of input sequence!\n";
		}
	}
	elsif ($_=~/^\@/)
	{
		$_=<IN>;
		if ($_=~/^[ACGTN]+$/)
		{
			$answer=2;
		}
		else
		{
			die "\nunidentified format of input sequence!\n";
		}
	}
	else
	{
		die "\nunidentified format of input sequence!\n";
	}
	close IN;
	return $answer;
}

sub MergeSeq
{
	my ($seq1,$seq2,$qual1,$qual2)=@_;
	my ($mo,$seq,$qual);
	my $seq2_rc=reverse($seq2);
	$seq2_rc=~tr/ACGTacgt/TGCAtgca/;
	for (my $i=0;$i<=length($seq2_rc)-$Seed;$i++)
	{
		my $seedSeq=substr($seq2_rc,$i,$Seed);
		if ($seq1=~/.+?($seedSeq)/i)
		{
			my $pos=$-[1];
			if ($pos<length($seq1)-$Seed)
			{
				my $mo=($pos+1)."-".($pos+$Seed)." ".(length($seq2_rc)-$i-$Seed)."-".(length($seq2_rc)-$i).":rc";
				my $f5seq=substr($seq1,0,$pos+$Seed);
				my $f3seq=substr($seq1,$pos+$Seed,length($seq2_rc)-$i-$Seed);
				my $r3seq=substr($seq2_rc,$i+$Seed,length($seq2_rc)-$i-$Seed);
				$seq="$f5seq$r3seq";
				if (defined $qual1 && defined $qual2)
				{
					my $qual2_r=reverse($qual2);
					my $comapre_qual=MaxQual(substr($qual1,$pos+1,$Seed),substr($qual2_r,$i,$Seed));
					return 0 if (!defined $comapre_qual);
					$qual=(substr($qual1,0,$pos)).$comapre_qual.(substr($qual2_r,$i+$Seed,length($seq2_rc)-$i-$Seed));
				}
				if (length($seq)<=length($seq1))
				{
					if (uc($f3seq) ne uc($r3seq))
					{
						my ($diff,$n1,$n2)=(0,0,0);
						my @f3=split "",uc($f3seq);
						my @r3=split "",uc($r3seq);
						for (my $j=0;$j<@f3;$j++)
						{
							$diff++ if ($f3[$j] ne $r3[$j] && $f3[$j] ne "N" && $r3[$j] ne "N");
							$n1++ if ($f3[$j] eq "N");
							$n2++ if ($r3[$j] ne "N");
							last if ($diff>$Diffbase);
						}
						return 0 if ($diff>$Diffbase);
						if (defined $qual1 && defined $qual2)
						{
							my $qual2_r=reverse($qual2);
							my $f3qual=substr($qual1,$pos+$Seed,length($seq2_rc)-$i-$Seed);
							my $r3qual=substr($qual2_r,$i+$Seed,length($seq2_rc)-$i-$Seed);
							my $compare_qual=MaxQual($f3qual,$r3qual);
							return 0 if (!defined $compare_qual);
							if ($compare_qual eq $f3qual)
							{
								$seq="$f5seq$f3seq";
								$qual=substr($qual1,0,$pos+length($seq2_rc)-$i);
							}
							else
							{
								$qual=(substr($qual1,0,$pos+$Seed)).$r3qual;
							}
						}
						elsif ($n2>$n1)
						{
							$seq = "$f5seq$f3seq";
						}
					}
					else
					{
						if (defined $qual1 && defined $qual2)
						{
							my $qual2_r=reverse($qual2);
							my $f3qual=substr($qual1,$pos+$Seed,length($seq2_rc)-$i-$Seed);
							my $r3qual=substr($qual2_r,$i+$Seed,length($seq2_rc)-$i-$Seed);
							my $compare_qual=MaxQual($f3qual,$r3qual);
							return 0 if (!defined $compare_qual);
							if ($compare_qual eq $f3qual)
							{
								$qual=substr($qual1,0,$pos+length($seq2_rc)-$i);
							}
							else
							{
								$qual=(substr($qual1,0,$pos+$Seed)).$r3qual;
							}
						}
					}
				}
				else
				{
					my $f_ol=substr($seq1,$pos,length($seq1)-$pos);
					my $r_ol=substr($seq2_rc,$i,length($seq1)-$pos);
					if (uc($f_ol) ne uc($r_ol))
					{
						my ($diff,$n1,$n2)=(0,0,0);
						my @fo=split "",uc($f_ol);
						my @ro=split "",uc($r_ol);
						for (my $j=0;$j<@fo;$j++)
						{
							$diff++ if ($fo[$j] ne $ro[$j] && $fo[$j] ne "N" && $ro[$j] ne "N");
							$n1++ if ($fo[$j] eq "N");
							$n2++ if ($ro[$j] ne "N");
							last if ($diff>$Diffbase);
						}
						return 0 if ($diff>$Diffbase);
						if (defined $qual1 && defined $qual2)
						{
							my $qual2_r=reverse($qual2);
							my $foqual=substr($qual1,$pos,length($seq1)-$pos);
							my $roqual=substr($qual2_r,$i,length($seq1)-$pos);
							my $compare_qual=MaxQual($foqual,$roqual);
							return 0 if (!defined $compare_qual);
							if ($compare_qual eq $foqual)
							{
								$seq=$seq1.(substr($seq2_rc,$i+length($seq1)-$pos,length($seq2_rc)-$i+length($seq1)+$pos));
								$qual=$qual1.(substr($qual2_r,$i+length($qual1)-$pos,length($qual2_r)-$i+length($qual1)+$pos));
							}
							else
							{
								$qual=(substr($qual1,0,$pos)).substr($qual2_r,$i,length($qual2_r)-$i);
							}
						}
						elsif ($n2>$n1)
						{
							$seq=$seq1.(substr($seq2_rc,$i+length($seq1)-$pos,length($seq2_rc)-$i+length($seq1)+$pos));
						}
					}
					else
					{
						if (defined $qual1 && defined $qual2)
						{
							my $qual2_r=reverse($qual2);
							my $foqual=substr($qual1,$pos,length($seq1)-$pos);
							my $roqual=substr($qual2_r,$i,length($seq1)-$pos);
							my $compare_qual=MaxQual($foqual,$roqual);
							return 0 if (!defined $compare_qual);
							if ($compare_qual eq $foqual)
							{
								$qual=$qual1.(substr($qual2_r,$i+length($qual1)-$pos,length($qual2_r)-$i+length($qual1)+$pos));
							}
							else
							{
								$qual=(substr($qual1,0,$pos)).substr($qual2_r,$i,length($qual2_r)-$i);
							}
						}
					}
				}
				return ($mo,$seq,$qual);
			}
		}
	}
}

sub MaxQual
{
	my ($q1,$q2)=@_;
	my $q="";
	my @qual1=split "",$q1;
	my @qual2=split "",$q2;
	my $qq1=0;
	my $qq2=0;
	my $lowq1=0;
	my $lowq2=0;
	for (my $i=0;$i<@qual1;$i++)
	{
		my $qval1=ord($qual1[$i])-$qualth{$QFormat};
		my $qval2=ord($qual2[$i])-$qualth{$QFormat};
		$qq1+=$qval1;
		$lowq1++ if ($qval1<=$LowQ);
		$qq2+=$qval2;
		$lowq2++ if ($qval2<=$LowQ);
	}
	if (@qual2==0)
	{
		return $q1;
	}
	elsif (@qual1==0)
	{
		return $q2;
	}
	elsif (($qq1/@qual1>=$Qth || $qq2/@qual2>=$Qth)&&($lowq1<$MaxLoq || $lowq2<$MaxLoq))
	{
		$q=($qq1>$qq2)?$q1:$q2;
		return $q;
	}
}