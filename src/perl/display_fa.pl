#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my %opts;
my ($Len_th,$Ctg_len,$Display,$Seq_list,$Lst_other,$Substr,$Pattern,$Match,$Redundancy,$Reverse,$Complement,$RC,$Region,$Shift,$FullName,$Help);
GetOptions(\%opts,"len_th:s"=>\$Len_th,"ctg_len:s"=>\$Ctg_len,"display_bases:s"=>\$Display,"seq_list:s"=>\$Seq_list,"list_other:s"=>\$Lst_other,"sub:s"=>\$Substr,"region:s"=>\$Region,"shift:s"=>\$Shift,"pattern:s"=>\$Pattern,"match:s"=>\$Match,"fr:s"=>\$Redundancy,"fullname"=>\$FullName,"reverse"=>\$Reverse,"complement"=>\$Complement,"RC"=>\$RC,"help"=>\$Help);
my $USAGE="USAGE: display_fa.pl <IN:FASTA> [OPTIONS] > out.fa
	-len_th <int-int>		length thesthold min-max
	-ctg_len <int-int>		total consensus length of same contig/scaffold thesthold min-max
	-display_bases <int>		display sequence length in each line, default: 100
	-seq_list <infile>		file was stored sequence name's list
	-list_other <outfile>		file was use to store sequence not exists in the seq_list
	-region <outfile:int[-int][-int]-int-int>	region_file:chr-id-strand-S-E input the column number for \"chr,[id,][strand,]S,E\",
					if defined both \"id\" and \"strand\" column it will join together regions with the same id
	-shift <int-int>		shift posion leftshift-rightshift
	-sub <int-int>			substr sequence, start_pos-end_pos, the left threshold is start from 0
	-pat <str>			pattern string in fasta's name
	-match <str>			find match sequences in fasta and report the location
	-fr <file>			filter redundancy of input sequence, and output redundancy list
	-fullname 			keep full name of sequence
	-reverse			reverse sequence
	-complement			complement sequence
	-RC				reverse & complement

Example:
display_fa.pl ref.fa -len_th 500 -region gene.gff3:0-8-6-3-4 -shift 1

Author: BENM <binxiaofeng\@gmail.com> V1.6.3, 2013-03-01
\n";
die $USAGE if (@ARGV==0||$Help);

print STDERR "Searching $Pattern...\n" if (defined $Pattern);
my ($Len_min,$Len_max) = split /\-/,$Len_th if (defined $Len_th);
my ($Ctg_min,$Ctg_max) = split /\-/,$Ctg_len if (defined $Ctg_len);
$Len_min = (defined $Len_min) ? $Len_min : 0;
$Ctg_min = (defined $Ctg_min) ? $Ctg_min :0;
$Display ||= 100;
my ($S,$E)=(0,0);
my ($LeftShift,$RightShift)=split /\-/,$Shift if (defined $Shift);
$LeftShift = (defined $LeftShift) ? $LeftShift : 0;
$RightShift = (defined $RightShift) ? $RightShift : 0;
if (defined $Substr)
{
	($S,$E) = (split/\-/,$Substr);
	$S+=$RightShift-$LeftShift;
	$E+=$RightShift-$LeftShift;
	$S = 0 if ($S eq "");
	$E = 0 if ($E eq "");
}
my $end = ($E==0) ? "End" : $E;
print STDERR "Searching from $S to $end\n" if (defined $Substr);

my %region;
my $join = 0;
if (defined $Region)
{
	my @R=();
	$R[0]=(split /\:/,$Region)[0];
	push @R,(split /\-/,(split /\:/,$Region)[1]);
	if (@R==6)
	{
		%region=get_region($R[0],$R[1],$R[4],$R[5],$R[3],$R[2]);
		$join = 1;
	}
	elsif (@R==5)
	{
		%region=get_region($R[0],$R[1],$R[3],$R[4],$R[2]);
	}
	elsif (@R==4)
	{
		%region=get_region($R[0],$R[1],$R[2],$R[3]);
	}
	else
	{
		die "-region setting error, Example: regions.txt:0-1-2-3\n";
	}
}

my %SeqList;
if (defined $Seq_list)
{
	open (IN,$Seq_list) || die "Can't open Sequence list file: $Seq_list for reading\n";
	while(<IN>)
	{
		if (/^(\S+)/)
		{
			my $name=$1;
			$name=~s/^\>//;
			$SeqList{$name}=1;
		}
	}
	close IN;
}

my %SeqName=();
my %Lengc=();
my ($chr,$seqname,$seq,$len,$gc_content,$ctg_len)=("","","",0,0,0);
open (LO,">$Lst_other") if (defined $Lst_other);
open (FR,">$Redundancy") if (defined $Redundancy);
while(<>)
{
	s/\s*$//;
	next if ($_ eq "" && !eof);
	if (/>/||(eof))
	{
		$seq .= $_ if ((eof)&&($_!~/^\>/)&&($_!~/^\s/));
		($len,$gc_content,$ctg_len)=cal_lengc($seq) if ($seq ne "");
		if ( ($seq ne "") && ($len>0) && (($len>=$Len_min)&&((!defined $Len_max)||($len<=$Len_max))) && (($ctg_len>=$Ctg_min)&&((!defined $Ctg_max)||($ctg_len<=$Ctg_max))) )
		{
			if ((!defined $Pattern)||((defined $Pattern)&&($chr =~ m/$Pattern/)))
			{
				if (defined $Substr)
				{
					if ($S>=$len)
					{
						if (defined $FullName)
						{
							print STDERR ("\'$seqname\' length $len < subStart $S\n");
						}
						else
						{
							print STDERR ("$chr length $len < subStart $S\n");
						}
						$seq="";
						$len=0;
						if (/^\>(\S+)/)
						{
							$chr=$1;
							$seqname=$_;
						}
						chomp $seqname;
						next;
					}
					$S=((defined $S)&&($S>=0))?$S:0;
					$E=((defined $E)&&($E>0))?$E:$len-1;
					$seq=substr($seq,$S,$E-$S+1);
					$chr.=":$S-$E";
					$seqname.=":$S-$E";
				}
				if ((defined $Redundancy)&&((exists $SeqName{$chr})||(exists $Lengc{"$len-$gc_content"})))
				{
					if ((exists $SeqName{$chr})&&($len!=$SeqName{$chr}))
					{
						if (defined $FullName)
						{
							print FR "different length in \'$seqname\', previous: $SeqName{$chr}, current: $len\n";
						}
						else
						{
							print FR "different length in $chr, previous: $SeqName{$chr}, current: $len\n";
						}
					}
					if (exists $Lengc{"$len-$gc_content"})
					{
						if (defined $FullName)
						{
							print FR "Redundancy: $seqname\n";
						}
						else
						{
							print FR "Redundancy: $chr\n";
						}
						$seq="";
						$len=0;
						if (/^\>(\S+)/)
						{
							$chr=$1;
							$seqname=$_;
							chomp $seqname;
						}
						next;
					}
				}
				if ((defined $Region)&&(exists $region{$chr}))
				{
					if ($join==1)
					{
						my $region_chr=$region{$chr};
						foreach my $k (sort keys %$region_chr)
						{
							my $joinseq="";
							my @region_list="";
							for (my $i=0;$i<@{$region{$chr}{$k}};$i++)
							{
								my ($s,$e,$strand)=@{${$region{$chr}{$k}}[$i]};
								$strand ||= "+";
								push @region_list,"$strand:$s..$e";
								my $subseq=substr($seq,$s,$e-$s+1);
								if ((defined $strand)&&($strand eq "-"))
								{
									$subseq=reverse($subseq);
									$subseq=~tr/ACGTacgt/TGCAtgca/;
								}
								$joinseq.=$subseq;
							}
							my $region_name=join ";",@region_list;
							$region_name=~s/^\;//;
							$region_name=~s/\;$//;
							if (defined $Match)
							{
								if (defined $FullName)
								{
									find_match($joinseq,"$seqname|$region_name|$k");
								}
								else
								{
									find_match($joinseq,"$chr|$region_name|$k");
								}
							}
							else
							{
								if (defined $FullName)
								{
									print "$seqname|$region_name|$k\n";
								}
								else
								{
									print ">$chr|$region_name|$k\n";
								}
								Display(\$joinseq,$Display);
							}
						}
					}
					else
					{
						for (my $i=0;$i<@{$region{$chr}};$i++)
						{
							my ($s,$e,$strand)=@{${$region{$chr}}[$i]};
							$strand ||= "+";
							my $subseq=substr($seq,$s,$e-$s+1);
							if ((defined $strand)&&($strand eq "-"))
							{
								$subseq=reverse($subseq);
								$subseq=~tr/ACGTacgt/TGCAtgca/;
							}
							if (defined $Match)
							{
								if (defined $FullName)
								{
									find_match($subseq,"$seqname:$strand:$s-$e");
								}
								else
								{
									find_match($subseq,"$chr:$strand:$s-$e");
								}
							}
							else
							{
								if (defined $FullName)
								{
									print "$seqname:$strand:$s-$e\n";
								}
								else
								{
									print ">$chr:$strand:$s-$e\n";
								}
								Display(\$subseq,$Display);
							}
						}
					}
				}
				elsif (!defined $Region)
				{
					if (defined $Seq_list)
					{
						my $name=$1 if ($chr=~/(\S+)/);
						if (exists $SeqList{$name})
						{
							if (defined $Match)
							{
								if (defined $FullName)
								{
									find_match($seq,$seqname);
								}
								else
								{
									find_match($seq,$chr);
								}
							}
							else
							{
								if (defined $FullName)
								{
									print "$seqname\n";
								}
								else
								{
									print ">$chr\n";
								}
								Display(\$seq,$Display);
							}
						}
						else
						{
							if (defined $Lst_other)
							{
								if (defined $FullName)
								{
									print LO "$seqname\n$seq\n";
								}
								else
								{
									print LO ">$chr\n$seq\n";
								}
							}
						}
					}
					else
					{
						if (defined $Match)
						{
							if (defined $FullName)
							{
								find_match($seq,$seqname);
							}
							else
							{
								find_match($seq,$chr);
							}
						}
						else
						{
							if (defined $FullName)
							{
								print "$seqname\n";
							}
							else
							{
								print ">$chr\n";
							}
							Display(\$seq,$Display);
						}
					}
				}
				$SeqName{$chr}=$len;
				$Lengc{"$len-$gc_content"}=1;
			}
		}
		if (/^\>(\S+)/)
		{
			$chr=$1;
			$seqname=$_;
			chomp $seqname;
		}
		$seq="";
		$len=0;
	}
	else
	{
		$seq.=$_;
	}
}
close LO if (defined $Lst_other);
close FR;

sub Display 
{
	my $seq_p=shift;
	$$seq_p=reverse($$seq_p) if ((defined $Reverse)||(defined $RC));
	$$seq_p=~ tr/ACGTacgt/TGCAtgca/ if ((defined $Complement)||(defined $RC));
	my $num ||= (@_) ? shift : 100;
	my $disp;
	$$seq_p =~ s/\s$//;
	for (my $i=0;$i<length($$seq_p);$i+=$num) {
		$disp = substr($$seq_p,$i,$num)."\n";
		print $disp;
	}
}

sub get_region
{
	my ($file,$c,$s,$e,$m,$d)=@_;
	my %hash;
	open (IN,$file) || die $!;
	print STDERR "Getting regions: $Region\n";
	while(<IN>)
	{
		chomp;
		my @col=split /\t+/,$_;
		if ((defined $d)&&(!defined $col[$d]))
		{
			my @col=split /\s+/,$_;
		}
		next if ((!defined $col[$s])||(!defined $col[$e])||($col[$s]!~/^\d+$/)||($col[$e]!~/^\d+$/));
		if ($col[$s]>$col[$e])
		{
			my $tmp=$col[$s];
			$col[$s]=$col[$e];
			$col[$e]=$tmp;
		}
		$col[$s]+=$RightShift-$LeftShift;
		$col[$e]+=$RightShift-$LeftShift;
		if (defined $d)
		{
			push @{$hash{$col[$c]}{$col[$d]}},[@col[$s,$e,$m]];
		}
		elsif (defined $m)
		{
			push @{$hash{$col[$c]}},[@col[$s,$e,$m]];
		}
		else
		{
			push @{$hash{$col[$c]}},[@col[$s,$e]];
		}
	}
	close IN;
	return %hash;
}

sub cal_lengc
{
	my $string=shift;
	my $length=length($string);
	$string=~s/[Nn\-]//g;
	my $len_ctg=length($string);
	$string=~s/[ATat]//g;
	my $len_gc=length($string);
	return ($length,$len_gc*100/$length,$len_ctg);
}

sub find_match
{
	my ($seq,$name)=@_;
	my @MatchSeq=();
	while ($Match=~/([ACGTNacgtn]+)/g)
	{
		my $fragment=$1;
		push @MatchSeq,$fragment;
	}
	foreach my $match_seq(@MatchSeq)
	{
		my $Match_len=length($match_seq);
		while($seq=~/($match_seq)/gi)
		{
			my $p=$-[0];
			print ("$name\t+\t$p\t",($p+$Match_len-1),"\t$match_seq\n");
		}
		my $rc_match_seq=reverse($match_seq);
		$rc_match_seq=~tr/ACGTacgt/TGCAtgca/;
		while($seq=~/($rc_match_seq)/gi)
		{
			my $p=$-[0];
			print ("$name\t\-\t$p\t",($p+$Match_len-1),"\t$match_seq\n");
		}
	}
}
