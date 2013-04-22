#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
my $program_name = $1 if ( $0 =~ /([^\/]+)$/);
my $usage = <<USAGE;
$program_name -- Gene Expression Error Correction

Author: BENM <binxiaofeng\@gmail.com>
Version: 0.1 alpha
Date: 2010-08-15
Update: 2010-08-18 v0.2 alpha, update data structure
        2010-12-07 v0.3 alpha, add \'--tablelist\'

$program_name [options] [functions]
	--low <int>		set corret CopyNumber, recommand lower than 3, default: 3
	--tablelist <file>	input correct table list file which stored two samples in similar condition
	--hashtable <file>	input hash table file
	--raw <dir>		input directory which stored raw fastq file
	--q_cut <int>		quality cutoff, default: 20
	--correct_bases <int>	correct maximum number of bases[1|2], default: 1
	--enzyme <str>		set digested recognition site with the anchoring enzyme, from 5\' to 3\', default \(17bp, NlaIII\): \'CATG\'; \(16bp, DpnII\): \'GATC\'
	--help			help & usuage

Example:

perl $program_name --hashtable HastTable.xls --raw RawData --correct_bases 2

USAGE

my ($Low,$TableList,$Hashtable,$Raw,$Qual_cutoff,$Correct_bases,$Enzyme,$Help);
my %opts;
GetOptions
(
	\%opts,
	"low:i"=>\$Low,
	"tablelist:s"=>\$TableList,
	"hashtable:s"=>\$Hashtable,
	"raw:s"=>\$Raw,
	"q_cut:i"=>\$Qual_cutoff,
	"correct_bases:i"=>\$Correct_bases,
	"enzyme"=>\$Enzyme,
	"help"=>\$Help
);
die $usage if ((!defined $Hashtable)||(!defined $Raw)||($Help));

$Low ||= 3;
$Qual_cutoff ||= 20;
$Correct_bases ||= 1;
my $tag_len ||= 17;

my @base = qw(A C G T);

my %low_conflict;
my %common;
my %hash;
my @name;
my @total;
my %Name;
my %Stat;
my %Table_list;
if (defined $TableList)
{
	open (TL,$TableList) || die $!;
	while(<TL>)
	{
		s/\s*$//;
		my @t=split /\s+/,$_;
		if (@t==2)
		{
			$Table_list{$t[0]}{$t[1]}=1;
			$Table_list{$t[1]}{$t[0]}=1;
		}
	}
	close TL;
}
open (CK,$Hashtable) || die "Can't open $Hashtable for reading\n";
$_=<CK>;
$_=<CK>;
my @col=split /\s+/,$_;
$tag_len = length($col[0]);
$Enzyme = ((!defined $Enzyme)&&($tag_len == 17)) ? "CATG" : "GATC";
close CK;
open (CON,$Hashtable) || die "Can't open $Hashtable for reading\n";
$_=<CON>;
my @title=split /\t+/,$_;
for (my $ti=1;$ti<@title;$ti++)
{
	if($title[$ti]=~/(\S+)\(CopyNumber/)
	{
		push @name,$1;
		$Name{$1}=$ti-1;
	}
}

while(<CON>)
{
	chomp;
	my $all=$_;
	my @t=split /\t+/,$all;
	for (my $i=1;$i<@t;$i++)
	{
		$t[$i]=~s/[\(\)]//g;
		my ($copy,$order,$tpm)=split /[\;\,]/,$t[$i];
		if (defined $TableList)
		{
			if (exists $Table_list{$name[$i-1]})
			{
				my $T_hash=$Table_list{$name[$i-1]};
				foreach my $k(keys %$T_hash)
				{
					if ((exists $Name{$k})&&($t[$Name{$k}+1]=~/\(0\;/))
					{
						if (($copy>0)&&($copy<=$Low))
						{
							$low_conflict{$i-1}{$t[0]}=$copy;
							$Stat{$name[$i-1]}{"low_conflict1"}++;
						}
					}
					else
					{
						if (!defined $low_conflict{$i-1}{$t[0]})
						{
							$common{$i-1}{$t[0]}=$copy;
							$Stat{$name[$i-1]}{"common1"}++;
						}
					}
				}
			}
		}
		else
		{
			if ($all=~/\(0\;/)
			{
				if (($copy>0)&&($copy<=$Low))
				{
					$low_conflict{$i-1}{$t[0]}=$copy;
					$Stat{$name[$i-1]}{"low_conflict1"}++;
				}
			}
			else
			{
				$common{$i-1}{$t[0]}=$copy;
				$Stat{$name[$i-1]}{"common1"}++;
			}
		}
		$hash{$i-1}{$t[0]}=$copy;
		if ($hash{$i-1}{$t[0]}>0)
		{
			$total[$i-1][0]++;
			$Stat{$name[$i-1]}{"tags_num1"}++;
		}
		$total[$i-1][1]+=$hash{$i-1}{$t[0]};
		$Stat{$name[$i-1]}{"copy_num1"}+=$hash{$i-1}{$t[0]};
	}
}
close CON;

# read RawData
foreach my $sample(@name)
{
	if (-f "$Raw/$sample.fq")
	{
		open (RAW,"$Raw/$sample.fq") || die "Can't find $Raw/$sample.fq for open-reading\n";
	}
	elsif (-f "$Raw/$sample.fq")
	{
		open (RAW,"$Raw/$sample.fastq") || die "Can't find $Raw/$sample.fq for open-reading\n";
	}
	my ($seq,$qual);
	open (C,">$sample.Correct.Tags") || die $!;
	open (F,">$sample.Filter.Tags") || die $!;
	while(<RAW>)
	{
		$seq=<RAW>;
		chomp $seq;
		$seq=$Enzyme.uc(substr($seq,0,$tag_len));
		<RAW>;
		$qual=<RAW>;
		chomp $qual;
		$qual="hhhh".substr($qual,0,$tag_len);
		if (exists $low_conflict{$Name{$sample}}{$seq})
		{
			my @q=split "", $qual;
			my @lowQ=();
			my $q_tot=0;
			for (my $i=0;$i<@q;$i++)
			{
				$q_tot+=ord("$q[$i]")-64;
				if (ord("$q[$i]")-64<$Qual_cutoff)
				{
					push @lowQ,$i;
				}
			}
			next if (($q_tot/@q<$Qual_cutoff)||(@lowQ==0));
			if (@lowQ<=2)
			{
				for (my $i=0;$i<@lowQ;$i++)
				{
					my $tmp=$seq;
					my $origin_base=substr($tmp,$lowQ[$i],1);
					for (my $j=0;$j<@base;$j++)
					{
						next if ($base[$j] eq $origin_base);
						substr($tmp,$lowQ[$i],1,$base[$j]);
						if (exists $common{$Name{$sample}}{$tmp})
						{
							my $correct_low_common_tag_info=("$origin_base:".($lowQ[$i]+1)."-$q[$lowQ[$i]] => $base[$j]");
							print C ("$sample\t$seq\t=>\t",(join "\t",($tmp,$common{$Name{$sample}}{$tmp},$correct_low_common_tag_info)),"\n");
							$Stat{$sample}{"correct1"}++;
							$hash{$Name{$sample}}{$seq}--;
							if ($hash{$Name{$sample}}{$seq}<0)
							{
								$hash{$Name{$sample}}{$seq}=0;
								next;
							}
							$hash{$Name{$sample}}{$tmp}++;
							if ($hash{$Name{$sample}}{$seq}==0)
							{
								$total[$Name{$sample}][0]--;
							}
						}
					}
					if (($Correct_bases==2)&&($i==0)&&(@lowQ==2))
					{
						my $Secondorigin_base=substr($tmp,$lowQ[$i+1],1);
						$tmp=$seq;
						for (my $j=0;$j<@base;$j++)
						{
							next if ($base[$j] eq $origin_base);
							substr($tmp,$lowQ[$i],1,$base[$j]);
							for (my $k=0;$k<@base;$k++)
							{
								next if ($base[$k] eq $Secondorigin_base);
								substr($tmp,$lowQ[$i+1],1,$base[$k]);
								if (exists $common{$Name{$sample}}{$tmp})
								{
									my $correct_low_common_tag_info=("$origin_base:".($lowQ[$i]+1)."-$q[$lowQ[$i]] => $base[$j]; $Secondorigin_base:".($lowQ[$i+1]+1)."-".$q[$lowQ[$i+1]]." => $base[$k]");
									print C ("$sample\t$seq\t=>\t",(join "\t",($tmp,$common{$Name{$sample}}{$tmp},$correct_low_common_tag_info)),"\n");
									$Stat{$sample}{"correct2"}++;
									$hash{$Name{$sample}}{$seq}--;
									if ($hash{$Name{$sample}}{$seq}<0)
									{
										$hash{$Name{$sample}}{$seq}=0;
										next;
									}
									$hash{$Name{$sample}}{$tmp}++;
									if ($hash{$Name{$sample}}{$seq}==0)
									{
										$total[$Name{$sample}][0]--;
									}
								}
							}
						}
					}
				}
			}
			else
			{
				my $low_q;
				for (my $j=0;$j<@lowQ;$j++)
				{
					$low_q.=($lowQ[$j]+1).":$q[$lowQ[$j]] ";
				}
				$Stat{$sample}{"filter"}++;
				print F ("$sample\t$seq\t",(join "; ",$low_q),"\n");
				$hash{$Name{$sample}}{$seq}--;
				if ($hash{$Name{$sample}}{$seq}<0)
				{
					$hash{$Name{$sample}}{$seq}=0;
					next;
				}
				if ($hash{$Name{$sample}}{$seq}==0)
				{
					$total[$Name{$sample}][0]--;
				}
				$total[$Name{$sample}][1]--;
			}
		}
	}
	close RAW;
	close C;
	close F;
}

# Re-Calculate
my %Hash;
foreach my $n(keys %hash)
{
	my $hash_p=$hash{$n};
	my $order=1;
	my $pre=0;
	foreach my $tags(sort{$hash_p->{$b}<=>$hash_p->{$a}}keys %$hash_p)
	{
		if ($hash{$n}{$tags}==0)
		{
			$Hash{$tags}{$n}="\(0;-1;0\)";
		}
		else
		{
			$Stat{$name[$n]}{"tags_num2"}++;
			$Hash{$tags}{$n}="\($hash{$n}{$tags}\;$order\;".(sprintf("%.2f",(1000000*$hash{$n}{$tags}/$total[$n][1])))."\)";
			$order++ unless ($pre==$hash{$n}{$tags});
			$pre=$hash{$n}{$tags};
		}
		$Stat{$name[$n]}{"copy_num2"}+=$hash{$n}{$tags};
	}
}

# Ouput Correction result
open (OUT,">Corrected_$Hashtable") || die $!;
print OUT "Tag-Seq";
foreach (@name)
{
	print OUT "\t$_\(CopyNumber\;Order\;TPM\)";
}
print OUT "\n";
foreach my $i(sort keys %Hash)
{
	my $out=$i;
	my $zero=0;
	for (my $j=0;$j<@name;$j++)
	{
		if ((!exists $Hash{$i}{$j})||($hash{$j}{$i}==0))
		{
			$zero++;
			$out.="\t"."\(0;-1;0\)";
		}
		else
		{
			$out.="\t".$Hash{$i}{$j};
		}
	}
	if ($zero>0)
	{
		map{$Stat{$_}{"low_conflict2"}++ if (($hash{$Name{$_}}{$i}>0)&&($hash{$Name{$_}}{$i}<$Low));}@name;
	}
	else
	{
		map{$Stat{$_}{"common2"}++}@name;
	}
	if ($zero<@name)
	{
		print OUT "$out\n";
	}
}
close OUT;

# Report Statistical report
open (ST,">Corrected_$Hashtable.stat.log") || die $!;
foreach my $sample(@name)
{
	print ST "$sample:\n";
	print ST ("Previous TagsNumber: ",$Stat{$sample}{"tags_num1"}," After Corrected: ",$Stat{$sample}{"tags_num2"},"\n");
	print ST ("Previous Total CopyNumber: ",$Stat{$sample}{"copy_num1"}," After Corrected: ",$Stat{$sample}{"copy_num2"},"\n");
	print ST ("Previous Common TagsNumber: ",$Stat{$sample}{"common1"}," After Corrected: ",$Stat{$sample}{"common2"},"\n");
	print ST ("Previous Low Conflict TagsNumber: ",$Stat{$sample}{"low_conflict1"}," After Corrected: ",$Stat{$sample}{"low_conflict2"},"\n");
	if ((exists $Stat{$sample}{"correct1"})||(exists $Stat{$sample}{"correct2"}))
	{
		$Stat{$sample}{"correct1"}=0 if (!exists $Stat{$sample}{"correct1"});
		$Stat{$sample}{"correct2"}=0 if (!exists $Stat{$sample}{"correct2"});
		print ST ("Corrected one bases: ",$Stat{$sample}{"correct1"}," two bases: ",$Stat{$sample}{"correct2"},"\n");
	}
	else
	{
		print ST ("Corrected one bases: ",$Stat{$sample}{"correct1"}," two bases: ",0,"\n") if (exists $Stat{$sample}{"correct1"});
		print ST ("Corrected one bases: ",0," two bases: ",$Stat{$sample}{"correct2"},"\n") if (exists $Stat{$sample}{"correct2"});
	}
	print ST ("Filtered low quality and low confilt TagsNumber: ",$Stat{$sample}{"filter"},"\n");
	print ST ("Double Check: ",$Stat{$sample}{"tags_num2"},"=",$total[$Name{$sample}][0],"\t",$Stat{$sample}{"copy_num2"},"=",$total[$Name{$sample}][1],"\n");
}
close ST;
