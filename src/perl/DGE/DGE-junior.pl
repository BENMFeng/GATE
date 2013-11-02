#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $program_name = $1 if ( $0 =~ /([^\/]+)$/);
my $usage = <<USAGE;
$program_name -- Gene Expression Hash Table and Differentially Expressed Comparision Analysis

Author: BENM <binxiaofeng\@gmail.com>
Version: 0.4 alpha
Date:   0.1 alpha 2010-07-07
Update: 0.2 alpha 2010-08-12
        0.3 alpha 2010-12-08 add "--compare_list"
        0.4 alpha 2010-1228 add "--plot"

$program_name [IN:copynumber files] [options]	[functions]
		--dir <directory>	input directory which contains CopyNumber(with suffix of ".CopyNumber") files
		--compare_list <file>	input file which contain samples in pairs (split by TAB) which need to compare
		--min <int>		set minimum copy number cutoff, default: 0
		--max <int>		set maximum copy number cutoff, default: 10000
		--out <file>		output hash table files
		--order <int>		output hash table sort by order
		--hash_table <file>	input hash table which want to be sorted
		--diff			output differential anaylsis files for each two TagSeq items
		--ranking <int>		ranking diff analysis files, type:1,2,3,4; default:2
		--sort <col1:col2:col3...>	sort hash table by different TagSeq items, it should require msort
		--plot			plot expression enrichment relativity figure
		--help			help & usuage

Example:
perl $program_name 1.CopyNubmer 2.CopyNumber 3.CopyNumber --out TagSeq.HashTable.xls
perl $program_name --dir ./ --min 1 --out TagSeq.HashTable.xls --diff
perl $program_name --sort 1,2,3 --hash_table TagSeq.HashTable.xls --diff

USAGE
my ($Dir,$CompList,$Min,$Max,$Out,$Order,$Sort,$Diff,$Rank,$Hashtable,$Plot,$Help);
my %opts;
GetOptions
(
	\%opts,
	"dir:s"=>\$Dir,
	"compare_list:s"=>\$CompList,
	"min:i"=>\$Min,
	"max:i"=>\$Max,
	"out:s"=>\$Out,
	"order"=>\$Order,
	"sort:s"=>\$Sort,
	"ranking:i"=>\$Rank,
	"hash_table:s"=>\$Hashtable,
	"diff"=>\$Diff,
	"plot"=>\$Plot,
	"help"=>\$Help
);

die $usage if (((@ARGV==0)&&(!defined $Dir)&&(!defined $Diff))||($Help));
$Min ||= 0;
$Max ||= 10000;
$Rank ||= 2;
die "--ranking type is one of \"1,2,3\"\n" if (($Rank!=1)&&($Rank!=2)&&($Rank!=3)&&($Rank!=4));
$Rank += 2;
my @File;

if (@ARGV!=0)
{
	map{push @File,$_ if (-f $_)}@ARGV;
}
else
{
	if (defined $Dir)
	{
		@File=glob("$Dir/*CopyNumber");
	}
}

my %Compare;
if (defined $CompList)
{
	open (CP,$CompList) || die $!;
	while(<CP>)
	{
		chomp;
		my @t=split;
		if (@t==2)
		{
			$Compare{$t[0]}{$t[1]}=1;
		}
	}
	close IN;
}

if (@File>0)
{
	&copynumber2hashtable;
}

if (defined $Sort)
{
	msort($Out) if (defined $Out);
	msort($Hashtable) if (defined $Hashtable);
}

if (defined $Diff)
{
	diff($Out) if (defined $Out);
	diff($Hashtable) if (defined $Hashtable);
}

sub copynumber2hashtable
{
	my %hash;
	my $n=0;
	foreach my $f(@File)
	{
		open (IN,"<$f") || die $!;
		while (<IN>)
		{
			chomp;
			my @t=split;
			next if (($t[0]=~/[^ACGTN]/i)||($t[1]!~/\d+/)||($t[1]<$Min)||($t[1]>$Max));
			$hash{$t[0]}{$n} = "(".(join (";",@t[1..(@t-1)])).")";
		}
		close IN;
		$n++;
	}
	
	for(my $k=0;$k<@File;$k++)
	{
		$File[$k]=(split /\//,$File[$k])[-1];
		$File[$k]=~s/txt//;
		$File[$k]=~s/[\_\-\.]$//;
		$File[$k]=~s/CopyNumber//;
		$File[$k]=~s/[\_\-\.]$//;
		$File[$k]=~s/Tag//;
		$File[$k]=~s/[\_\-\.]$//;
	}
	open (OUT,">$Out") || die $!;
	my $title=("Tag-Seq\t".((join "(CopyNumber;Order;TPM)\t",@File)."(CopyNumber;Order;TPM)"));
	print OUT "$title\n";
	foreach my $str(sort keys %hash)
	{
		my $out=$str;
		for (my $j=0;$j<@File;$j++)
		{
			if (!exists $hash{$str}{$j})
			{
				$out.="\t(0;-1;0)";
			}
			else
			{
				$out.="\t".$hash{$str}{$j};
			}
		}
		print OUT "$out\n";
	}
	close OUT;
}

sub msort
{
	my $file=shift;
	my @sort_col=split /[\:\;\-\_\,]/,$Sort;
	my $k;
	map{$k.="mr".($_*2).","}@sort_col;
	$k=~s/\,$//;
	my $msort=qq(msort -t '\t\(\,\;\)' -k '$k' $file > $file.tmp && mv $file.tmp $file);
	system $msort;
}

sub diff
{
	my $file=shift;
	open (IN,$file) || die $!;
	$_=<IN>;
	chomp;
	my @col=split /\s+/,$_;
	if ($col[0]!~/Tag/)
	{
		warn "no title\n";
		exit();
	}
	close IN;
	for (my $i=1;$i<@col;$i++)
	{
		my $A=$1 if ($col[$i]=~/(\S+)[\(\_]/);
		next if ((defined $CompList)&&(!exists $Compare{$A}));
		for (my $j=1;$j<@col;$j++)
		{
			next if ($i==$j);
			my $B=$1 if ($col[$j]=~/(\S+)[\(\_]/);
			next if ((defined $CompList)&&(!exists $Compare{$A}{$B}));
			open (OUT1,">$A\_ConflictExp\_$B.xls") || die $!;
			print OUT1 "$col[0]\t$A\_CopyNumber\t$B\_CopyNumber\t\($A-$B\)\t$A\/$B\t\($A-$B\)\/$A\tlog\($A\)-log\($B\)\n";
			open (OUT2,">$A\_CommonExp\_$B.xls") || die $!;
			print OUT2 "$col[0]\t$A\_CopyNumber\t$B\_CopyNumber\t\($A-$B\)\t$A\/$B\t\($A-$B\)\/$A\tlog\($A\)-log\($B\)\n";
			open (OUT3,">$A\_DiffExp\_$B.xls") || die $!;
			print OUT3 "$col[0]\t$A\_CopyNumber\t$B\_CopyNumber\t\($A-$B\)\t$A\/$B\t\($A-$B\)\/$A\tlog\($A\)-log\($B\)\n";
			my (@con,@com,@dif);
			open (IN,$file) || die $!;
			<IN>;
			while(<IN>)
			{
				chomp;
				my @t=split /\s+/,$_;
				next if (($t[$i]=~/\-1/)&&($t[$j]=~/\-1/));
				my ($I,$J)=($t[$i],$t[$j]);
				$I=~s/[\(\)]//g;
				$J=~s/[\(\)]//g;
				my ($C1,$O1,$T1)=split /[\,\;]/,$I;
				my ($C2,$O2,$T2)=split /[\,\;]/,$J;
				next if ((($C1>0)&&(($C1<$Min)||($C1>$Max)))||(($C2>0)&&(($C2<$Min)||($C2>$Max))));
				$O1=($O1==-1)?10000000:$O1;
				$O2=($O2==-1)?10000000:$O2;
				$T1=(!defined $T1 || $T1==0)?1/10000000:$T1;
				$T2=(!defined $T2 || $T2==0)?1/10000000:$T2;
				if (($t[$i]=~/\(0\;/)||($t[$j]=~/\(0\;/))
				{
					next if (($t[$i]=~/\(0\;/)&&($t[$j]=~/\(0\;/));
					if (defined $Order)
					{
						push @con,[$t[0],$C1,$C2,($O1-$O2),$O1/$O2,($O1-$O2)/$O1,(log($O1)-log($O2))];
						push @dif,[$t[0],$C1,$C2,($O1-$O2),$O1/$O2,($O1-$O2)/$O1,(log($O1)-log($O2))];
					}
					else
					{
						push @con,[$t[0],$C1,$C2,($T1-$T2),$T1/$T2,($T1-$T2)/$T1,(log($T1)-log($T2))];
						push @dif,[$t[0],$C1,$C2,($T1-$T2),$T1/$T2,($T1-$T2)/$T1,(log($T1)-log($T2))];
					}
				}
				else
				{
					if (defined $Order)
					{
						push @com,[$t[0],$C1,$C2,($O1-$O2),$O1/$O2,($O1-$O2)/$O1,(log($O1)-log($O2))];
						push @dif,[$t[0],$C1,$C2,($O1-$O2),$O1/$O2,($O1-$O2)/$O1,(log($O1)-log($O2))];
					}
					else
					{
						push @com,[$t[0],$C1,$C2,($T1-$T2),$T1/$T2,($T1-$T2)/$T1,(log($T1)-log($T2))];
						push @dif,[$t[0],$C1,$C2,($T1-$T2),$T1/$T2,($T1-$T2)/$T1,(log($T1)-log($T2))];
					}
				}
			}
			close IN;
			if (defined $Order)
			{
				foreach my $x(sort{$a->[$Rank]<=>$b->[$Rank]}@con)
				{
					print OUT1 (join "\t",@{$x},"\n");
				}
				foreach my $y(sort{$a->[$Rank]<=>$b->[$Rank]}@com)
				{
					print OUT2 (join "\t",@{$y},"\n");
				}
				foreach my $z(sort{$a->[$Rank]<=>$b->[$Rank]}@dif)
				{
					print OUT3 (join "\t",@{$z},"\n");
				}
			}
			else
			{
				foreach my $x(sort{$b->[$Rank]<=>$a->[$Rank]}@con)
				{
					print OUT1 (join "\t",@{$x},"\n");
				}
				foreach my $y(sort{$b->[$Rank]<=>$a->[$Rank]}@com)
				{
					print OUT2 (join "\t",@{$y},"\n");
				}
				foreach my $z(sort{$b->[$Rank]<=>$a->[$Rank]}@dif)
				{
					print OUT3 (join "\t",@{$z},"\n");
				}
			}
			close OUT1;
			close OUT2;
			close OUT3;
			plot_relativity("$A\_DiffExp\_$B.xls") if (defined $Plot);
		}
	}
}

sub plot_relativity
{
	my $file=shift;
	my ($A,$B)=(split /[\_\.]/,$file)[0,2];
	my $cmd=qq(awk \'\$1!=\"Tag-Seq\"{print \$2"\t"\$3}\' $file |sort -n -k 1> relativity.dat);
	system "$cmd";
	my $fh="$file.gnuplot.dat";
open(PL, ">$fh") || die $!;
print PL qq(
reset;
set title '$A and $B expression enrichment relativity dd plot';
set xlabel '$A';
set ylabel '$B';
set nokey;
plot "relativity.dat" u 1:2 w dots;
set terminal svg;
set output "$A\_$B.ExpRelativity.svg";
replot "relativity.dat" u 2:1 w dots;
);
close (PL);
system "gnuplot $fh";
}
