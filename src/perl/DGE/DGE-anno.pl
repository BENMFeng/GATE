#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
my $program_name = $1 if ( $0 =~ /([^\/]+)$/);
my $usage = <<USAGE;
$program_name -- Gene Expression Annotation

Author: BENM <binxiaofeng\@gmail.com>
Version: 0.1 alpha
Date: 2010-07-14
Update: 0.2 alpha 2010-12-30	Add \"gene_format\";

$program_name [options] [functions]
	--len <int>		set tags reads length, default: 17
	--i <file>		input file, tag reads should be located in first column
	--ref <file>		input reference sequence file
	--gene_format <str>	set input reference gene format [oases]
	--aln <file>		input alignment file
	--a2 <file>		input alignment file of which tags were trimmed CATG, uniq mapping
	--f <str>		set alignment format, default: MAQ
	--enzyme_site <str>	set enzyme site, default: CATG
	--unanno		output unannotation tags
	--help			help & usuage

Example:

perl $program_name --len 17 --i CopyNumber.xls --ref ref.fa --aln aln.mapview --f MAQ --enzyme_site CATG > Anno-CopyNumber.xls

USAGE

my ($Len,$Input,$Ref,$Gene_format,$Aln,$A2,$Format,$Site,$Unanno,$Help);
my %opts;
GetOptions
(
	\%opts,
	"len:i"=>\$Len,
	"i:s"=>\$Input,
	"ref:s"=>\$Ref,
	"gene_format:s"=>\$Gene_format,
	"aln:s"=>\$Aln,
	"a2:s"=>\$A2,
	"f:s"=>\$Format,
	"enzyme_site"=>\$Site,
	"unanno"=>\$Unanno,
	"help"=>\$Help
);
die $usage if ((!defined $Input)||(!defined $Ref)||((!defined $Aln)&&(!defined $A2))||($Help));


$Format ||= "MAQ";
$Site ||= "CATG";
$Len =($Site=~/CATG/) ? 17 : 16;

#Alignment	MAQ	SOAP2	BOWTIE	BWA
#chromosome	1	7	2 	2
#position	2	8	3 	3
#strand		3	6	1 	(1)
#sequence	14	1	4 	9

my ($A,$B,$C,$D);
if ($Format =~ /MAQ/i)
{
	($A,$B,$C,$D)=(1,2,3,14);
}
elsif ($Format =~ /SOAP/i)
{
	($A,$B,$C,$D)=(7,8,6,1);
}
elsif ($Format =~ /BOWTIE/i)
{
	($A,$B,$C)=(2,3,1,4);
}

my %Seq=();
read_aln($Aln,\%Seq) if (defined $Aln);
my %Tag=();
read_aln($A2,\%Tag) if (defined $A2);

my ($chr,$seq)=("","");
open (REF,$Ref) || die $!;
my %hash;
while(<REF>)
{
	if (/^\>(\S+)/)
	{
		if ((length($seq)>0)&&((exists $Seq{$chr})||(exists $Tag{$chr})))
		{
			if (exists $Seq{$chr})
			{
				foreach my $reads(keys %{$Seq{$chr}})
				{
					my @hit=split /\;/,$Seq{$chr}{$reads};
					for (my $i=0;$i<@hit;$i++)
					{
						next if ($hit[$i] !~ /\t/);
						my ($p,$s)=split /\t/,$hit[$i];
						$hash{$reads}.="$chr, "."$p, $s, ".anno($seq,$s,$p)."; ";
					}
					delete $Seq{$chr}{$reads};
				}
			}
			if (exists $Tag{$chr})
			{
				foreach my $reads(keys %{$Tag{$chr}})
				{
					if (!exists $Seq{$chr}{$reads})
					{
						my @hit=split /\;/,$Tag{$chr}{$reads};
						for (my $i=0;$i<@hit;$i++)
						{
							next if ($hit[$i] !~ /\t/);
							my ($p,$s)=split /\t/,$hit[$i];
							$hash{$reads}.="$chr, "."$p, $s, ".anno($seq,$s,$p)."; ";
						}
					}
					delete $Tag{$chr}{$reads};
				}
			}
		}
		$chr=$1;
		$seq="";
	}
	else
	{
		s/\s*//g;
		$seq.=$_;
	}
}
if ((length($seq)>0)&&(exists $Seq{$chr}))
{
	foreach my $reads(keys %{$Seq{$chr}})
	{
		my @hit=split /\;/,$Seq{$chr}{$reads};
		for (my $i=0;$i<@hit;$i++)
		{
			next if ($hit[$i] !~ /\t/);
			my ($p,$s)=split /\t/,$hit[$i];
			$hash{$reads}="$chr, "."$p, $s, ".anno($seq,$s,$p)."; ";
		}
		delete $Seq{$chr}{$reads};
	}
}
close REF;

open (IN,$Input) || die $!;
$_=<IN>;
chomp;
print "$_\tMark\tMatchGene, Location, Strand, TagsReads, Nth-$Site-from-3\'-end-of-Gene, Distance-Tag-from-3\'-end-of-Gene; ...\n";
while(<IN>)
{
	chomp;
	my @t=split;
	$t[0]=(length($t[0])>$Len) ? substr($t[0],length($t[0])-$Len,$Len) : $t[0];
	if (exists $hash{$t[0]})
	{
		my @semicolon=split /\;\s+/,$hash{$t[0]};
		my $q;
		if (@semicolon>1)
		{
			if (defined $Gene_format)
			{
				if ($Gene_format eq "oases")
				{
					my $j=0;
					my $str=$1 if ($semicolon[0]=~/Locus\_(\d+)/);
					if ($semicolon[0]=~/\s+\+\,\s+$Site/)
					{
						for (my $i=0;$i<@semicolon;$i++)
						{
							if ($semicolon[$i]=~/Locus\_(\d+)/)
							{
								if (($str ne $1)||($semicolon[$i]!~/\s+\+\,\s+$Site/))
								{
									$j++;
								}
							}
						}
						$q=($j==0)?"+":"*";
					}
					elsif ($semicolon[0]=~/\s+\-\,\s+$Site/)
					{
						for (my $i=0;$i<@semicolon;$i++)
						{
							if ($semicolon[$i]=~/Locus\_(\d+)/)
							{
								if (($str ne $1)||($semicolon[$i]!~/\s+\-\,\s+$Site/))
								{
									$j++;
								}
							}
						}
						$q=($j==0)?"-":"*";
					}
					else
					{
						$q="*";
					}
				}
				else
				{
					$q="*";
				}
			}
			else
			{
				$q="*";
			}
		}
		else
		{
			if ($semicolon[0]=~/\s+\+\,\s+$Site/)
			{
				$q="+";
			}
			elsif ($semicolon[0]=~/\s+\-\,\s+$Site/)
			{
				$q="-";
			}
			else
			{
				$q="?";
			}
		}
		print (join "\t",@t,"$q\t$hash{$t[0]}\n");
	}
	else
	{
		print (join "\t",@t,"null\t\-\n") if ($Unanno);
	}
}
close IN;

sub read_aln
{
	my ($file,$hash)=@_;
	open (ALN,$file) || die $!;
	while(<ALN>)
	{
		next if (/^\#/);
		my @t=split;
		my ($seq,$pos)=(length($t[$D])>$Len)?(substr($t[$D],length($t[$D])-$Len,$Len),$t[$B]-1+length($t[$D])-$Len):($t[$D],$t[$B]-1);
		$$hash{$t[$A]}{$seq}.="$pos\t$t[$C];"; #position	strand
	}
	close ALN;
}

sub anno
{
	my ($seq,$s,$p)=@_;
	my $tag="";
	my $subseq="";
	my $n=0;
	my $loc_3=0;
	if ($s eq "+")
	{
		$loc_3=length($seq)-$p-$Len;
		$tag=substr($seq,$p-length($Site),length($Site)+$Len);
		$subseq=substr($seq,$p+$Len,$loc_3);
		while ($subseq=~/$Site/g)
		{
			$n++;
		}
	}
	elsif ($s eq "-")
	{
		$loc_3=$p;
		$tag=reverse(substr($seq,$p,length($Site)+$Len));
		$tag=~tr/ACGTacgt/TGCAtgca/;
		$subseq=reverse(substr($seq,0,$loc_3));
		$subseq=~tr/ACGTacgt/TGCAtgca/;
		while ($subseq=~/$Site/g)
		{
			$n++;
		}
	}
	$tag=($tag =~ /^$Site/) ? $tag : substr($tag,length($tag)-$Len,$Len);
	return "$tag, $n, $loc_3";
}
