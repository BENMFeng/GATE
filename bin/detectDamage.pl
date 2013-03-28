#!/usr/bin/perl -w
## Detect DNA damage from aDNA sequencing alignment and varation calling.
## Author: BENM <binxiaofeng@gmail.com>
## v1.3, 2013-03-23
use strict;
use warnings;
use Getopt::Long;
use constant G=>0.6180339887498948482045868343656;
my ($mapQ,$BAQ,$DP,$QUAL,$Haplo,$SnV,$Help);
my %opts;
GetOptions(
	\%opts,
	"mapQ:i"=>\$mapQ,
	"BAQ:i"=>\$BAQ,
	"DP:i"=>\$DP,
	"haplo"=>\$Haplo,
	"correctSNV:s"=>\$SnV,
	"help"=>\$Help
);

die "perl $0 <IN:gatk.vcf> <IN:realigned.baq.bam> <IN:ref.fa> <chr> [-correctSNV snp.xls]\n" if (@ARGV<4);
$mapQ ||= 0;
$BAQ ||= 13;
$DP ||= 3;

my ($Vcf,$Bam,$Ref,$Chr)=@ARGV;
my %Ref;
my %Mutation;
my %HET;
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  M13-62
#MT	73	.	A	G	735	.	AC=2;AF=1.00;AN=2;DP=35;Dels=0.00;FS=0.000;HaplotypeScore=0.9469;MLEAC=2;MLEAF=1.00;MQ=34.27;MQ0=0;QD=21.00;SB=-3.640e+02	GT:AD:DP:GQ:PL	1/1:0,33:35:60:768,60,0
#HaplotypeScore
#Consistency of the site with two (and only two) segregating haplotypes.
#Higher scores are indicative of regions with bad alignments, often leading to artifactual SNP and indel calls. Note that the Haplotype Score is only calculated for sites with read coverage.
open (VCF,$Vcf) if ($Vcf=~/vcf$/);
open (VCF,"tabix $Vcf $Chr|") if ($Vcf=~/vcf\.gz$/);
while(<VCF>)
{
	next if (/^\#/ || /INDEL/);
	my @col=split /\t+/,$_;
	next if ($col[0] ne $Chr || $col[4] !~ /[ACGT]/);
	$Ref{$col[1]}=$col[3];
	$Mutation{$col[1]}=$col[4];
	if ($col[9]!~/1\/1/)
	{
		my $GQ=(split /\:/,$col[9])[-2];
		if ($col[8]=~/HaplotypeScore\=([^\;\s]+)\;/)
		{
			next if ($1<1 || $GQ<30);
			next if (($col[3] eq "C" && $col[4] eq "T") || ($col[3] eq "G" && $col[4] eq "A") );
			$HET{$col[1]}=$col[4] if ($1>1);
		}
	}
}
close VCF;

my %hash;
my @prestrand=();
open (IN,"samtools view -u $Bam $Chr | samtools mpileup -q $mapQ -Q $BAQ -f $Ref -O - |") || die $!;
#MT      1736    A       31      ,.G.GG....G.,.,,G..,G...,..gG,^],       HJ:J64HJEJ3GHJFF5EJE4BIHBH@98;> 88,58,56,55,54,51,49,44,41,39,38,54,36,32,30,29,44,25,22,21,18,16,14,10,10,9,22,5,5,3,1
while(<IN>)
{
	chomp;
	my ($chr,$coor,$refbase,$depth,$cns,$qual,$pos)=split /\t+/,$_;
	next if (!exists $Mutation{$coor});
	$cns=~s/[\+\-]\d+[ACGTacgt]+//g;
	$cns=~s/[^\,\.ACGTacgt]+//g;
	if ($depth<$DP || $cns!~/[ACGTacgt]/)
	{
		push @{$hash{$coor}},($refbase,$Mutation{$coor},"$refbase\>$Mutation{$coor}",0,0,0.0000,"-");
		next;
	}
	my %ntStrandFreq;
	my %ntFreq;
	my %ntStrandPos;
	my %ntPos=();
	my @allpos=split /\,/,$pos;
	while($cns=~/([\,\.ACGTacgt])/g)
	{
		my ($b,$poisonreads)=($1,$-[0]);
		$b=$refbase if ($b=~/[\,\.]/);
		if (defined $prestrand[$poisonreads])
		{
			$ntStrandFreq{uc($b)}{$prestrand[$poisonreads]}++;
			$ntStrandPos{uc($b)}{$prestrand[$poisonreads]}{$allpos[$poisonreads]}++;
		}
		else
		{
			$ntStrandFreq{uc($b)}{'.'}++;
			$ntStrandPos{uc($b)}{'.'}{$allpos[$poisonreads]}++;
		}
	}
	foreach my $nt (keys %ntStrandFreq)
	{
		if (exists $ntStrandFreq{$nt}{'.'} && exists $ntStrandFreq{$nt}{','})
		{
			$ntFreq{$nt}=($ntStrandFreq{$nt}{'.'} > $ntStrandFreq{$nt}{','})?$ntStrandFreq{$nt}{'.'}:$ntStrandFreq{$nt}{','};
			%{$ntPos{$nt}}=(keys %{$ntStrandPos{$nt}{'.'}} > keys %{$ntStrandPos{$nt}{','}})?%{$ntStrandPos{$nt}{'.'}}:%{$ntStrandPos{$nt}{','}};
		}
		elsif (exists $ntStrandFreq{$nt}{'.'})
		{
			$ntFreq{$nt}=$ntStrandFreq{$nt}{'.'};
			%{$ntPos{$nt}}=%{$ntStrandPos{$nt}{'.'}};
		}
		elsif (exists $ntStrandFreq{$nt}{','})
		{
			$ntFreq{$nt}=$ntStrandFreq{$nt}{','};
			%{$ntPos{$nt}}=%{$ntStrandPos{$nt}{','}};
		}
	}
	@prestrand=();
	while($cns=~/([\,\.ACGTacgt])/g)
	{
		my $s=$1;
		my $tail=$';
		if ($s !~ /[\,\.]/)
		{
			if (defined $prestrand[-1])
			{
				$s = $prestrand[-1];
			}
			elsif ($tail=~/^([\,\.])/)
			{
				$s = $1;
			}
			else
			{
				$s = '.';
			}
		}
		push @prestrand,$s;
	}
	if (($ntFreq{$Mutation{$coor}}/length($cns) >= 0.8) && ( (keys %{$ntPos{$Mutation{$coor}}} >= 2) || ($Mutation{$coor}!~/[AT]/i) ) )
	{
		$Mutation{$coor}.="|$depth|".sprintf("%.4f",$ntFreq{$Mutation{$coor}}/length($cns));
		next;
	}
	my @ntbase=sort{$ntFreq{$b}<=>$ntFreq{$a}}keys %ntFreq;
	if ( @ntbase >= 2 && ($ntFreq{$ntbase[1]}>=2 || $ntFreq{$ntbase[1]}/length($cns) >= 0.3) && $ntFreq{$ntbase[1]}/length($cns) < 0.8 )
	{
		if ( (@ntbase > 2 && $ntFreq{$ntbase[2]}>=2) || (keys %{$ntPos{$ntbase[0]}} <= 3) || 
		   (keys %{$ntPos{$ntbase[0]}} > 3 && keys %{$ntPos{$ntbase[1]}} >= 3) )
		{
			my $raw=$ntbase[0];
			my $j=1;
			for (my $i=0;$i<@ntbase;$i++)
			{
				if ($ntbase[$i] =~ /([AT])/i)
				{
					if ($ntbase[$i]=~/A/ && exists $ntFreq{'G'} && ($ntFreq{'G'}>=2 || $ntFreq{'G'}/length($cns)>=1/3))
					{
						$raw="G";
					}
					elsif ($ntbase[$i]=~/T/ && exists $ntFreq{'C'} && ($ntFreq{'C'}>=2 || $ntFreq{'C'}/length($cns)>=1/3))
					{
						$raw="C";
					}
					$j=$i;
				}
			}
			push @{$hash{$coor}},($refbase,$Mutation{$coor},"$raw\>$ntbase[$j]");
			#print "$chr\t$coor\t$refbase\t$Mutation{$coor}\t$raw\>$ntbase[$j]\t";
			if ( (($raw eq "C" && $ntbase[$j] eq "T") || ($raw eq "G" && $ntbase[$j] eq "A")) && (keys %{$ntPos{$ntbase[$j]}} <= 3) )
			{
				#print "TRUE";
				push @{$hash{$coor}},1;
			}
			else
			{
				#print "FALSE";
				push @{$hash{$coor}},0;
			}
			#print "\t$ntFreq{$ntbase[$j]}\t";
			push @{$hash{$coor}},$ntFreq{$ntbase[$j]};
			push @{$hash{$coor}},sprintf ("%.4f",$ntFreq{$ntbase[$j]}/length($cns));
			#printf ("%.4f",$ntFreq{$ntbase[$j]}/length($cns));
			#print "\t";
			my @posary=();
			foreach my $p(sort{$a<=>$b}keys %{$ntPos{$ntbase[$j]}})
			{
				foreach (1..$ntPos{$ntbase[$j]}{$p})
				{
					push @posary,$p;
				}
			}
			#print ((join ",",@posary),"\n");
			push @{$hash{$coor}},(join ",",@posary);
		}
	}
}
close IN;
print "#CHR\tPOS\tREF\tSNV\tChangedBase\tSourceFromDamage\tScore\tFrequency\tPresentRate\tPos-in-reads\n";
my ($last,$next);
my @coor=sort{$a<=>$b}keys %hash;
my $Score=1;
for (my $p=0;$p<@coor;$p++)
{
	my $tmp=$coor[$p];
	if ( (defined $last) && ($tmp-$last<=100) )
	{
		if ( (${$hash{$tmp}}[3] | ${$hash{$last}}[3])==0)
		{
			$Score-=${$hash{$tmp}}[5]+${$hash{$last}}[5];
		}
		else
		{
			$Score-=${$hash{$last}}[5] if (${$hash{$last}}[3]==0);
		}
	}
	$next=$coor[$p+1];
	if (($p+1<@coor) && ($next-$tmp<=100) )
	{
		if ( (${$hash{$tmp}}[3] | ${$hash{$next}}[3])==0)
		{
			$Score-=${$hash{$tmp}}[5]+${$hash{$next}}[5];
		}
		else
		{
			$Score-=${$hash{$next}}[5] if (${$hash{$next}}[3]==0);
		}
	}
	if ( (!defined $last || $tmp-$last>100) && ($p+1>=@coor || $next-$tmp>100) && (${$hash{$tmp}}[3]==0))
	{
		$Score-=${$hash{$tmp}}[5];
	}
	$Score=0 if ($Score<0);
	$last=$tmp;
	print "$Chr\t$tmp\t";
	print (join "\t",@{$hash{$tmp}}[0..2]);
	if (${$hash{$tmp}}[3]==1 && $Score==1)
	{
		print "\tTRUE\t$Score\t";
		my ($raw,$change)=split /\>/,${$hash{$tmp}}[2];
		if ($change eq ${$hash{$tmp}}[1] && $raw eq ${$hash{$tmp}}[0])
		{
			$Mutation{$tmp}.="|TRUE|".sprintf("%.4f",${$hash{$tmp}}[4]/${$hash{$tmp}}[5]);
			if (($Ref{$tmp} eq "C" && $Mutation{$tmp} eq "T") || ($Ref{$tmp} eq "G" && $Mutation{$tmp} eq "A") )
			{
				undef $Mutation{$tmp};
				delete $Mutation{$tmp};
			}
		}
	}
	elsif (${$hash{$tmp}}[3]==0 && $Score==0 && exists $HET{$tmp})
	{
		print "\tHET\t$Score\t";
		$Mutation{$tmp}.="|HET|".sprintf("%.4f",${$hash{$tmp}}[4]/${$hash{$tmp}}[5]);
	}
	else
	{
		print "\tFALSE\t$Score\t";
		if (${$hash{$tmp}}[5]>0)
		{
			$Mutation{$tmp}.="|FALSE|".sprintf("%.4f",${$hash{$tmp}}[4]/${$hash{$tmp}}[5]);
			undef $Mutation{$tmp};
			delete $Mutation{$tmp};
		}
		else
		{
			$Mutation{$tmp}.="|FALSE|0.0000";
			undef $Mutation{$tmp};
			delete $Mutation{$tmp};
		}
	}
	print ((join "\t",@{$hash{$tmp}}[4..6]),"\n");
	$Score=1;
}

if (defined $SnV)
{
	open (OUT,">$SnV");
	foreach my $p(sort{$a<=>$b}keys %Mutation)
	{
		next if (!exists $Mutation{$p} || !defined $Mutation{$p});
		if ($Mutation{$p}=~/([^\|\s]+)\|(\S+)/)
		{
			print OUT "$Chr\t$p\t$Ref{$p}\t$1\t$2\n";
		}
		else
		{
			print OUT "$Chr\t$p\t$Ref{$p}\t$Mutation{$p}\n";
		}
	}
	close OUT;
}