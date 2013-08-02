#!/usr/bin/perl -w
## Author: BENM <binxiaofeng@gmail.com>
## v1.0, 2013-01-14
use strict;
use Data::Dumper;
use Getopt::Long;

my ($bam,$chr,$ref,$prefix,$vcf)=@ARGV;
my ($mapQ,$BAQ,$DP,$QUAL,$haplo,$Help);
my %opts;
GetOptions(
	\%opts,
	"mapQ:i"=>\$mapQ,
	"BAQ:i"=>\$BAQ,
	"DP:i"=>\$DP,
	"QUAL:i"=>\$QUAL,
	"haplo"=>\$haplo,
	"help"=>\$Help
);

die "perl $0 <aln.bam> <chromosome_id> <ref.fa> <out_prefix> [gatk.vcf] [-haplo] [-mapQ 0] [-BAQ 13] [-DP 3] [-QUAL 50]\n" if (@ARGV<4 || $Help);

my %abbrev=(
	"M"=>"A/C",
	"R"=>"A/G",
	"W"=>"A/T",
	"S"=>"C/G",
	"Y"=>"C/T",
	"K"=>"G/T",
);

my %HET=();
while(my ($key,$value)=each %abbrev)
{
	$HET{$value}=$key;
}

$mapQ ||= 0;
$BAQ ||= 13;
$DP ||= 3;
$QUAL ||= 50;

my $fa="$prefix.$chr.var.fa";

if (!defined $vcf || !-f $vcf)
{
	$vcf="$prefix.$chr.vcf";
	system "samtools view -u $bam $chr | samtools mpileup -q $mapQ -Q $BAQ -uf $ref - | bcftools view -cg - > $vcf";
	system qq(vcfutils.pl vcf2fq $vcf | perl -ne 'if (/^\\@\\S+/){print ">$prefix\_MT\\n"}elsif (/^\\+/){exit}else{print}' > $fa);
}
else
{
	system qq(samtools view -u $bam $chr | samtools mpileup -q $mapQ -Q $BAQ -uf $ref - | bcftools view -cg - | vcfutils.pl vcf2fq  - | perl -ne 'if (/^\\@\\S+/){print ">$prefix\_MT\\n"}elsif (/^\\+/){exit}else{print}' > $fa); 
}

my %cns;
my @polymorphisms=();
open (IN,$vcf) if ($vcf=~/\.vcf$/);
open (IN,"tabix $vcf $chr|") if ($vcf=~/\.vcf.gz$/);
while(<IN>)
{
	next if (/^\#/ || /INDEL/);
	my @col=split /\t+/,$_;
	if ($col[0] eq $chr)
	{
		if ($haplo && $col[4]=~/[ACGT]/ && $col[4] ne $col[3])
		{
			my $depth=$1 if ($col[7]=~/DP\=(\d+)\;/);
			my @base=($col[3],(split /\,/,$col[4]));
			if ($depth>=$DP)
			{
				my %allel=();
				for (my $i=9;$i<@col;$i++)
				{
					my $PL;
					if ($col[8] eq "GT:PL:GQ")
					{
						$PL=(split /\:/,$col[$i])[1];
					}
					elsif ($col[8] eq "GT:AD:DP:GQ:PL")
					{
						$PL=(split /\:/,$col[$i])[-1];
					}
					if ($PL !~ /0\,0\,0/)
					{
						my @PLary=split /\,/,$PL;
						my $n=0;
						for (my $i=0;$i<@base;$i++)
						{
							for (my $j=0;$j<=$i;$j++)
							{
								if ($base[$i] eq $base[$j])
								{
									if (exists $allel{$base[$i]})
									{
										$allel{$base[$i]}*=10**(-$PLary[$n]/10);
									}
									else
									{
										$allel{$base[$i]}=10**(-$PLary[$n]/10);
									}
								}
								else
								{
									my $het=join "/",sort($base[$i],$base[$j]);
									if (exists $HET{$het})
									{
										if (exists $allel{$HET{$het}})
										{
											$allel{$HET{$het}}*=10**(-$PLary[$n]/10);
										}
										else
										{
											$allel{$HET{$het}}=10**(-$PLary[$n]/10);
										}
									}
									else
									{
										print STDERR "$col[1]\t$het\n";
										if (exists $allel{$het})
										{
											$allel{$het}*=10**(-$PLary[$n]/10);
										}
										else
										{
											$allel{$het}=10**(-$PLary[$n]/10);
										}
									}
								}
								$n++;
							}
						}
					}
				}
				my $alt=$col[4];
				if (keys %allel>=2)
				{
					my ($fst,$sec)=(sort{$allel{$b}<=>$allel{$a}}keys %allel)[0,1];
					if ($allel{$fst}==$allel{$sec})
					{
						if ($fst ne $col[3])
						{
							$alt=$fst;
						}
						if ($sec ne $col[3])
						{
							$alt=$sec;
						}
					}
					else
					{
						$alt=$fst;
					}
				}
				elsif (keys %allel>0)
				{
					$alt=(sort{$allel{$b}<=>$allel{$a}}keys %allel)[0];
				}
				else
				{
					$alt=($col[5]>=$QUAL)?$col[4]:$col[3];
				}
				if (defined $haplo && exists $abbrev{$alt})
				{
					if (($col[5]<$QUAL)&&(($col[3] eq "C" && $col[4] eq "T") || ($col[3] eq "G" && $col[4] eq "A")))
					{
						$alt=$col[3];
					}
					elsif ($col[4] !~ /\,/)
					{
						$alt=$col[4];
					}
				}
				$cns{$col[1]}=$alt;
			}
			else
			{
				$cns{$col[1]}=($col[5]>=$QUAL)?$col[4]:$col[3];
			}
			push @polymorphisms,"$col[1]$cns{$col[1]}" if ($cns{$col[1]} ne $col[3] && (defined $haplo && !exists $abbrev{$cns{$col[1]}}));
		}
		else
		{
			$cns{$col[1]}=$col[3];
		}
	}
}
close IN;

open (IN,$fa) || die $!;
open (OUT,">$prefix.MT.cns.fa") || die $!;
my $loc=1;
while(<IN>)
{
	if (/\>/)
	{
		print OUT ">$prefix\_MT\n";
	}
	else
	{
		if ($_=~/[^ACGTN]/i)
		{
			s/\s*$//;
			my @base=split '',$_;
			for(my $i=0;$i<@base;$i++)
			{
				my $j=$loc+$i;
				if (exists $abbrev{uc($base[$i])})
				{
					my @het=split /\//,$abbrev{uc($base[$i])};
					foreach my $b(@het)
					{
						next if ($b eq uc($cns{$j}));
						$base[$i]=$b;
					}
				}
			}
			print OUT ((join '',@base),"\n");
		}
		else
		{
			print OUT $_;
		}
		$loc+=length($_);
	}
}
close IN;
close OUT;

if (defined $haplo)
{
	my $seq;
	open (IN,"$prefix.MT.cns.fa") || die $!;
	while(<IN>)
	{
		if (!/^\>/)
		{
			s/\s*//g;
			$seq.=$_;
		}
	}
	close IN;
	my $region="";
	while($seq=~/([^Nn]+)/gi)
	{
		$region.=($-[0]+1)."-".($-[0]+length($1)).";";
	}
	print qq($prefix\t"$region"\t\?\t);
	print (join "\t",@polymorphisms);
	print "\n";
}