#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my %opts;
my ($Reads1,$Reads2,$Index,$Qual,$Prefix,$Help);
GetOptions(%opts,"fastq1:s"=>\$Reads1,"fastq2:s"=>\$Reads2,"fastq:s"=>\$Reads1,"index:s"=>\$Index,"qual:s"=>\$Qual,"prefix:s"=>\$Prefix,"help"=>\$Help);

die qq(perl $0 -fastq1 reads1.fastq -fastq2 reads2.fastq -index idx1:idx2 [-prefix prefixName] [-qual 30 or "?"]\n) if (!defined $Reads1 || !defined $Index || defined $Help);

$Qual ||= 30;
$Qual = ord($Qual)-33 unless ($Qual=~/^\d+/ && $Qual>=10);
my ($Idx1,$Idx2)=split /\:/,$Index;
my ($Out1,$Out2)=("","");
if (defined $Prefix)
{
	$Out1="$Prefix.$Index\_1.fastq";
	$Out2="$Prefix.$Index\_2.fastq"
}
else
{
	$Out1=(split /\//,$Reads1)[-1];
	if ($Out1=~/([^\/\s]+)\.fastq$/ || $Out1=~/([^\/\s]+)\.fastq\.gz$/)
	{
		$Out1="$1\.$Idx1.fastq";
	}
	if (defined $Reads2 && $Reads2 ne "")
	{
		$Out2=(split /\//,$Reads2)[-1];
		my $Idx=(defined $Idx2 && $Idx2 ne "")?$Idx2:$Idx1;
		if ($Out2=~/([^\/\s]+)\.fastq$/ || $Out2=~/([^\/\s]+)\.fastq\.gz$/)
		{
			$Out2="$1.$Idx.fastq";
		}
	}
}


open (IN1,$Reads1) if ($Reads1=~/fastq$/i || $Reads1=~/fq$/i);
open (IN1,"gzip -cd $Reads1|") if ($Reads1=~/fastq\.gz$/i || $Reads1=~/fq\.gz$/i);
open (IN2,$Reads2) if ((defined $Reads2 && $Reads2 ne "") && ($Reads2=~/fastq$/i || $Reads2=~/fq$/i) );
open (IN2,"gzip -cd $Reads2|") if ((defined $Reads2 && $Reads2 ne "") && ($Reads2=~/fastq\.gz$/i || $Reads2=~/fq\.gz$/i) );

open (OUT1,">$Out1") || die $!;
open (OUT2,">$Out2") if (defined $Reads2 && $Reads2 ne "");
while(<IN1>)
{
	my $out1=$_;
	$_=<IN1>;
	my $withIdx=0;
	if ($_=~/^$Idx1(\S+)/) {
		$out1.="$1\n";
		$withIdx=1;
	}
	$_=<IN1>;
	$out1.=$_;
	$_=<IN1>;
	chomp;
	my $qual1=substr($_,0,length($Idx1));
	$out1.=($withIdx==1)?substr($_,length($Idx1),length($_)-length($Idx1))."\n":"$_\n";
	my $qualcheck=check_qual($qual1);
	if (defined $Reads2 && $Reads2 ne "")
	{
		$_=<IN2>;
		my $out2=$_;
		$_=<IN2>;
		if (defined $Idx2 && $Idx2 ne "")
		{
			if ($_=~/^$Idx2(\S+)/)
			{
				$out2.="$1\n";
			}
			else
			{
				$withIdx=0;
			}
		}
		else
		{
			$out2.=$_;
		}
		$_=<IN2>;
		$out2.=$_;
		my $qual2=substr($_,0,length($Idx2)) if (defined $Idx2);
		$qualcheck=check_qual($qual2) if (defined $qual2 && $qual2 ne "");
		$_=<IN2>;
		chomp;
		$out2.=($withIdx==1 && defined $Idx2 && $Idx2 ne "")?substr($_,length($Idx2),length($_)-length($Idx2))."\n":"$_\n";
		print OUT2 $out2 if ($withIdx == 1 && $qualcheck==1);
	}
	print OUT1 $out1 if ($withIdx == 1 && $qualcheck==1);
}
close IN1;
close IN2;
close OUT1;
close OUT2;

sub check_qual{
	my $qual=shift;
	my @q=split //,$qual;
	for (my $i=0;$i<@q;$i++)
	{
		return 0 if (ord($q[$i])-33<$Qual);
	}
	return 1;
}

__END__