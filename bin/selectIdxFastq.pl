#!/usr/bin/perl -w
# Authore: BENM <binxiaofeng\@gmail.com>
# Version: v1.1, 2013-05-03
use strict;
use Getopt::Long;
use Data::Dumper;

my %opts;
my ($Reads1,$Reads2,$Barcode,$Index,$Qual,$Mismatch,$Prefix,$Debug,$Stat,$Help);
GetOptions(%opts,"fastq1:s"=>\$Reads1,"fastq2:s"=>\$Reads2,"fastq:s"=>\$Reads1,"index:s"=>\$Index,"barcode:s"=>\$Barcode,"mis:i"=>\$Mismatch,"qual:s"=>\$Qual,"prefix:s"=>\$Prefix,"debug"=>\$Debug,"stat"=>\$Stat,"help"=>\$Help);

die qq(perl $0 -fastq1 reads1.fastq -fastq2 reads2.fastq [-index index] [-barcode bc1:bc2] [-prefix prefixName] [-mis 1] [-qual 30 or "?"]\n) if ((!defined $Reads1 || (!defined $Barcode && !defined $Index) || defined $Help)&&(!defined $Debug));

$Qual ||= 30;
$Qual = ord($Qual)-33 unless ($Qual=~/^\d+/ && $Qual>=10);
$Mismatch ||= 1 if (!defined $Mismatch);
my ($Bar1,$Bar2)=split /\:/,$Barcode if (defined $Barcode);
my ($Out1,$Out2,$StatOut)=("","","");
my %Idx=();
my $Index_len=(defined $Index)?length($Index):0;
if (defined $Index && $Mismatch>0) {
	my @idx_ary=split //,$Index;
	my @idx_idx=();
	for (my $i=0;$i<$Index_len;$i++) {
		$idx_idx[$i]=$i;
	}
	my @COMB=combinatorialSelect($Mismatch,\@idx_idx);
	$Idx{$Index}=1;
	my @base=("A","C","G","T","N");
	for (my $j=0;$j<@COMB;$j++) {;
		my @var_loc=@{$COMB[$j]};
		#print ((join "",@var_loc),"\n") if ($Debug);
		for (my $k=0;$k<5**scalar(@var_loc);$k++)
		{
			my @tmp=@idx_ary;
			my $m = $k % 5;
			my $left=$k/5;
			while ($left>=1) {
				$m .= $left % 5;
				$left = $left/5;
			}
			if (length($m)<$Mismatch) {
				$m .= "0"x($Mismatch-length($m));
			}
			my @loc_idx=split "",$m;
			#print ((join " ",@loc_idx),"\n") if ($Debug);
			for (my $i=0;$i<@loc_idx;$i++)
			{
				$tmp[$var_loc[$i]]=$base[$loc_idx[$i]];
			}
			my $idx_p=join "",@tmp;
			print "$idx_p\n" if ($Debug);
			$Idx{$idx_p}=1;
		}
	}
} elsif ($Mismatch==0) {
	$Idx{$Index}=1;
}
#print Dumper %Idx if ($Debug);
exit() if ($Debug);

if (defined $Prefix)
{
	$Out1=$Prefix;
	$StatOut="$Prefix.idx.stat";
	$Out1.="\.$Index" if (defined $Index && $Index ne "");
	$Out1.="\.$Bar1" if (defined $Bar1 && $Bar1 ne "");
	$Out1.="\_1.fastq";
	$Out2=$Prefix;
	my $bar=(defined $Bar2 && $Bar2 ne "")?$Bar2:$Bar1;
	$Out2.="\.$Index" if (defined $Index && $Index ne "");
	$Out2.="\.$bar" if (defined $bar && $bar ne "");
	$Out2.="\_2.fastq";
}
else
{
	$Out1=(split /\//,$Reads1)[-1];
	if ($Out1=~/([^\/\s]+)\.fastq$/ || $Out1=~/([^\/\s]+)\.fastq\.gz$/)
	{
		$Out1= $1;
		$StatOut="$Out1.idx.stat";
		$Out1.="\.$Index" if (defined $Index && $Index ne "");
		$Out1.="\.$Bar1" if (defined $Bar1 && $Bar1 ne "");
		$Out1.="\.fastq";
	}
	if (defined $Reads2 && $Reads2 ne "")
	{
		$Out2=(split /\//,$Reads2)[-1];
		my $bar=(defined $Bar2 && $Bar2 ne "")?$Bar2:$Bar1;
		if ($Out2=~/([^\/\s]+)\.fastq$/ || $Out2=~/([^\/\s]+)\.fastq\.gz$/)
		{
			$Out2=$1;
			$Out2.="\.$Index" if (defined $Index && $Index ne "");
			$Out2.="\.$bar" if (defined $bar && $bar ne "");
			$Out2.="\.fastq";
		}
	}
}


open (IN1,$Reads1) if ($Reads1=~/fastq$/i || $Reads1=~/fq$/i);
open (IN1,"gzip -cd $Reads1|") if ($Reads1=~/fastq\.gz$/i || $Reads1=~/fq\.gz$/i);
open (IN2,$Reads2) if ((defined $Reads2 && $Reads2 ne "") && ($Reads2=~/fastq$/i || $Reads2=~/fq$/i) );
open (IN2,"gzip -cd $Reads2|") if ((defined $Reads2 && $Reads2 ne "") && ($Reads2=~/fastq\.gz$/i || $Reads2=~/fq\.gz$/i) );

my %StatHash=();
open (OUT1,">$Out1") || die $!;
open (OUT2,">$Out2") if (defined $Reads2 && $Reads2 ne "");
while(<IN1>)
{
	my $out1=$_;
	my $withBar=(defined $Barcode)?0:1;
	my $withIdx=(defined $Index)?0:1;
	if (/\:([ACGTN]{$Index_len})/) {
		my $reads_idx=$1;
		$StatHash{'Idx'}{$reads_idx}++ if (defined $Stat);
		$withIdx=1 if (exists $Idx{$reads_idx});
	}
	$_=<IN1>;
	if (defined $Bar1 && $_=~/^$Bar1(\S+)/) {
		$out1.="$1\n";
		$StatHash{'bar'}{$Bar1}++ if (defined $Stat);
		$withBar=1;
	} elsif (!defined $Bar1) {
		$out1.=$_;
	}
	$_=<IN1>;
	$out1.=$_;
	$_=<IN1>;
	chomp;
	my $qual1=substr($_,0,length($Bar1)) if (defined $Bar1);
	my $qualcheck=(defined $Barcode)?check_qual($qual1):1;
	$out1.=($withBar==1 && defined $Bar1)?substr($_,length($Bar1),length($_)-length($Bar1))."\n":"$_\n";
	if (defined $Reads2 && $Reads2 ne "")
	{
		$_=<IN2>;
		my $out2=$_;
		$_=<IN2>;
		if (defined $Bar2 && $Bar2 ne "")
		{
			if ($_=~/^$Bar2(\S+)/)
			{
				$out2.="$1\n";
				$StatHash{'bar'}{$Bar2}++ if (defined $Stat);
			}
			else
			{
				$withBar=0;
			}
		}
		else
		{
			$out2.=$_;
		}
		$_=<IN2>;
		$out2.=$_;
		my $qual2=substr($_,0,length($Bar2)) if (defined $Bar2);
		$qualcheck=check_qual($qual2) if (defined $qual2 && $qual2 ne "");
		$_=<IN2>;
		chomp;
		$out2.=(defined $Barcode && defined $Bar2 && $Bar2 ne "")?substr($_,length($Bar2),length($_)-length($Bar2))."\n":"$_\n";
		if ($withIdx==1 && $withBar == 1 && $qualcheck==1) {
			print OUT2 $out2;
			$StatHash{'sel'}{'R2'}++ if (defined $Stat);
		}
	}
	if ($withIdx==1 && $withBar == 1 && $qualcheck==1) {
		print OUT1 $out1;
		$StatHash{'sel'}{'R1'}++ if (defined $Stat);
	}
}
close IN1;
close IN2;
close OUT1;
close OUT2;
if (defined $Stat)
{
	open (ST,">$StatOut") ;
	my $totalReads=0;
	if (exists $StatHash{'Idx'}) {
		print ST "Index\tReads\n";
		foreach my $index(sort{$StatHash{'Idx'}->{$b}<=>$StatHash{'Idx'}->{$a}}keys %{$StatHash{'Idx'}}) {
			print ST qq($index\t$StatHash{'Idx'}{$index}\n);
			$totalReads+=$StatHash{'Idx'}{$index};
		}
	}
	print ST "\n";
	if (exists $StatHash{'bar'}) {
		print ST "Barcode\tReads\n";
		foreach my $barcode(keys %{$StatHash{'bar'}}) {
			print ST qq($barcode\t$StatHash{'bar'}{$barcode}\n);
		}
	}
	print ST "\n";
	print ST "Total Reads: $totalReads\n";
	print ST qq(R1 selected reads: $StatHash{'sel'}{'R1'}\n) if (exists $StatHash{'sel'}{'R1'});
	print ST qq(R2 selected reads: $StatHash{'sel'}{'R2'}\n) if (exists $StatHash{'sel'}{'R2'});
	close ST;
}

sub check_qual{
	my $qual=shift;
	my @q=split //,$qual;
	for (my $i=0;$i<@q;$i++)
	{
		return 0 if (ord($q[$i])-33<$Qual);
	}
	return 1;
}


sub combinatorialSelect
{
	my ($num,$array_p) = @_;
	my @COMB=();
	my @array=sort{$a<=>$b}@$array_p;
	my @vernier;
	for my $i(0..($num-1))
	{
		$vernier[$i]=$i;
	}
	push @COMB,[@array[@vernier]];
	my $first=$vernier[0];
	SEED:
	if ($first!=@array-$num)
	{
		Shift_vernier((scalar @array),\@vernier);
		$first=$vernier[0];
		push @COMB,[@array[@vernier]];
		if (($vernier[-1]==@array-1)&&($first!=@array-$num))
		{
			Back_vernier(\@vernier);
			$first=$vernier[0];
			push @COMB,[@array[@vernier]];
		}
		goto SEED;
	}
	return @COMB;
}

sub Shift_vernier
{
	my ($size,$array_p)=@_;
	for (my $i=@$array_p-1;$i>=0;$i--)
	{
		if ($$array_p[$i]!=$size-@$array_p+$i)
		{
			$$array_p[$i]++;
			if (($i==0)||(($i>=1)&&($i<@$array_p-2)&&($$array_p[$i+1]>=$$array_p[$i]+1)))
			{
				for (my $j=$i+1;$j<@$array_p;$j++)
				{
					$$array_p[$j]=$$array_p[$j-1]+1;
				}
			}
			last;
		}
	}
}

sub Back_vernier
{
	my ($array_p)=@_;
	my $pre=$$array_p[@$array_p-1];
	for (my $i=@$array_p-2;$i>=0;$i--)
	{
		if ($pre-$$array_p[$i]>=2)
		{
			$$array_p[$i]++;
			if ($i<@$array_p)
			{
				for (my $j=$i+1;$j<@$array_p;$j++)
				{
					$$array_p[$j]=$$array_p[$j-1]+1;
				}
			}
			last;
		}
		$pre=$$array_p[$i];
	}
}

__END__