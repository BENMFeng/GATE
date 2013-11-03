#!/usr/bin/perl -w
#Author BENM <binxiaofeng\@gmail.com>
#Version: v1.0, 2011-01-13
use strict;
use warnings;
use Getopt::Long;

my %opts;
my ($PE,$N,$Q,$SD,$LenCutoff,$Trim,$Format,$Prefix,$Outfile,$Help);
GetOptions(\%opts,"PE"=>\$PE,"n:s"=>\$N,"q:i"=>\$Q,"sd:s"=>\$SD,
"f:s"=>\$Format,"len:i"=>\$LenCutoff,"trim"=>\$Trim,"prefix:s"=>\$Prefix,"out:s"=>\$Outfile,"help"=>\$Help);

my $usage = qq(
fastqcut.pl -- quality control for Solexa/illumina FASTQ
Usage: perl $0 [option] <Input files|pipein>
        --q <int>			set the minimal mean quality, which small than this value will be discarded
        --f <str>			quality format [SXIJL] follow FASTQ quality instruction, default: L
        --n <float|int>			set N's reads contained percentage (0~1) or length, default: 0.1 i.e. 10% length of reads are ambiguous bases
        --PE				input reads are in PE format, R1 shuffled with R2
        --len <int>			set the minimal bases pair of the reads, default: 35
        --sd <float>			set quality variance S.D.(STDEV), default: 5
        --trim				trim edge low qualtiy bases, recommand for retrieving high quality sequence
        --prefix <str>			set output filter/discard reads prefix name, if null it won't generate *.qcut.out file
        --out <file>			output high quality sequences after fastqcut, if null it will be STDOUT
        --help				help information
Example: perl solexa-qc.pl seq1.fastq -f L -n 5 -len 35 -sd 0.2 -trim -prefix seq1
Author : BENM <BinxiaoFeng\@gmail.com>
Version: 0.2.2 alpha, 2012-11-27

  FASTQ quality
  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
  !"#$%&'()*+,-./0123456789:;<=>?\@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |                              |                     |
 33                        59   64       73                            104                   126

 S - Sanger        Phred+33,  raw reads typically (0, 40)
 X - Solexa        Solexa+64, raw reads typically (-5, 40)
 I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
    (Note: See discussion above).
 L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

);

$N ||= 0.1;
$Q ||= 20;
$SD ||= 5;
$LenCutoff ||= 35;
$Format ||= "L";
my %qualth=('S'=>33,'X'=>59,'I'=>64,'J'=>66,'L'=>33);

die $usage if ($Help || !exists $qualth{$Format});
open (FT,">$Prefix.qcut.out") if (defined $Prefix);
open (OUT,">$Outfile") if (defined $Outfile);

while(<>) #@
{
	chomp;
	my $name1=$_;
	$_=<>; #seq1
	chomp;
	my $seq1=$_;
	my $tmp1=$seq1;
	$tmp1=~tr/N/\-/c;
	$tmp1=~s/\-//g;
	if ((($N<=1)&&(length($tmp1)>=length($seq1)*$N))||(($N>1)&&(length($tmp1)>=$N)))
	{
		$_=<>; #+
		$_=<>; #qual1
		if (defined $Prefix)
		{
			print FT ("$name1 N-len:",length($tmp1),"\n$seq1\n+\n$_");
		}
		if (defined $PE)
		{
			$_=<>; #@
			print FT $_ if (defined $Prefix);
			$_=<>; #seq2
			print FT $_ if (defined $Prefix);
			$_=<>; #+
			print FT $_ if (defined $Prefix);
			$_=<>; #qual2
			print FT $_ if (defined $Prefix);
		}
	}
	else
	{
		$_=<>; #+
		$_=<>; #qual1
		chomp;
		my $qual1=$_;
		my @q1=split "",$qual1;
		($seq1,$qual1)=trim_qcut($seq1,\@q1) if (defined $Trim);
		if (length($qual1)>=$LenCutoff)
		{
			@q1=split "",$qual1;
			my ($q_tot1,$add_square1,$q_over1)=(0,0,0);
			map{my $n=ord($_)-$qualth{$Format};$q_tot1+=$n;if ($n>=$Q){$q_over1++;}$add_square1+=$n*$n;}@q1;
			my $q_avg1=int($q_tot1/scalar(@q1));
			my $sd1=sqrt(($add_square1-scalar(@q1)*$q_avg1*$q_avg1)/(scalar(@q1)-1));
			if (($q_over1<scalar(@q1)/2)||(($q_avg1<$Q)&&($sd1>$SD)))   #quality control
			{
				if (defined $Prefix)
				{
					print FT ("$name1 q_avg:$q_avg1 Q$Q\_num:$q_over1 q_sd:$sd1\n$seq1\n+\n$qual1\n");
				}
				if (defined $PE)
				{
					$_=<>; #@
					print FT $_ if (defined $Prefix);
					$_=<>; #seq2
					print FT $_ if (defined $Prefix);
					$_=<>; #+
					print FT $_ if (defined $Prefix);
					$_=<>; #qual2
					print FT $_ if (defined $Prefix);
				}
			}
			else
			{
				if (defined $PE)
				{
					$_=<>; #@
					my $name2=$_;
					chomp;
					$_=<>; #seq2
					chomp;
					my $seq2=$_;
					my $tmp2=$seq2;
					$tmp2=~tr/N/\-/c;
					$tmp2=~s/\-//g;
					if ((($N<=1)&&(length($tmp2)>=length($seq2)*$N))||(($N>1)&&(length($tmp2)>=$N)))
					{
						if (defined $Prefix)
						{
							print FT ("$name1\n$seq1\n+\n$qual1\n$name2 N-len:",length($tmp2),"\n$seq2\n");
						}
						$_=<>; #+
						print FT $_ if (defined $Prefix);
						$_=<>; #qual2
						print FT $_ if (defined $Prefix);
					}
					else
					{
						$_=<>; #+
						$_=<>; #qual2
						chomp;
						my $qual2=$_;
						my @q2=split "",$qual2;
						if (defined $Trim)
						{
							($seq2,$qual2)=trim_qcut($seq2,\@q2);
							@q2=split "",$qual2;
						}
						my ($q_tot2,$add_square2,$q_over2)=(0,0,0);
						map{my $n=ord($_)-$qualth{$Format};if ($n>=$Q){$q_over2++}$q_tot2+=$n;$add_square2+=$n*$n;}@q2;
						my $q_avg2=int($q_tot2/scalar(@q2));
						my $sd2=sqrt(($add_square2-scalar(@q2)*$q_avg2*$q_avg2)/(scalar(@q2)-1));
						if (($q_over2>scalar(@q1)/2)||(($q_avg1<$Q)&&($sd1>$SD)))
						{
							if ((length($qual1)>=$LenCutoff)&&(length($qual2)>=$LenCutoff))
							{
								if (defined $Outfile)
								{
									print OUT "$name1\n$seq1\n+\n$qual1\n$name2\n$seq2\n+\n$qual2\n";
								}
								else
								{
									print "$name1\n$seq1\n+\n$qual1\n$name2\n$seq2\n+\n$qual2\n";
								}
							}
							elsif (defined $Prefix)
							{
								print FT ("$name1\n$seq1\n+\n$qual1\n$name2 seq_len:",length($qual2),"\n$seq2\n+\n$qual2\n");
							}
						}
						else
						{
							print FT "$name1\n$seq1\n+\n$qual1\n$name2 q_avg:$q_avg2 Q$Q\_num:$q_over2 q_sd:$sd2\n$seq2\n+\n$qual2\n";
						}
					}
				}
				else
				{
					if (length($qual1)>=$LenCutoff)
					{
						if (defined $Outfile)
						{
							print OUT "$name1\n$seq1\n+\n$qual1\n";
						}
						else
						{
							print "$name1\n$seq1\n+\n$qual1\n";
						}
					}
					elsif (defined $Prefix)
					{
						print FT ("$name1 seq_len:",length($qual1),"\n$seq1\n+\n$qual1\n");
					}
				}
			}
		}
		else
		{
			if (defined $Prefix)
			{
				print FT ("$name1 seq_len:",length($qual1),"\n$seq1\n+\n$qual1\n");
			}
			if (defined $PE)
			{
				$_=<>; #@
				print FT $_ if (defined $Prefix);
				$_=<>; #seq2
				print FT $_ if (defined $Prefix);
				$_=<>; #+
				print FT $_ if (defined $Prefix);
				$_=<>; #qual2
				print FT $_ if (defined $Prefix);
			}
		}
	}
}
close FT;
close OUT;

sub trim_qcut
{
	my ($seq_str,$qual_array)=@_;
	my ($j,$k)=(0,(length($seq_str)-1));
	my @seq_array=split "",$seq_str;
	for (my $i=0;$i<@seq_array;$i++)
	{
		$j=$i;
		last if ($seq_array[$i]!~/N/i && ($i<@seq_array-1 && $seq_array[$i+1]!~/N/i) && ($i<@seq_array-2 && $seq_array[$i+2]!~/N/i) );
	}
	for (my $i=@seq_array-1;$i>=0;$i--)
	{
		$k=$i;
		last if ($seq_array[$i]!~/N/i && ($i>=1 && $seq_array[$i-1]!~/N/i) && ($i>=2 && $seq_array[$i-2]!~/N/i));
	}
	my $seq=substr($seq_str,$j,($k-$j+1));
	my @qual_trim=@$qual_array[$j..$k];
	my $qual=join "",@qual_trim;
	
	($j,$k)=(0,(length($seq)-1));
	for (my $i=0;$i<@qual_trim;$i++)
	{
		$j=$i;
		last if (ord($qual_array->[$i])-$qualth{$Format}>=$Q);
	}
	for (my $i=@qual_trim-1;$i>=0;$i--)
	{
		$k=$i;
		last if (ord($qual_array->[$i])-$qualth{$Format}>=$Q);;
	}
	$seq=substr($seq,$j,($k-$j+1));
	$qual=join "",@qual_trim[$j..$k];
	return ($seq,$qual);
}