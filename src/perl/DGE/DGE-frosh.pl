#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $program_name = $1 if ( $0 =~ /([^\/]+)$/);
my $usage = qq(
$program_name -- TagSeq junior analysis script

Author: BENM <binxiaofeng\@gmail.com>
Version: v1.0 Mar-8-2010
 Update: v1.2 Aug-17-2010
         v1.3 Dec-9-2010

Usage: perl $program_name <IN:fastq_file|fasta_file|RawData_dir> [Option]
	-adaptor <str>		input adapter sequence [TCGTATGCCGTCTTCTGCTTG], default: it can be found by prgram
	-min <int>		minimum copy number, default: 1
	-max <int>		maximum copy number, default: 100000
	-qual <int>		quality type [solex|std], default: solexa
	-q_cut <int>		quality cutoff, default: 20
	-out <str>		output file, default its name is the same as rawfile
	-dir <str>		output dir, default: ./
	-enzyme <str>		set digested recognition site with the anchoring enzyme, from 5\' to 3\', default \(17bp, NlaIII\): \'CATG\'; \(16bp, DpnII\): \'GATC\'
	-trim <int-int>		trim 5\' or 3\' bases
	-polt			plot copy number distribution
	-help			help & usuage

Example:
perl $program_name RawData/

);
my ($Adaptor,$Min,$Max,$Qual,$Qual_cutoff,$Output,$Dir,$Enzyme,$Trim,$Plot,$Help);
my %opts;
GetOptions
(
	\%opts,
	"adaptor:s"=>\$Adaptor,
	"min:i"=>\$Min,
	"max:i"=>\$Max,
	"qual:s"=>\$Qual,
	"q_cut:s"=>\$Qual_cutoff,
	"out:s"=>\$Output,
	"dir:s"=>\$Dir,
	"enzyme:s"=>\$Enzyme,
	"trim:s"=>\$Trim,
	"plot"=>\$Plot,
	"help"=>\$Help
);
die($usage) if ((@ARGV==0)||($Help));

my ($left_trim,$right_trim) = (0,0);
if (defined $Trim)
{
	if ($Trim=~/-/)
	{
		($left_trim,$right_trim) = split /\-/,$Trim;
	}
	else
	{
		$right_trim = $1 if ($Trim=~/(\d+)/);
	}
}

$Min ||= 1;
$Max ||= 100000;
$Qual ||= "solexa";
$Qual_cutoff ||= 20;
$Dir ||= ".";
$Dir =~ s/\/$//;

my $norSize ||= 1000000;
my $tag_len ||= 17;

system "mkdir -p $Dir" if (-d $Dir);

# Solexa->Sanger quality conversion table
my @conv_table;
for (-64..64) {
	$conv_table[$_+64] = chr(int(33 + 10*log(1+10**($_/10.0))/log(10)+.499));
}

my @FileList;
if (-f $ARGV[0])
{
	push @FileList,$ARGV[0];
}
elsif (-d $ARGV[0])
{
	@FileList=glob("$ARGV[0]/*.fq");
	@FileList=glob("$ARGV[0]/*.fastq") if (@FileList==0);
}

my $Format ||= 0;  #0: fastq; 1: fasta
my $ADAPTOR;
my $File;
for (@FileList)
{
	$File=$_;
	$Output=(split /\//,$File)[-1];
	my @outname = split /\./,$Output;
	$Output = (join ".", @outname[0..(@outname-2)]) if (@outname>=2);
	$Output .= ".CopyNumber";
	if (defined $Adaptor)
	{
		$ADAPTOR=$Adaptor;
		chomp $ADAPTOR;
	}
	else
	{
		$ADAPTOR=get_adaptor($File);       #find adaptor sequence
	}
	if (!$ADAPTOR){
		my $ada=(split /\//,$File)[-1];
		open (AD,">$Dir/$ada.Adaptor") || die $!;
		print AD "Warning: No Adaptor can be found!\nUsing default adaptor sequence TCGTATGCCGTCTTCTGCTTG!\n";
		close AD;
		$ADAPTOR="TCGTATGCCGTCTTCTGCTTG" ;
	}
	$tag_len = ($ADAPTOR eq "TCGTATGCCGTCTTCTGCTTG") ? 17 : 16;
	if (defined $Enzyme)
	{
		$Enzyme = "" if ($Enzyme =~ /N/i);
	}
	else
	{
		$Enzyme = ($tag_len == 17) ? "CATG" : "GATC";
	}
	tag_hash($File);
}

sub get_adaptor
{
	my $file=shift;
	my $adaptor="";
	open (IN,$file) || die $!;
	$_=<IN>;
	if ($_=~/^\@/)
	{
		$Format = 0;
	}
	else
	{
		$Format = 1;
	}
	my $Seq=<IN>;
	chomp $Seq;
	my $Length =length($Seq);
	if ($Length<=$tag_len)
	{
		return $adaptor;
	}
	my %Sixteen;
	my %Seventeen;
	$Sixteen{uc(substr($Seq,16))}+=1;
	$Seventeen{uc(substr($Seq,17))}+=1;
	my $reads=0;
	if ($Format==0)
	{
		<IN>;<IN>;
	}
	while(<IN>)
	{
		$Seq=<IN>;
		chomp $Seq;
		my $sixteen=uc(substr($Seq,16));
		my $seventeen=uc(substr($Seq,17));
		$Sixteen{$sixteen}+=1;
		$Seventeen{$seventeen}+=1;
		if ($Format==0)
		{
			<IN>;<IN>;
		}
		$reads++;
		last if ($reads>999999);
	}
	close IN;
	my @array=(sort{$Sixteen{$b}<=>$Sixteen{$a}}keys %Sixteen);
	my $A=$array[0];
	my ($B,$C)=(sort{$Seventeen{$b}<=>$Seventeen{$a}}keys %Seventeen);
	if ((!defined $C)||(($Seventeen{$B}>1.2*$Seventeen{$C})&&($B=~/$A$/)&&($C=~/$A$/)))
	{
		$adaptor=$B;
	}
	else
	{
		$adaptor=$A;
	}
	my $ada=(split /\//,$File)[-1];
	open (AD,">$Dir/$ada.adaptor") || die $!;
	print AD "Find adaptor: $adaptor\n";
	close AD;
	return $adaptor;
}

sub tag_hash
{
	my $file=shift;
	open (IN,$file) || die $!;
	my %Tag_hash;
	my ($reads,$noAdap,$NsReads,$lowQ)=(0,0,0,0);
	my $out=(split /\//,$file)[-1];
	open (LOW,">$Dir/$out.Lowquality") || die $!;
	while(<IN>)
	{
		my $name=$_;
		chomp $name;
		my $Seq=<IN>;
		chomp $Seq;
		my $tag=uc($Seq);
		if (length($tag)>$tag_len)
		{
			substr($tag,0,$left_trim,"") if ($left_trim>0);
			substr($tag,length($tag)-$right_trim,$right_trim,"") if ($right_trim>0);
			if (length($tag)>$tag_len)
			{
				if ($tag!~/$ADAPTOR/)
				{
					$noAdap++;
					if ($Format==0)
					{
						<IN>;<IN>;
					}
					print LOW "$name\t$Seq\t*\t$_";
					next;
				}
				$tag=~s/$ADAPTOR$//;
			}
			if ((length($tag)<$tag_len)||($tag=~/N/i))      #filter N's reads
			{
				$NsReads++;
				if ($Format==0)
				{
					<IN>;$_=<IN>;
				}
				print LOW "$name\t$Seq\t-\t$_";
				next;
			}
		}
		if ($Format==0)
		{
			<IN>;$_=<IN>;
			chomp;
			my @q=split "",substr($_,0,length($tag));
			my ($q_tot,$q_over)=(0,0);
			for (my $i=0;$i<length($tag);$i++)
			{
				my $j=$q[$i];
				if ($Qual=~/solexa/i)
				{
					$q_tot += (ord($j)-64);
					$q_over++ if ((ord($j)-64)>=$Qual_cutoff);
				}
				else
				{
					$q_tot += (ord($j)-33);
					$q_over++ if ((ord($j)-64)>=$Qual_cutoff);
				}
			}
			my $q_avg=int($q_tot/length($tag));
			if (($q_avg<$Qual_cutoff)||($q_over<length($tag)/2))   #quality control
			{
				$lowQ++;
				print LOW "$name\t$Seq\t$q_avg\t$_\n";
				next;
			}
		}
		$Tag_hash{$tag}++;
		$reads++;
	}
	close IN;
	close LOW;
	my $order=0;
	my $pre=0;
	my %Stat_CopyNumber;
	my ($tot_a,$tot_v,$min,$min_v,$max,$max_v,$avg,$avg_v,$median,$median_v,$n50,$n50_v)=(0,0,-1,0,0,0,0,0,0,0,0,0);

	my $cleanNor = $reads - $noAdap - $NsReads - $lowQ;
	open (ST,">$Dir/$out.Stat.log") || die $!;
	print ST ("$out - Reads: ",Thousandths($reads),"\n");
	print ST ("$out - Clean Tag-reads: ",Thousandths($cleanNor),"\n");
	print ST ("$out - No Adapters Reads: ",Thousandths($noAdap),"\n");
	print ST ("$out - Ns Reads: ",Thousandths($NsReads),"\n");
	print ST ("$out - LowQuality Reads: ",Thousandths($lowQ),"\n");
	print ST ("$out - Total filterred Reads: ",Thousandths($noAdap+$NsReads+$lowQ),"\n");
	#print "Tag\tCopyNumber\tOrder\tTPM(Transcript Per Million)\n";
	open (OUT,">$Dir/$Output") || die $!;
	print OUT "Tag\tCopyNumber\tOrder\tNormalized_Frequency\n";

	foreach my $seq(sort{$Tag_hash{$b}<=>$Tag_hash{$a}}keys %Tag_hash)
	{
		next if (($Tag_hash{$seq}<$Min)||($Tag_hash{$seq}>$Max));
		$order++ unless ($pre==$Tag_hash{$seq});
		print OUT ("$Enzyme$seq\t$Tag_hash{$seq}\t$order\t",sprintf("%.2f",  $norSize / $cleanNor * $Tag_hash{$seq}),"\n");
		$Stat_CopyNumber{$Tag_hash{$seq}}++;
		$tot_v+=$Tag_hash{$seq};
		$pre=$Tag_hash{$seq};
	}
	close OUT;

	my ($acc,$aca)=(0,0);
	$tot_a=(keys %Stat_CopyNumber);
	$avg_v=$tot_v/$tot_a;
	my $copy_distribution;
	foreach my $x(sort{$a<=>$b} keys %Stat_CopyNumber)
	{
		$aca+=$Stat_CopyNumber{$x};
		$acc+=$tot_a*$Stat_CopyNumber{$x};
		if($min==-1)
		{
			$min=$x;
			$min_v=$Stat_CopyNumber{$x};
		}
		$max=$x;
		$max_v=$Stat_CopyNumber{$x};
		if (($median==0)&&($aca>=$tot_a/2))
		{
			$median=$x;
			$median_v=$Stat_CopyNumber{$x};
		}
		if (($n50==0)&&($acc>=$tot_v/2))
		{
			$n50=$x;
			$n50_v=$Stat_CopyNumber{$x};
		}
		$avg=$x*$Stat_CopyNumber{$x}/$tot_v;
		$copy_distribution.=("$x\t$Stat_CopyNumber{$x}\t".sprintf("%.4f",(100*$x*$Stat_CopyNumber{$x}/$tot_v))."%\n");
	}
	print ST ("$out - Total Identifiable(Unique) Tags: ",Thousandths(scalar(keys %Tag_hash)),"\n");
	print ST ("$out - Total Raised CopyNumber: ",Thousandths($tot_a)," Total Frequency: ",Thousandths($tot_v),"\n");
	print ST ("$out - The Lowest CopyNumber: ",Thousandths($min)," Frequency: ",Thousandths($min_v),"\n");
	print ST ("$out - The Highest CopyNumber: ",Thousandths($max)," Frequency: ",Thousandths($max_v),"\n");
	print ST ("$out - The Average CopyNumber: ",Thousandths($avg)," Frequency: ",Thousandths($avg_v),"\n");
	print ST ("$out - The Median CopyNumber: ",Thousandths($median)," Frequency: ",Thousandths($median_v),,"\n");
	print ST ("$out - The N50 CopyNumber: ",Thousandths($n50)," Frequency: ",Thousandths($n50_v),"\n");
	print ST ("------------------------------------------------------------\n");
	print ST ("CopyNumber\tFrequency\tPercentage\n");
	print ST $copy_distribution;
	close ST;
	if (defined $Plot)
	{
		open (PL,">$Dir/$out.CopyNumber.distribution");
		print PL $copy_distribution;
		close PL;
		gnuplot(1,2,$out,"$Dir/$out.CopyNumber.distribution");
	}
}

sub Thousandths
{
	my $num=shift;
	my $thousands=$num;
	next if ($thousands=~/\,/);
	if (sprintf("%.4f",$thousands)=~/(\d+)\.(\d+)/)
	{
		my $udecimal=reverse($1);
		my $ddecimal=$2;
		$udecimal=~s/(\d{3})/$&.","/eg;
		$udecimal=~s/\,$//;
		$udecimal=reverse($udecimal);
		$thousands=($ddecimal eq "0000")? $udecimal : $udecimal.".".$ddecimal;
	}
	else
	{
		$thousands=reverse($thousands);
		$thousands=~s/(\d{3})/$&.","/eg;
		$thousands=~s/\,$//;
		$thousands=reverse($thousands);
	}
	return($thousands);
}

sub gnuplot
{
	my ($x,$y,$name,$file)=@_;
	my $fh="$file.gnuplot.dat";
open(PL, ">$fh") || die $!;
print PL qq(
reset;
set title '$name CopyNumber Distribution';
set xlabel 'CopyNumber';
set ylabel 'Frequency';
set nokey;
set terminal svg;
set output "$file.svg";
plot "$file" u $x:$y w boxes;
);
close (PL);
system "gnuplot $fh";
}
