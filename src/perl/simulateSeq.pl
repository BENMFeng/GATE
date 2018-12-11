#!/usr/bin/perl -w

=head1 Name

Copyright (c) 2009-2011 BENM(Binxiao) Feng                            
All Rights Reserved                                                   
Send all comments to BENM - BinxiaoFeng@gmail.com                     

This program is free software: you can redistribute it and/or modify  
it under the terms of the GNU General Public License as published by  
the Free Software Foundation, either version 3 of the License, or     
(at your option) any later version.                                   
This program is distributed in the hope that it will be useful,       
but WITHOUT ANY WARRANTY; without even the implied warranty of        
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
GNU General Public License for more details.                          
You should have received a copy of the GNU General Public License     
along with this program.  If not, see <http://www.gnu.org/licenses/>. 

simulateSeq.pl -- PERL script for simulating NGS's sequence with variation & paired-end reads (single-end reads supported);
			it can output Solexa, SOLiD, 454, 3730, std(sanger) all types format seqencing reads & quality data.

=head1 Description

This are three main functions in this program:

(1) This program will be used to make random diversity dna sequence with a set
happening frequecy It can make an offspring(testing) sequence with mutation
according to a ancestor(reference) sequence, which could either be given by the
input file, or randomly generated; 

(2)Simulate sequence with some structure vairation cases
The structure vairation(SV) can be comprised of Deletion, Insertion, Inversion,
Duplication, Translocation, each of these 5 main SV elements, or their combination.
And it can be simulate sequence with SNPs, Indel, base bias too;

(3)Simulate singled-end and paired-end (mate pairs) reads
Meanwhile it can output the SE or PE reads according to the simulated testing sequence,
and it can output with fasta or fastq format by modified the option of "--type"

More detail and options you can look for below Usage and Exmple.

=head1 Version
Author: BENM <binxiaofeng@gmail.com>
Birthday: 2007-12-26
Update: May,26th,2011
Version: 0.2.2

=head1 Usage

	perl simulateSeq.pl [options] [ seq.fa | STDIN ] [>SV.log]

  ***** Sequence parameters *****
	--Ref <outfile>    make the reference(ancestor) sequence:[chrS.fasta]
	--Seq <outfile>    make the testing(offspring) diploid sequence:[Seq.fasta]
	--num <int>        set number of testing sequences [1]
	--SeqName <str>    set simulated testing sequence name: [chrS]
	--Reverse          simulate the testing reversed sequence compared to reference
	--Complement       simulate the testing complement sequence compared to reference
	--RC               simulate the reversed and complmented sequence compared to reference
	--len <int>        set length of the sequence, default: whole sequence
	--BAC <num1-num2-num3>       make a number of BAC, num1 is BAC's number, num2 is BAC's length, num3 is sd of BAC's lenth
	--vector <str>     add vector sequence to BAC clone sequencing
	--enzyme <str>     restriction digestive enzyme site sequence [File or STDIN, example EcoRI: AATT]
	--qual_th <num1-num2-num3-num4>    set higher quality scores theshold of 5' end bases: <3' reads split position-5'endquality_theshold-3'end quality theshold-quality variance>, default: 5-30-20-2

  ***** Mutation parameters *****
	--div <float>      set divergency rate of the offspring sequences, it is used for SNPs and Indels [0]
	--changeb <str> only change the base by setting, Case-sensitive, it is used for methylation
	--Del <int[-start:int-end:int-size:int]>    set Deletion vairation number 
	--Ins <int[-start:int-end:int-size:int]>    set Insertion vairation number 
	--Inv <int[-start:int-end:int-size:int]>    set Inversion vairation number 
	--Dup <int[-start:int-end:int-size:int]>    set Duplication vairation number
	--Tran <int[-start:int-end:int-size:int]>   set Translocation vairation number
	--size <num>        if  the above vairation seting hasn't  defined the size, it can use the common
	                    vairation size setting here.
	                    For example: ~1000, "~" means around 1000, has 20 deviation, which set by the "ISSD"
	                    option, if no this symbol, it will be precise size.
	--mutation          Output muation report

  ***** Read parameters *****
	--PE <int|numX>     make a number of pair-end reads for each sequence [1 pairs = 2 reads]
	--ISPE <int>        set insert size (outer coordinate) of the pair-end reads [135]
	--ISSD <float>      set insert size standard deviation [10]
	--RL <int-int>      set reads' length of each end of pairs,(if "null" or "0" means SE, both of num
	                    can't be 0) [35-35]
	--ReadsSD <float>   set reads length standard deviation, it used to long reads[0]
	--name <str>        set the pair-end reads name, (In fa format, it will give the general name),
	                    default: BENM_2011_0526_022
	--type < solexa[fa|fq|stdfq] || solid[fa|fq|stdfq] || 454[fa|fq|stdfq] || 3730[fa|fq|stdfq] || stdfq>
	                    output the sequence with simulate fasta or fastq format [stdfq]
	--header <str>      Solid header base [G]
	--adapter <infile>  Adapter sequence file, sequence name should use 5' and 3' to mark different adapter 
	--normal_strand     set the normal strands of two reads of PE [FR],SOLiD type:[RF], 454 type[RF];
	--circle            it set for circle chromosome(Chloroplast Mitochondrion and BAC) generate PE data

  ***** Display parameters *****
	--outdir <str>      set the output directory for storing  result [Simulated-result]
	--display <int>     set the fasta file displayed how many bases in each lines [50]
	--Example           output Example to screen
	--verbose           output verbose information to screen
	--help              output help information to screen

=head1 Exmple

	check out "perl simulateSeq.pl --Example"
 
=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use constant PI=>3.14159265;

##get options from command line into variables and set default values
my ($Reference,$Sequence,$SeqName,$Reverse,$Complement,$RC,$Length,$Divergency,$ChangeB,$Number,$Qual_th,
	$Deletion,$Insertion,$Inversion,$Duplication,$Translocation,$Mutation,$Size,$BAC,$Vector,
	$InsertSizePE,$InsertSizeSD,$ReadsSD,$Type,$Header,$Adapter,$PEreads,$PE,$Name,$NormalStrand,
	$Enzyme,$Circle,$Outdir,$Display,$Example,$Verbose,$Help);
my %opts;
GetOptions
(
	\%opts,
	"Ref:s"=>\$Reference,
	"Seq:s"=>\$Sequence,
	"SeqName:s"=>\$SeqName,
	"Reverse"=>\$Reverse,
	"Complement"=>\$Complement,
	"RC"=>\$RC,
	"len:i"=>\$Length,
	"div:f"=>\$Divergency,
	"changeb:s"=>\$ChangeB,
	"num:i"=>\$Number,
	"qual_th:s"=>\$Qual_th,
	"Del:s"=>\$Deletion,
	"Ins:s"=>\$Insertion,
	"Inv:s"=>\$Inversion,
	"Dup:s"=>\$Duplication,
	"Tran:s"=>\$Translocation,
	"mutation"=>\$Mutation,
	"size:s"=>\$Size,
	"ISPE:i"=>\$InsertSizePE,
	"ISSD:s"=>\$InsertSizeSD,
	"ReadsSD:s"=>\$ReadsSD,
	"RL:s"=>\$PEreads,
	"PE:s"=>\$PE,
	"name:s"=>\$Name,
	"type:s"=>\$Type,
	"header:s"=>\$Header,
	"adapter:s"=>\$Adapter,
	"normal_strand:s"=>\$NormalStrand,
	"circle"=>\$Circle,
	"BAC:s"=>\$BAC,
	"vector:s"=>\$Vector,
	"enzyme:s"=>\$Enzyme,
	"outdir:s"=>\$Outdir,
	"display:i"=>\$Display,
	"Example"=>\$Example,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
if (defined $Example)
{
    my $example = qq{
Example:
    0.easy make sequences
    
    \$ perl simulateSV-PE.pl --Seq Seq.fa ACGTCGCTTGGCNNGTACCGTA CGTAGTCCGTAAAGTTCC ACGTAGTCAGTCC
    
    1.make one sequence with 500000bp bases according to reference fasta file
 
    \$ perl simulateSV-PE.pl ./test-data/chrS.fa --Seq Seq.fasta --len 500000
 
    2.simulate two sequence with 500000bp bases and with divergency rate of 0.1
    according to reference fasta file, and two sequence will be stored in the
    same file
 
    \$ perl simulateSV-PE.pl ./test-data/chrS.fa --div 0.01 --Seq Seq.fasta --len 500000 --num 2
 
    3.simulate one reference sequence and two testing sequence
 
    \$ perl simulateSV-PE.pl --Ref chrS.fa --div 0.01 --Seq Seq.fa --len 500000 --num 2
 
    4.simulate one reference sequence and one testing sequence and the testing
    sequence has "GC" bias, and display 100 bases in each lines
 
    \$ perl simulateSeq.pl --Ref chrS.fa --div 0.05 --Seq Seq.fa --len 500000 --num 1 --changeb ATat --display 100
 
    5.simulate a number[8X] of paired-end reads with set insert size[440] and
    the insert size will be distributed with a set standard deviation[20]
    according to approximative Gauss distribution model
 
    \$ perl simulateSeq.pl ./test-data/chrS.fa --Seq Seq.fa --PE 8X --ISPE 440 --ISSD 20 --type stdfq
 
    6.simulate sequence with some structure varirations such as deletion,
    insertion, inversion cases and output 50000 paired-end reads with fastq
    format and the information of each case will be ouput in the logfile
 
    \$ perl simulateSeq.pl ./test-data/chrS.fa --Seq Seq.fa --div 0 --Del 1 --Ins 1 --Inv 1 --PE 50000 --size 1000 --type stdfq > SV.log

    7.simulate sequence with some structure varirations such as deletion,
    insertion, inversion cases with precise position and size and output 50000
    paired-end reads with fastq format and the information of each case will be
    ouput in the logfile
 
    \$ perl simulateSeq.pl ./test-data/chrS.fa --Seq Seq.fa --div 0 --Del 1-start:100-end:1100 --Ins 1-start:2000-size:1000 --Inv 1-end:50000-size:1000 --Dup 1 --Tran 1-size:2000 --PE 50000 --type stdfq > SV.log
 
    8.simulate sequence with some structure varirations and output 3X pair-ends
    reads with normal strand of "RF"
 
    \$ perl simulateSeq.pl --Ref chrS.fa --Seq Seq.fa --div 0 --len 5000000 --Del 1-start:10000-end:14000 --Ins 1-start:2000-size:1500 --Inv 3 --Dup 5 --size ~1300 --PE 3X --normal_strand RF >SV.log
 
    9.simulate 3X single-end reads in 454fa format
 
    \$ perl simulateSeq.pl --Ref chrS.fa --div 0 --len 5000000 -PE 3X --RL 0-450 --ReadsSD 10 --name 454 --type 454fa --outdir SE
    
    10.simulate reads with adapter sequence with random length(1~10)
    
    \$ perl simulateSeq.pl --Ref chrS.fa --len 5000000 --adapter adapter.fa -PE 3X --RL 450-450 --ReadsSD 20 --type 3730stdfq
    };
    die "$example\n";
}

##default settings
$SeqName ||= "chrS";
$Number ||= 1;
$Type ||= "stdfq";
$InsertSizePE ||= 135;
$InsertSizeSD ||= 10;
$PEreads ||="35-35";
$Name ||= "BENM_2011_0526_022";
$Outdir ||= "Simulated-result";
$Display ||= 50;
$NormalStrand ||= (($Type=~/solid/)||($Type=~/454/)||($InsertSizePE>600)) ? "RF" : "FR";
$Divergency ||= 0;
$ReadsSD ||= 0;
my @qual_th=split/[\-\;\,\.\_]/,$Qual_th if (defined $Qual_th);
$qual_th[0] = 5 unless ((defined $qual_th[0])&&($qual_th[0]>=0));
$qual_th[1] = 30 unless ((defined $qual_th[1])&&($qual_th[1]>=0));
$qual_th[2] = 20 unless ((defined $qual_th[2])&&($qual_th[2]>=0));
$qual_th[3] = 2 unless ((defined $qual_th[3])&&($qual_th[3]>=0));


# Solexa->Sanger quality conversion table
my @conv_table;
for (-64..64) {
	$conv_table[$_+64] = chr(int(33 + 10*log(1+10**($_/10.0))/log(10)+.499));
}

# Solid color code
my @code = ([0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]);
my @bases = qw(A C G T);
my @othercode = (".",4,5,6);
my %colcode = ();
my %decode = ();
foreach my $i(0..3)
{
	foreach my $j(0..3)
	{
		$decode{$code[$i]->[$j]} -> {$bases[$i]} = $bases[$j];
	}
}

foreach my $i(0..3)
{
	foreach my $j(0..3)
	{
		$colcode{"$bases[$i]$bases[$j]"} = $code[$i]->[$j];
	}
}

my @base= ("A","G","C","T"); 
for (my $i=0; $i<10*(1-$Divergency); $i++)
{
	push @base,("A","G","C","T");
}
#push @base,("a","c","g","t","a","c","g","t","N","n","-"," "); ## indel: "" stand for deletion, "-" stand for insertion
push @base,("N","-"); ## indel: "" stand for deletion, "-" stand for insertion

if (defined $Length)
{
	if ($Length=~/(\S+)k/i)
	{
		$Length=$1*1000;
	}
	if ($Length=~/(\S+)M/i)
	{
		$Length=$1*1000*1000;
	}
	if ($Length=~/(\S+)G/i)
	{
		$Length=$1*1000*1000*1000;
	}
}

my ($RL1,$RL2);
if (defined $PEreads)
{
	chomp $PEreads;
	($RL1,$RL2)=split /\-/,$PEreads;
	$RL1 = 0 if ($RL1 !~ /\d+/);
	$RL2 = 0 if ($RL2 !~ /\d+/);
	die `pod2text $0` if (($RL1==0)&&($RL2==0));
}
if ($InsertSizePE<$RL1+$RL2)
{
	$InsertSizePE+=$RL1+$RL2;
}

my $Enzyme_r;
if (defined $Enzyme)
{
	if (-f $Enzyme)
	{
		open (E,$Enzyme) || die "Can't open Enzyme file for reading\n";
		while(<E>)
		{
			if ($_!~/^>/)
			{
				s/\s+//g;
				$Enzyme=$_;
			}
		}
		close E;
	}
	$Enzyme_r=reverse($Enzyme);
}

my @Vector_seq;
if (defined $Vector)
{
	my $vector_seq='';
	open (V,$Vector) || die "Can't open vector sequence\n";
	while(<V>)
	{
		if (/>/)
		{
			(push @Vector_seq,$vector_seq) if ($vector_seq ne '');
			$vector_seq='';
		}
		else
		{
			s/\s+//g;
			$vector_seq.=$_;
		}
	}
	my $vector_seq_rc=reverse($vector_seq);
	$vector_seq_rc=~tr/ACGTacgt/TGCAtgca/;
	if (defined $Enzyme)
	{
		($vector_seq=$Enzyme.$2.$1) if ($vector_seq=~/(\w+)$Enzyme(\w+)/);
		($vector_seq_rc=$Enzyme.$2.$1) if ($vector_seq_rc=~/(\w+)$Enzyme(\w+)/);
	}
	(push @Vector_seq,$vector_seq) if ($vector_seq ne '');
	(push @Vector_seq,$vector_seq_rc) if ($vector_seq_rc ne '');
	if (@Vector_seq!=2)
	{
		die "vector or enzyme erro\n";
	}
	if (defined $Verbose)
	{
		Display_seq(\$vector_seq);
		print STDERR ">Vector+\n$vector_seq\n";
		Display_seq(\$vector_seq_rc);
		print STDERR ">Vector-\n$vector_seq_rc\n";
	}
	close V;
}

my ($BAC_num,$BAC_length,$BAC_sd);
my $BAC_count=1;
my @BAC_len;
if (defined $BAC)
{
	chomp $BAC;
	($BAC_num,$BAC_length,$BAC_sd)=split /\-/,$BAC;
	my $bac_num = ($BAC_num=~/\d+X/i) ? 100000 : $BAC_num;
	if ((($BAC_num=~/\d+/)||($BAC_num=~/\d+X/i))&&($BAC_length=~/\d+/)&&($BAC_sd=~/\d+/))
	{
		@BAC_len=gauss($bac_num,$BAC_length,$BAC_sd);
	}
	else
	{
		die "BAC setting illegal : --BAC <num1-num2-num3>\n";
	}
}

die `pod2text $0` if ($Help);
die `pod2text $0` unless ((@ARGV>0)||(defined $Reference));

##global variables
my $seq;

my $RefFile = (defined $Reference) ? $Reference : "chrS.fasta";
my $SeqFile =  (defined $Sequence) ? $Sequence : "Seq.fasta";

my $print_title=0;

if ($Divergency>=1)
{
	print STDERR "Oops!!!\nIt is God's forbidden domination!\n\t\t\t\t  --BENM\n";
	sleep(10);
	system "cat $0";
	system "clear";
	system "telnet towel.blinkenlights.nl";
	exit;
}

#######################################################################################
#---------------------------------- Main  Function ----------------------------------#
######################################################################################

##simulate reference or testing sequence and ouput the files in the ./simulated-result/ directory
`mkdir -p $Outdir` unless (-d $Outdir);

# adapter join length

my @rand_len=(0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5,5,6,6,7,8,9,10);
#read adapter file
my (@adapter1,@adapter2);
if (defined $Adapter)
{
	open (IN,$Adapter) || die "Fail to open adapter file: $Adapter for reading\n";
	my ($adapter_num,$adapter_len,$adapter_sd);
	if ((defined $RL1)&&(defined $RL2))
	{
		$adapter_len=(sort{$a<=>$b}($RL1,$RL2))[0];
	}
	elsif ((defined $RL1)&&(!defined $RL2))
	{
		$adapter_len=$RL1;
	}
	elsif ((!defined $RL1)&&(defined $RL2))
	{
		$adapter_len=$RL2;
	}
	my ($adapter_name,$adapter_seq)=('','');
	while(<IN>)
	{
		if ($_=~/^\>(.*)/)
		{
			if ((!defined $RL1)&&(!defined $RL2))
			{
				$adapter_len=int(length($adapter_seq)/2.5);
			}
			$adapter_num=int(length($adapter_seq))/$adapter_len*2;
			$adapter_sd=int(length($adapter_seq)/6+0.5);
			if ($adapter_seq ne '')
			{
				if (($adapter_name=~/5/)||($adapter_name=~/left/i))
				{
					push @adapter1,$adapter_seq;
					push @rand_len, gauss($adapter_num,$adapter_num,$adapter_sd);
				}
				elsif (($adapter_name=~/3/)||($adapter_name=~/right/i))
				{
					push @adapter2,$adapter_seq;
					push @rand_len, gauss($adapter_num,$adapter_num,$adapter_sd);
				}
			}
			$adapter_name=$1;
			$adapter_seq='';
		}
		else
		{
			s/\s+//g;
			$adapter_seq.=$_;
		}
		if (eof)
		{
			if ((!defined $RL1)&&(!defined $RL2))
			{
				$adapter_len=int(length($adapter_seq)/2.5);
			}
			$adapter_num=int(length($adapter_seq))/$adapter_len*2;
			$adapter_sd=int(length($adapter_seq)/6+0.5);
			if ($adapter_seq ne '')
			{
				if (($adapter_name=~/5/)||($adapter_name=~/left/i))
				{
					push @adapter1,$adapter_seq;
					push @rand_len, gauss($adapter_num,$adapter_num,$adapter_sd);
				}
				elsif (($adapter_name=~/3/)||($adapter_name=~/right/i))
				{
					push @adapter2,$adapter_seq;
					push @rand_len, gauss($adapter_num,$adapter_num,$adapter_sd);
				}
			}
		}
	}
	close IN;
}

#make random seqequence and PE or read from file
my $chr_num=0;
if((@ARGV == 0)||(defined $Reference))
{
	$seq=make_seq($Length);
	my $BAC_number=$BAC_num;
	make_BAC($seq,$BAC_number) if (defined $BAC);
}
else
{
	my %BAC_seq;
	my $w=0;
	
	for (my $i=0;$i<@ARGV;$i++)
	{
		if (-f $ARGV[$i])
		{
			my $length_cutoff=0;
			$Length=0 if (!defined $opts{len});
			open (IN,$ARGV[$i]) || die $!;
			$/=">"; <IN>; $/="\n";
			while (<IN>) {
				$/=">";
				$SeqName=$_;
				$SeqName=~s/\s+$//g;
				$seq = <IN>;
				chomp $seq;
				$seq =~ s/\s//g;
				next if ($seq eq "");
				$length_cutoff += length($seq);
				$Length += length($seq) if (!defined $opts{len});
				if ((defined $Length)&&($length_cutoff>$Length))
				{
					substr($seq,0,($length_cutoff-$Length),"");
				}
				push @{$BAC_seq{length($seq)}},$SeqName;
				$w++;
				if ((defined $Sequence)||($Divergency>0))
				{
					makeSeq($SeqName,$seq);
				}
				else
				{
					if (defined $PE)
					{
						makePE($PE,$seq,Complement($seq),$SeqName);
					}
				}
				$/="\n";
				last if ((!defined $opts{len})&&($length_cutoff>$Length));
			}
			close IN;
		}
		else
		{
			$seq=$ARGV[$i];
			chomp $seq;
			$Length += length($seq) if (!defined $Length);
			my $seq_name = (@ARGV>1) ? $SeqName."_".($i+1) : $SeqName;
			push @{$BAC_seq{length($seq)}},$seq_name;
			$w++;
			if ((defined $Sequence)||($Divergency>0))
			{
				makeSeq($seq_name,$seq);
			}
			else
			{
				makePE($PE,$seq,Complement($seq),$seq_name) if (defined $PE);
			}
		}
	}
	$chr_num=$w;
	if ((defined $BAC)&&(keys %BAC_seq>0))
	{
		my %BAC_seq_num;
		my $seq_name;
		foreach my $l(sort{$a<=>$b} keys %BAC_seq)
		{
			if ($l<$BAC_length)
			{
				$w--;
				next;
			}
			for (my $j=0;$j<@{$BAC_seq{$l}};$j++)
			{
				$seq_name=${$BAC_seq{$l}}[$j];
				chomp $seq_name;
				if ($BAC_num=~/(\S+)X/i)
				{
					$BAC_seq_num{$seq_name}=int($l*$1/$BAC_length+0.5);
				}
				else
				{
					$BAC_seq_num{$seq_name}=int ($BAC_num/$w);
				}
			}
		}
		$BAC_seq_num{$seq_name}+=($BAC_num-(int($BAC_num/$w))*$w) if ($BAC_num-(int($BAC_num/$w))*$w>0);
		
		for (my $i=0;$i<@ARGV;$i++)
		{
			if (-f $ARGV[$i])
			{
				my $length_cutoff=0;
				open (IN,$ARGV[$i]);
				$/=">"; <IN>; $/="\n";
				while (<IN>) {
					$/=">";
					$SeqName=$_;
					$SeqName=~ s/\s+$//g;
					$seq = <IN>;
					chomp $seq;
					$seq=~s/\s//g;
					next if ($seq eq "");
					if (exists $BAC_seq_num{$SeqName})
					{
						print STDERR "simulating BAC sequence of $SeqName...\n" if($Verbose);
						make_BAC($seq,$BAC_seq_num{$SeqName});
					}
					$/="\n";
				}
				close IN;
			}
			else
			{
				$seq=$ARGV[$i];
				chomp $seq;
				my $seq_name = (@ARGV>1) ? $SeqName."_".($i+1) : $SeqName;
				if (exists $BAC_seq_num{$seq_name})
				{
					print STDERR "simulating BAC sequence of $seq_name...\n" if($Verbose);
					make_BAC($seq,$BAC_seq_num{$seq_name});
				}
			}
		}
	}
}

$Length += length($seq) if (!defined $Length);
my $sequences=($Number>1) ? "sequences" : "sequence";
print STDERR "-- $Number $sequences, $chr_num chromosomes, total length: $Length\n" if($Verbose);

if(defined $Reference)
{
	if (length($seq)>$Length)
	{
		$seq = substr ($seq,0,$Length);
	}
	elsif (length($seq)<$Length)
	{
		$seq = substr ($seq,0,length($seq));
		$seq .= make_seq($Length-length($seq));
	}
	
	$seq = reverse ($seq) if ( defined ($Reverse || $RC) );
	$seq =~ tr/ACGTacgt/TGCAtgca/ if ( defined ($Complement || $RC) );
	
	open OUT1,">$Outdir/$RefFile" or die "Can not create file:$!";
	my $str = $seq;
	Display_seq(\$str,$Display); 
	print OUT1 (">$SeqName\n".$str);
	my $ref_length=length($seq);
	print STDERR "simulated reference sequence done!\tLength:$ref_length\n" if($Verbose);
	makeSeq($SeqName,$seq) if ((defined $Sequence)||(defined $PE));
}
close OUT1;

system "sort -n -k 1,2 $Outdir/Mutation.txt > $Outdir/tmp.BENM && mv $Outdir/tmp.BENM $Outdir/Mutation.txt" if ((defined $Mutation)&&(-f "$Outdir/Mutation.txt"));

#####################################################################################
#----------------------------------  Subroutine ---------------------------------- #
#####################################################################################

##simulate Sequence, SV & PE;
## usage: makeSeq($seq_name,$seq);
#############################################
sub makeSeq
{
	my ($seq_name,$seq)=@_;
	chomp $seq;
	if ((defined $Sequence)||(defined $PE))
	{
		open OUT2,">>$Outdir/$SeqFile" if (defined $Sequence);
		
		my $seq_length=length($seq);
		my @double_helix=($seq,Complement($seq));
		my $i;
		my $het_hom;
		if (defined $Insertion)
		{
			$i = int rand(2);
			$het_hom = int rand(2);
			my ($o,$s1,$s2) = ($i==0) ? (1,"+","-") : (0,"-","+");
			my @ins_array=read_num($seq_length,$Insertion);
			make_Insertion(\$double_helix[$i],\@ins_array,$s1);
			make_Insertion(\$double_helix[$o],\@ins_array,$s2) if ($het_hom==1);
		}
		if (defined $Inversion)
		{
			$i = int rand(2);
			$het_hom = int rand(2);
			my ($o,$s1,$s2) = ($i==0) ? (1,"+","-") : (0,"-","+");
			my @inv_array=read_num($seq_length,$Inversion);
			make_Inversion(\$double_helix[$i],\@inv_array,$s1);
			make_Inversion(\$double_helix[$o],\@inv_array,$s2) if ($het_hom==1);
		}
		if (defined $Duplication)
		{
			$i = int rand(2);
			$het_hom = int rand(2);
			my ($o,$s1,$s2) = ($i==0) ? (1,"+","-") : (0,"-","+");
			my @dup_array=read_num($seq_length,$Duplication);
			make_Duplication(\$double_helix[$i],\@dup_array,$s1);
			make_Duplication(\$double_helix[$o],\@dup_array,$s2) if ($het_hom==1);
		}
		if (defined $Translocation)
		{
			$i = int rand(2);
			$het_hom = int rand(2);
			my ($o,$s1,$s2) = ($i==0) ? (1,"+","-") : (0,"-","+");
			my @tra_array=read_num($seq_length,$Translocation);
			make_Translocation(\$double_helix[$i],\@tra_array,$s1);
			make_Translocation(\$double_helix[$o],\@tra_array,$s2) if ($het_hom==1);
		}
		if (defined $Deletion)
		{
			$i = int rand(2);
			$het_hom = int rand(2);
			my ($o,$s1,$s2) = ($i==0) ? (1,"+","-") : (0,"-","+");
			my @del_array=read_num($seq_length,$Deletion);
			make_Deletion(\$double_helix[$i],\@del_array,$s1);
			make_Deletion(\$double_helix[$o],\@del_array,$s2) if ($het_hom==1);
		}
		if ($Divergency>0)
		{
			my $N=1;
			my @div;
			for my $j(0..1)
			{
				my $Seq="";
				#$seq = $double_helix[$j];
				#store the divergent sequences for SNPs and Indels
				$het_hom = int rand(2);
				my $div_rate = (($j==0)||($het_hom==0)) ? $Divergency : 0.0001;
				$double_helix[$j] = (($j==1)&&($het_hom==1)) ? Complement($div[0]) : $double_helix[$j];
				@div = make_div_seq($double_helix[$j],$div_rate,$Number,$j,$SeqName) if (($j==0)||(($j==1)&&($het_hom==1)));
				foreach (@div)
				{
					my $seQ = (($j==0)||(($j==1)&&($het_hom==1))) ? $_ :  Complement($_);
					chomp $seQ;
					if (defined $Sequence)
					{
						Display_seq(\$seQ,$Display);
						$seq_name=~ s/\_\d+//;
						$seq_name=$SeqName."_".$N if ((defined $Reference)&&($Number>0));
						my $length_seQ=length($seQ);
						print STDERR "simulated testing sequence plus strand done!\tLength:$Length, actural Length:$length_seQ\n" if(($Verbose)&&($j==0));
						print STDERR "simulated testing sequence minus strand done!\tLength:$Length, actural Length:$length_seQ\n" if(($Verbose)&&($j==1));
						print OUT2 (">$seq_name\n$seQ\n") if ($j==0);
						print OUT2 (">$seq_name\_complement\n$seQ\n") if ($j==1);
					}
					$Seq .= $seQ;
					$N++;
				}
				$Seq =~ s/\s+//g;
				$double_helix[0] = $Seq if ($j==0);
				$double_helix[1] = $Seq if ($j==1);
			}
		}
		else
		{
			if (defined $Sequence)
			{
				for (my $n=0;$n<$Number;$n++)
				{
					for my $j(0..1)
					{
						my $seQ = $double_helix[$j];
						Display_seq(\$seQ,$Display);
						if ((defined $Reference)&&($Number>0))
						{
							$seq_name=~ s/\_\d+//;
							$seq_name=$SeqName."_".($n+1);
						}
						my $length_seQ=length($seQ);
						print STDERR "simulated testing sequence done!\tLength:$Length, actural Length:$length_seQ\n" if(($Verbose)&&($j==0));
						print STDERR "simulated testing sequence minus strand done!\tLength:$Length, actural Length:$length_seQ\n" if(($Verbose)&&($j==1));
						print OUT2 (">$seq_name\n$seQ\n") if ($j==0);
						print OUT2 (">$seq_name\_RC\n",reverse($seQ),"\n") if ($j==1);
					}
				}
			}
		}
		makePE($PE,$double_helix[0],$double_helix[1],$seq_name) if defined ($PE);
	}
	close OUT2 if (defined $Sequence);
}

## make a random DNA BAC sequence
## usage: make_BAC($seq,$num)
#############################################
sub make_BAC
{
	my ($seq,$num)=@_;
	open (BAC,">>$Outdir/BAC.seq") || die $!;
	my $BAC_len_tot=0;
	for (my $i=0;$i<$num;$i++)
	{
		my ($BAClen,$BACpos,$BACseq,$BACseq_rc);
		my $j=0;
		BAC_LOOP:$BAClen=$BAC_len[int(rand(@BAC_len))];
		$BACpos=int(rand(length($seq)-$BAClen));
		$BACseq=uc(substr($seq,$BACpos,$BAClen));
		if (defined $Enzyme)
		{
			last if ($j>20);
			if ($BACseq!~/$Enzyme(\w+)$Enzyme_r/)
			{
				$j++;
				goto BAC_LOOP;
			}
			
			if ($BACseq =~ /$Enzyme(\w+)/)
			{
				$BACseq=$1;
			}
			else
			{
				$BACseq="";
			}
			
			my $BAC_r = reverse($BACseq);
			
			if ($BAC_r =~ /$Enzyme_r(\w+)/)
			{
				$BAC_r=$1;
			}
			else
			{
				$BAC_r="";
			}
			
			$BACseq = reverse($BAC_r);
			
			if (length($BACseq)<$BAC_length-3*$BAC_sd)
			{
				$j++;
				goto BAC_LOOP;
			}
		}
		$BACseq_rc=reverse($BACseq);
		$BACseq_rc=~tr/ACGTacgt/TGCAtgca/;
		my $BAC_rand=int(rand(2));
		$BACseq=($BAC_rand==0)?$BACseq:$BACseq_rc;
		$BACseq.=$Vector_seq[$BAC_rand];
		Display_seq(\$BACseq);
		print BAC ">BAC_$BAC_count\_$BAC_rand\n$BACseq";
		$BAC_len_tot+=$BAClen;
		$BAC_count++;
	}
	print STDERR "Number: $num\tTotal Length: $BAC_len_tot\n" if ($Verbose);
	close BAC;
}


## make a random DNA sequence with specified length
## usage: $seq=make_seq($opts{len})
#############################################
sub make_seq
{
	my $len=shift || 300;
	my $str;
	for (my $i=0; $i<$len; $i++) {
		my $order=int rand(scalar(@base)-3);
		$str.=$base[$order];
	}
	return $str;
}

## make diversity sequences of the ancestor sequences
## according to a specified rate
## usage: my @div=make_div_seq($seq,$opts{div},$opts{num},$j,$SeqName)
#############################################
sub make_div_seq
{
	my $str=shift; #store reference sequence
	my $div=shift; #diversity rate
	my $num=shift || 1; #number of output sequences
	my $j=shift;
	my $seq_name=shift;
	my $strand = ($j==1) ? "-" : "+";
	open (MU,">>$Outdir/Mutation.txt") if (defined $Mutation);
	print STDERR "Ouput $Outdir/Mutation.txt reports SNPs & Indels\n" if ((defined $Mutation)&&(defined $Verbose));
	my @out_ap; #store output in an array

	my $len=length($str);

	my $seq="";
	for (my $i=0; $i<$num; $i++)
	{
		if ($div>0)
		{
			my @pos_array=rand_pos($len,$div);
			$seq=$str;
			for (my $j=0; $j< @pos_array; $j++)
			{
				my $pos=$pos_array[$j];
			
				my $sub2=substr($seq,$pos,1);
				
				my $forebase=$sub2;
				my $change;
			
				if (defined $ChangeB)
				{
					next unless (($sub2 eq $ChangeB)||($sub2 =~ /[$ChangeB]/)||($ChangeB =~ /$sub2/i));
					do
					{
						$change=$base[int rand(scalar(@base))];
					}
					until (($change !~ /[$ChangeB]/)&&(($ChangeB !~ /$change/))&&($change ne "-")&&($sub2 ne $change));
				}
				else
				{
					do
					{
						$change=$base[int rand(scalar(@base))];
					}
					until($sub2 ne $change);
				}
				if ($change eq "-")  ## insert a random base
				{
					$sub2 = $base[int rand(@base-2)].$base[int rand(@base-2)]; #insertion
				}
				else
				{
					$sub2 = $change;
				}
				my $changebase = ($change eq " ") ? "-" : $sub2;
				print MU ("$seq_name\t",($pos+1),"\t$strand\t$forebase\t$changebase\n") if (defined $Mutation);
				my $sub1=substr($seq,0,$pos);
				my $sub3=substr($seq,$pos+1,$len-$pos-1);
				$seq=$sub1.$sub2.$sub3;
			}
			push (@out_ap,$seq);
		}
		else
		{
			push (@out_ap,$seq);
		}
	}
	close MU;
	return @out_ap;
}

## generate a array in order to stroe randem postion
## usage: @array = rand_pos($length,$divergency)
#############################################
sub rand_pos
{
	my ($length,$divergency)=@_;
	my $pos=0;
	my %POS=();
	if ($divergency>=1)
	{
		for (my $i=0;$i<$length;$i++)
		{
			$POS{$i}=1;
		}
	}
	else
	{
		for (my $i=0;$i<int($length*$divergency+0.5);$i++)
		{
			do
			{
				$pos=int rand($length);
			}
			until(!exists $POS{$pos});
			$POS{$pos}=1;
		}
	}
	return (sort{$a<=>$b} keys %POS); 
}

##simulate paired-end reads also store the file into the "simulated-result" directory
##usage: makePE($PE,$seq1,$seq2,$Name);
#############################################
sub makePE
{
	my ($PE_num,$seq1,$seq2,$seqname)=@_;
	$seqname ||= "chrS";
	$seq1 =~ s/\-//g;
	$seq2 =~ s/\-//g;
	my $length=length($seq1);
	print STDERR "simulating sequencing reads from $seqname\nSequence Length: $length\n" if (defined $Verbose);
	
	my $multimes=1;
	if ($PE_num =~ /(\d+)X/i)
	{
		$multimes = $1;
		$PE_num = int ($multimes*$length/($RL1+$RL2)+0.5);;
	}

	my ($read1,$read2,$qual1,$qual2);
	nameRead(\$read1,\$read2,\$qual1,\$qual2);
	
	open OUT3,">>$Outdir/$read1" or die "Can not create file:$read1\n" unless (($RL1==0)||($RL1 eq ""));
	open OUT4,">>$Outdir/$read2" or die "Can not create file:$read2\n" unless (($RL2==0)||($RL2 eq ""));
	open OUT5,">>$Outdir/$qual1" or die "Can not create file:$qual1\n" if (($RL1>0)&&($Type !~ /std/i)&&($Type=~/fq/i));
	open OUT6,">>$Outdir/$qual2" or die "Can not create file:$qual2\n" if (($RL2>0)&&($Type !~ /std/i)&&($Type=~/fq/i));
	if (($Type =~ /Solid/i)&&($Type !~ /std/i)&&($print_title==0))
	{
		my $title = ($Name eq "null") ? "BENM_2011_0526_Solid_" : $Name."_";
		print OUT3 "# Title: $title\n" unless (($RL1==0)||($RL1 eq ""));
		print OUT4 "# Title: $title\n" unless (($RL2==0)||($RL2 eq ""));
		print OUT5 "# Title: $title\n" if (($RL1>0)&&($Type !~ /std/i)&&($Type=~/fq/i));
		print OUT6 "# Title: $title\n" if (($RL1>0)&&($Type !~ /std/i)&&($Type=~/fq/i));
		$print_title=1;
	}
	
	my ($lane,$loop,$insert)=(0,0,0);
	my @Insert = gauss($PE_num,$InsertSizePE,$InsertSizeSD,$multimes);
	my @ReadLength1 = gauss($PE_num,$RL1,$ReadsSD,$multimes) if (($ReadsSD != 0)&&(($Type =~ /454/)||($Type =~ /3730/))&&(defined $RL1)&&($RL1 != 0));
	my @ReadLength2 = gauss($PE_num,$RL2,$ReadsSD,$multimes) if (($ReadsSD != 0)&&(($Type =~ /454/)||($Type =~ /3730/))&&(defined $RL2)&&($RL2 != 0));
	
	for (my $n=0;$n<$PE_num;$n++)
	{
		#for simulating heterozygous variation
		my $H=int rand (2);
		my $seq_h = ($H==0) ? $seq1: $seq2;
		$insert = $Insert[$n];
		if (($Type =~ /454/)||($Type =~ /3730/))
		{
			if ($ReadsSD != 0)
			{
				$RL1=abs($ReadLength1[$n]) if ((defined $RL1)&&($RL1 != 0));
				$RL2=abs($ReadLength2[$n]) if ((defined $RL2)&&($RL2 != 0));
			}
		}
		$lane = int (rand(8)) + 1;
		$loop++;
		my $odd_len = ($insert-$RL1>0) ? $insert : (($RL1>0)?$RL1:$RL2);

		RLOOP:my $pos= (defined $Circle) ? int rand(length($seq_h)) : int rand(length($seq_h)-$odd_len);
		my ($PE1,$PE2) = ("","");

		if ( (defined $Circle) && ($pos>length($seq_h)-$odd_len) )
		{
			if ($RL1>0)
			{
				if ($pos>length($seq_h)-$RL1)
				{
					my $part1=substr($seq_h,$pos,length($seq_h)-$pos);
					my $part2=substr($seq_h,0,$RL1-(length($seq_h)-$pos));
					$PE1="$part1$part2";
				}
				else
				{
					$PE1=substr($seq_h,$pos,$RL1);
				}
				die ("$PE1\t$pos\n") if (length($PE1)!=$RL1);
			}
			if ($RL2>0)
			{
				if ((length($seq_h)-($pos+$insert-$RL2)>0)&&(length($seq_h)-($pos+$insert-$RL2)<$RL2))
				{
					my $PE2_pos=$pos+$insert-$RL2;
					my $part1=substr($seq_h,$PE2_pos,length($seq_h)-($pos+$insert-$RL2));
					my $part2=substr($seq_h,0,$RL2-(length($seq_h)-($pos+$insert-$RL2)));
					$PE2="$part1$part2";
					die ("$PE2\t$PE2_pos\n") if (length($PE2)!=$RL2);
				}
				else
				{
					my $PE2_pos=($insert-(length($seq_h)-$pos)>0)?$insert-(length($seq_h)-$pos):0;
					$PE2=substr($seq_h,$PE2_pos,$RL2);
					die ("$PE2\t$PE2_pos\n") if (length($PE2)!=$RL2);
				}
			}
		}
		else
		{
			$PE1=substr($seq_h,$pos,$RL1) if ($RL1 > 0);
			$PE2=substr($seq_h,$pos+$insert-$RL2,$RL2) if ($RL2 > 0);
		}
		if ($H==1)
		{
			$PE1=Complement($PE1);
			$PE2=Complement($PE2);
		}
		goto RLOOP if (((defined $PE1)&&($PE1 =~ /N{3}/))||((defined $PE2)&&($PE2 =~ /N{3}/)));
		
		if (($Type =~ /solid/i)&&((($RL1>0)&&($PE1=~/X/i))||(($RL1>0)&&($PE2=~/X/i))))
		{
			$PE1 =~ s/X/N/ig if ($RL1>0);
			$PE2 =~ s/X/N/ig if ($RL2>0);
		}
		
		if (defined $Adapter)
		{
			my ($l,$r,$l_len,$r_len,$l_s,$r_s);
			if (@adapter1>0)
			{
				$l = int(rand(@adapter1));
				$l_len = $rand_len[int rand(@rand_len)];
				$l_len = (abs($l_len)>length($adapter1[$l])) ? length($adapter1[$l]) : abs($l_len);
				$l_s = substr($adapter1[$l],length($adapter1[$l])-$l_len,$l_len);
			}
			if (@adapter2>0)
			{
				$r = int(rand(@adapter2));
				$r_len = $rand_len[int rand(@rand_len)];
				$r_len = (abs($r_len)>length($adapter2[$r])) ? length($adapter2[$r]) : abs($r_len);
				$r_s = substr($adapter2[$r],0,$r_len);
			}
			if (($RL1>0)&&($RL2==0))
			{
				substr($PE1,0,$l_len,$l_s) if (($l_len>0)&&($RL1>$l_len));
				substr($PE1,$RL1-$r_len,$r_len,$r_s) if (($r_len>0)&&($RL1>$r_len));
			}
			elsif (($RL2>0)&&($RL1==0))
			{
				substr($PE2,0,$l_len,$l_s) if (($l_len>0)&&($RL2>$l_len));
				substr($PE2,$RL2-$r_len,$r_len,$r_s) if (($r_len>0)&&($RL2>$r_len));
			}
			else
			{
				substr($PE1,0,$l_len,$l_s) if ((defined $l_len)&&($RL1>$l_len));
				substr($PE2,$RL2-$r_len,$r_len,$r_s) if ((defined $r_len)&&($RL2>$r_len));
			}
		}
		
		($PE1,$PE2) = Relative_Orientation($PE1,$PE2);

		##simulate pair-end reads quality
		my $quality1 = "";
		my $quality2 = "";
		#$RL1=length($PE1) if ($PE1 ne "");
		#$RL2=length($PE2) if ($PE2 ne "");
		generate_qual(\$quality1,\$quality2,$RL1,$RL2) if ($Type =~ /fq/);
		
		#$PE1 =~ tr/ACGTNacgtn/0123.0123./ if (($Type =~ /solid/i)&&(defined $PE1)&&($RL1>0));
		#$PE2 =~ tr/ACGTNacgtn/0123.0123./ if (($Type =~ /solid/i)&&(defined $PE2)&&($RL2>0));
		base2col(\$PE1) if (($Type =~ /solid/i)&&(defined $PE1)&&($RL1>0));
		base2col(\$PE2) if (($Type =~ /solid/i)&&(defined $PE2)&&($RL2>0));
		
		my $print_name = ($Name eq "null") ? $seqname : $Name;
		
## standart Solexa paired-end read name format: @HWI-EAS6_2_FC12994_PE:8:1:678:777
## SOLiD Name format: >3_16_150_F3
## 454 Name format: @SRR000921.7  E7OAHBI01A9IFH  length=173
## 3730 Name format : read eeq03a02.b1 is univ fwd template: eeq03a02 library: eeq03
##              or : mgsaea0_000103.z1.scf CHROMAT_FILE: mgsaea0_000103.z1.scf PHD_FILE: mgsaea0_000103.z1.scf.phd.1 CHEM: unknown DYE: unknown TIME: Wed Dec 29 09:56:26 2004
		my ($PE1_name,$PE2_name);
		if (($Type =~ /Solexa/i)||($Type =~ /^std/i))
		{
			$PE1_name =$print_name.":$lane:$loop:$insert:$RL1\/1";
			$PE2_name =$print_name.":$lane:$loop:$insert:$RL2\/2";
		}
		elsif ($Type =~ /Solid/i)
		{
			$PE1_name = "$lane\_$loop\_$insert\_F3";
			$PE2_name = "$lane\_$loop\_$insert\_R3";
		}
		elsif ($Type =~ /454/)
		{
			$PE1_name = "$print_name"."F template=$print_name dir=F $Type\_$pos length=$RL1";
			$PE2_name = "$print_name"."R template=$print_name dir=R $Type\_$pos length=$RL2";
			if (($RL1>0)&&($RL2>0)) 
			{
				$PE1_name.=" library=pairlab";
				$PE2_name.=" library=pairlab";
			}
		}
		elsif ($Type =~ /3730/)
		{
			$PE1_name = "$print_name\_$Type$pos.y1.scf fwd template $print_name length=$RL1";
			$PE2_name = "$print_name\_$Type$pos.z1.scf rev template $print_name length=$RL2";
			if (($RL1>0)&&($RL2>0)) 
			{
				$PE1_name.=" library=pairlab";
				$PE2_name.=" library=pairlab";
			}
		}
		my $print_mark=($Type=~/std/i)? "@" : ">";
		Display_seq(\$PE1,$Display) if (($PE1 ne "")&&($Type !~ /std/i));
		Display_seq(\$PE2,$Display) if (($PE2 ne "")&&($Type !~ /std/i));
		if ($RL1>0)
		{
			print OUT3 ("$print_mark$PE1_name\n");
			print OUT3 "$PE1";
			print OUT3 "\n" if ($Type =~ /std/);
		}
		if ($RL2>0)
		{
			print OUT4 ("$print_mark$PE2_name\n");
			print OUT4 "$PE2";
			print OUT4 "\n" if ($Type =~ /std/);
		}
		if ($Type =~ /std/i)
		{
			if ($PE1 ne "")
			{
				print OUT3 ("+\n$quality1\n");
			}
			if ($PE2 ne "")
			{
				print OUT4 ("+\n$quality2\n");
			}
		}
		elsif ($Type =~ /fq/i)
		{
			if (($RL1>0)&&($PE1 ne ""))
			{
				if ($Type=~/solexa/i)
				{
					Display_seq(\$quality1,$Display);
					print OUT5 (">$PE1_name\n$quality1");
				}
				else
				{
					print OUT5 (">$PE1_name\n$quality1\n");
				}
			}
			if (($RL2>0)&&($PE2 ne ""))
			{
				if ($Type=~/solexa/i)
				{
					Display_seq(\$quality2,$Display);
					print OUT6 (">$PE2_name\n$quality2");
				}
				else
				{
					print OUT6 (">$PE2_name\n$quality2\n");
				}
			}
		}
	}

	my $PE1_bases=$PE_num*$RL1;
	my $PE2_bases=$PE_num*$RL2;
	my $PE_reads=$PE_num*2;
	print STDERR "Simulated paired-end reads done! Total bases:$PE1_bases+$PE2_bases=".($PE1_bases+$PE2_bases) if($Verbose);
	print STDERR " Total reads:$PE_num+$PE_num=".($PE_reads)."\n" if($Verbose);
	print STDERR "Real coverage: ".(($PE1_bases+$PE2_bases)/$length)."\n" if($Verbose);
	close OUT3;
	close OUT4;
	close OUT5;
	close OUT6;
}

## name reads
## usage: nameRead(\$read1,\$read2,\$qual1,\$qual2)
#############################################
sub nameRead
{
	my ($read1,$read2,$qual1,$qual2)=@_;
	if (($Type =~ /Solexa/i)||($Type =~ /^f[aq]$/i)||($Type =~ /^fast[aq]$/i))
	{
#s_3_0023_sequence.txt
#s_3_0023_qcal.txt
		if ($Type !~ /std/)
		{
			$$read1 = ($Name eq "null") ? "s_1_sequence.txt" : $Name."_1_sequence.txt";
			$$read2 = ($Name eq "null") ? "s_2_sequence.txt" : $Name."_2_sequence.txt";
			$$qual1 = ($Name eq "null") ? "s_1_qcal.txt" : $Name."_1_qcal.txt";
			$$qual2 = ($Name eq "null") ? "s_2_qcal.txt" : $Name."_2_qcal.txt";
		}
		else
		{
			$$read1 = ($Name eq "null") ? "s1.fastq" : $Name."_s1.fastq";
			$$read2 = ($Name eq "null") ? "s2.fastq" : $Name."_s2.fastq";
		}
	}
	elsif ($Type =~ /Solid/i)
	{
#BARB_20071114_2_YorubanMP-BC3_F3.csfasta
#BARB_20071114_2_YorubanMP-BC3_F3_QV.qual
		if ($Type !~ /std/)
		{
			$$read1 = ($Name eq "null") ? "F3.csfasta" : $Name."_F3.csfasta";
			$$read2 = ($Name eq "null") ? "R3.csfasta" : $Name."_R3.csfasta";
			$$qual1 = ($Name eq "null") ? "F3_QV.qual" : $Name."_F3_QV.qual";
			$$qual2 = ($Name eq "null") ? "R3_QV.qual" : $Name."_R3_QV.qual";
		}
		else
		{
			$$read1 = ($Name eq "null") ? "F3.fastq" : $Name."_F3.fastq";
			$$read2 = ($Name eq "null") ? "R3.fastq" : $Name."_R3.fastq";
		}
	}
	elsif ($Type =~ /454/i)
	{
#1.TCA.454Reads.fna
#1.TCA.454Reads.qual
		if ($Type !~ /std/)
		{
			$$read1 = ($Name eq "null") ? "1.454Reads.fna" : $Name."_1.454Reads.fna";
			$$read2 = ($Name eq "null") ? "2.454Reads.fna" : $Name."_2.454Reads.fna";
			$$qual1 = ($Name eq "null") ? "1.454Reads.qual" : $Name."_1.454Reads.qual";
			$$qual2 = ($Name eq "null") ? "2.454Reads.qual" : $Name."_2.454Reads.qual";
		}
		else
		{
			$$read1 = ($Name eq "null") ? "1.454Reads.fastq" : $Name."_1.454Reads.fastq";
			$$read2 = ($Name eq "null") ? "2.454Reads.fastq" : $Name."_2.454Reads.fastq";
		}
	}
	elsif ($Type =~ /3730/i)
	{
#Y34.sanger.reads.seq
#Y34.sanger.reads.seq.qual
		if ($Type !~ /std/)
		{
			$$read1 = ($Name eq "null") ? "1.sanger.reads.seq" : $Name."_1.sanger.reads.seq";
			$$read2 = ($Name eq "null") ? "2.sanger.reads.seq" : $Name."_2.sanger.reads.seq";
			$$qual1 = ($Name eq "null") ? "1.sanger.reads.seq.qual" : $Name."_1.sanger.reads.seq.qual";
			$$qual2 = ($Name eq "null") ? "2.sanger.reads.seq.qual" : $Name."_2.sanger.reads.seq.qual";
		}
		else
		{
			$$read1 = ($Name eq "null") ? "1.sanger.fastq" : $Name."_1.sanger.fastq";
			$$read2 = ($Name eq "null") ? "2.sanger.fastq" : $Name."_2.sanger.fastq";
		}
	}
	elsif ($Type =~ /^stdfq$/i)
	{
		$$read1 = ($Name eq "null") ? "PE1.fq" : $Name."_PE1.fq";
		$$read2 = ($Name eq "null") ? "PE2.fq" : $Name."_PE2.fq";
	}
}

## make distributing for PE insert size around the set standard deviation
## usage: @array = gauss($num,$mean,$sd)
#############################################
sub gauss
{
	my ($num,$mean,$sd,$multimes) = @_;
	$multimes ||= 1;
	my $rate1;
	my $rate2;
	my $random1;
	my $random2;
	my @random;
	my $below = $mean-3*$sd;
	my $above = $mean+3*$sd;
	for (my $i=0;$i<2.71828183*$num/$multimes;$i++)
	{
		$rate1 = rand(1);
		$rate2 = rand(1);
		$random1 = sqrt(-2*log($rate1))*sin(2*PI*$rate2);	#simulate Gauss distribution
		$random2 = sqrt(-2*log($rate2))*cos(2*PI*$rate1);
		$random1 = int ($mean+($random1*$sd));
		$random2 = int ($mean+($random2*$sd));
		($random1>$below) && ($random1<$above) && push @random,$random1;
		($random2>$below) && ($random2<$above) && push @random,$random2;
	}
	my @Random;
	for (my $j=0;$j<$multimes;$j++)
	{
		push @Random,@random;
	}
	return @Random;
}

## make distributing for PE insert size around the set standard deviation
## usage: ($PE1,$PE2)=Relative_Orientation($PE1,$PE2)
#############################################
sub Relative_Orientation
{
	my ($R1,$R2) = @_;
	if ($NormalStrand eq "FR")
	{
		$R2 = reverse($R2);
		$R2 = Complement($R2);
	}
	elsif ($NormalStrand eq "RF")
	{
		$R1 = reverse($R1);
		$R1 = Complement($R1);
	}
	elsif ($NormalStrand eq "RR")
	{
		$R1 = reverse($R1);
		$R1 = Complement($R1);
		$R2 = reverse($R2);
		$R2 = Complement($R2);
	}
	return ($R1,$R2);
}

## make Complement sequence
## usage: $complemental_seq=Complement($seq)
#############################################
sub Complement
{
	my $Temp_seq = shift;
	my $Comp_seq = $Temp_seq;
	$Comp_seq =~ tr/ACGTacgt/TGCAtgca/;
	return $Comp_seq;
}

## convert bases to color int reads
## usage: base2col(\$seq)
#############################################
sub base2col
{
	my $reads = shift;
	my @bases = split '',uc($$reads);
	my $col_code = (defined $Header) ? $Header : "G";
	my $current_base = $col_code;
	my $code = '';
	for(my $i=0;$i<@bases;$i++)
	{
		$code = (exists $colcode{"$current_base$bases[$i]"}) ? $colcode{"$current_base$bases[$i]"} : $othercode[int(rand(@othercode))];
		$col_code .= $code;
		$current_base = $bases[$i];
	}
	$$reads=$col_code;
}

## convert color int reads to bases reads
## usage: col2base(\$seq)
#############################################
sub col2base
{
	my $reads = shift;
	my $col = $reads;
	my @colors = split '',$col;
	my $string = $colors[0];
	if($string !~ /[acgt]/i)
	{
		warn "$col has no header base\n";
		return 0;
	}
	my $last_base = $string;
	my $current_base = '';
	for(my $i=1;$i<@colors;$i++)
	{
		if (($last_base=~/N/i)&&($colors[$i]==5))
		{
			$current_base = $bases[int(rand(@bases))];
		}
		else
		{
			$current_base = (exists $decode{$colors[$i]}->{$last_base}) ? $decode{$colors[$i]}->{$last_base} : "N";
		}
		$string .= $current_base;
		$last_base = $current_base;
	}
	$$reads = $string;
}

#
#
########
sub generate_qual
{
	my ($quality1,$quality2,$RL1,$RL2)=@_;
	
	#convert Dec code to ASCII code for random quality value (solexa or std/standard quality)
	if ((($Type =~ /\d+/)||($Type =~ /solid/i))&&($Type !~ /std/i))
	{
		my $q1=int(rand(40-$qual_th[1])+$qual_th[1]+1);
		my $q2=int(rand(40-$qual_th[2])+$qual_th[2]+1);
		if ($RL1>$qual_th[0])
		{
			for (my $i=0;$i<$RL1-$qual_th[0];$i++)
			{
				$$quality1 .= sprintf("%02d ",int(rand($qual_th[3])+1+$q1));
				$$quality1 =~ s/\s+$/\n/ if (($i+1)%$Display==0);
			}
			for (my $i=$RL1-$qual_th[0];$i<$RL1;$i++)
			{
				$$quality1 .= sprintf("%02d ",int(rand($qual_th[3])+1+$q2));
				$$quality1 =~ s/\s+$/\n/ if (($i+1)%$Display==0);
			}
		}
		if ($RL2>$qual_th[0])
		{
			for (my $i=0;$i<$RL2-$qual_th[0];$i++)
			{
				$$quality2 .= sprintf("%02d ",int(rand($qual_th[3])+1+$q1));
				$$quality2 =~ s/\s+$/\n/ if (($i+1)%$Display==0);
			}
			for (my $i=$RL2-$qual_th[0];$i<$RL2;$i++)
			{
				$$quality2 .= sprintf("%02d ",int(rand($qual_th[3])+1+$q2));
				$$quality2 =~ s/\s+$/\n/ if (($i+1)%$Display==0);
			}
		}
	}
	else
	{
		if ($RL1>$qual_th[0])
		{
			my $q1=int(rand(40-$qual_th[1])+$qual_th[1]+1);
			my $q2=int(rand(40-$qual_th[2])+$qual_th[2]+1);
			for (my $i=0;$i<$RL1-$qual_th[0];$i++)
			{
				if ($Type =~ /std/i)
				{
					$$quality1 .= $conv_table[int(rand($qual_th[3])+$q1)+64]; #std quality
				}
				elsif ($Type =~ /Solexa/i)
				{
					$$quality1 .= chr(int(rand($qual_th[3]))+$q1+64); #solexa quality
				}
			}
			for (my $i=$RL1-$qual_th[0];$i<$RL1;$i++)
			{
				if ($Type =~ /std/i)
				{
					$$quality1 .= $conv_table[int(rand($qual_th[3])+$q2)+64]; #std quality
				}
				elsif ($Type =~ /Solexa/i)
				{
					$$quality1 .= chr(int(rand($qual_th[3]))+$q2+64); #solexa quality
				}
			}
		}
		if ($RL2>$qual_th[0])
		{
			my $q1=int(rand(40-$qual_th[1])+$qual_th[1]+1);
			my $q2=int(rand(40-$qual_th[2])+$qual_th[2]+1);
			for (my $j=0;$j<$RL2-$qual_th[0];$j++)
			{
				if ($Type =~ /std/i)
				{
					$$quality2 .= $conv_table[int(rand($qual_th[3])+$q1)+64]; #std quality
				}
				elsif ($Type =~ /Solexa/i)
				{
					$$quality2 .= chr(int(rand($qual_th[3]))+$q1+64); #solexa quality
				}
			}
			for (my $j=$RL2-$qual_th[0];$j<$RL2;$j++)
			{
				if ($Type =~ /std/i)
				{
					$$quality2 .= $conv_table[int(rand($qual_th[3])+$q2)+64]; #std quality
				}
				elsif ($Type =~ /Solexa/i)
				{
					$$quality2 .= chr(int(rand($qual_th[3]))+$q2+64); #solexa quality
				}
			}
		}
	}
}

#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq
{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s$//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line)
	{
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}

#read the string of setting for the number and size of variation, and return an array
#This array has four elements: (number, start position, end position, length)
#usage: @num_array=read_num($seq_length,$string);
#############################################
sub read_num
{
	my ($seq_length,$string)=@_;
	chomp $string;
	my @NUM=();
	my @array=($string=~/\d+\-\w+/) ? (split /-/,$string) : ($string,0,0,0);
	if ($array[0] =~ /^\d+/) 
	{
		$NUM[0]=int($array[0]+0.5);
		$NUM[3]=($string=~/\-size:(\d+)$/) ? int($1+0.5) : 100+int(rand(1000));
		$NUM[3]=(defined $Size) ? ( ($Size=~/\~(\d+)/) && (defined $InsertSizeSD) ) ? int ($1+rand($InsertSizeSD)+0.5) : $Size : $NUM[3];
		$NUM[3]=($string=~/\-start:(\d+)\-end:(\d+)/) ? (($2-$1)>0) ? int($2-$1+1.5) : die `pod2text $0` : $NUM[3];
		if ($string=~/\-start:(\d+)/)
		{
			$NUM[1]= int($1-0.5) ;
		}
		else
		{
			$NUM[1]= ($string=~/\-end:(\d+)/ ) ? int($1-$NUM[3]+0.5) : int(rand($seq_length-2000));
		}
		if ($string=~/\-end:(\d+)/)
		{ 
			$NUM[2]= int($1-0.5); 
		}
		else
		{
			$NUM[2]= ($string=~/\-start:(\d+)/) ? int($1+$NUM[3]-0.5) : int($NUM[1]+$NUM[3]-0.5);
		}
	}
	else
	{
		die `pod2text $0`;
	}
	return @NUM;
}

#make Deletion variation and output the position
#usage: make_Deletion(\$string,\@num);
#############################################
sub make_Deletion
{
	my $seq_d=shift;
	my $num_d= shift;
	my $hel_d=shift;
	my $del;
	my $DelPos;
	my $DelLen;
	
	my ($num,$start,$end,$size)= @$num_d;
	
	$$seq_d =~ s/\s//g;
	for (my $x=0;$x<$num;$x++)
	{
		$DelLen = $size;
		if ($x>1)
		{
			$DelLen = (defined $Size) ? ( ($Size=~/\~(\d+)/) && (defined $InsertSizeSD) ) ? int ($1+rand($InsertSizeSD)+0.5) : $size : 100+int(rand(1000));
		}
		$DelPos = $start*($x+1);
		($DelPos = ( int (rand (length($$seq_d) )/$num) ) ) if (($DelPos<0) || $DelPos>length($$seq_d)-$DelLen-1);
		$del = substr($$seq_d,$DelPos,$DelLen);
		print "Deletion_$num $hel_d ".($x+1).": ".($DelPos+1)."\t".($DelPos+$DelLen)."\tLength:$DelLen\n$del\n";
		substr($$seq_d,$DelPos,$DelLen,"");
	}
}

#make Insertion variation and output the position
#usage: make_Insertion(\$string,\@num);
#############################################
sub make_Insertion
{
	my $seq_s=shift;
	my $num_s=shift;
	my $hel_s=shift;
	my $ins;
	my $InsPos;
	my $InsLen;
	
	my ($num,$start,$end,$size)=@$num_s;
	
	$$seq_s =~ s/\s//g;
	for (my $y=0;$y<$num;$y++)
	{
		$InsLen = $size;
		if ($y>1)
		{
			$InsLen = (defined $Size) ? ( ($Size=~/\~(\d+)/) && (defined $InsertSizeSD) ) ? int ($1+rand($InsertSizeSD)+0.5) : $size : 100+int(rand(1000));
		}
		$InsPos = $start*($y+1);
		($InsPos= ( int (rand (length($$seq_s) )/$num) ) ) if (($InsPos<0) || $InsPos>length($$seq_s)-$InsLen-1);
		$ins = make_seq($InsLen);
		substr($$seq_s,$InsPos,0,$ins);
		print "Insertion_$num $hel_s ".($y+1).": ".($InsPos+1)."\tLength:$InsLen\n$ins\n";
	}
}

#make Inversion variation and output the position
#usage: make_Inversion(\$string,\@num);
#############################################
sub make_Inversion
{
	my $seq_v=shift;
	my $num_v=shift;
	my $hel_v=shift;
	my $inv;
	my $InvPos;
	my $InvLen;
	
	my ($num,$start,$end,$size)=@$num_v;
	
	$$seq_v =~ s/\s//g;
	
	for (my $z=0;$z<$num;$z++)
	{
		$InvLen = $size;
		if ($z>1)
		{
			$InvLen = (defined $Size) ? ( ($Size=~/\~(\d+)/) && (defined $InsertSizeSD) ) ? int ($1+rand($InsertSizeSD)+0.5) : $size : 100+int(rand(1000));
		}
		$InvPos = $start*($z+1);
		($InvPos= ( int (rand (length($$seq_v) )/$num) ) ) if (($InvPos<0) || $InvPos>length($$seq_v)-$InvLen-1);
		$inv = substr($$seq_v,$InvPos,$InvLen);
		my $inv=reverse($inv);
		$inv = Complement($inv);
		substr($$seq_v,$InvPos,$InvLen,$inv);
		print "Inversion_$num $hel_v ".($z+1).": ".($InvPos+1)."\t".($InvPos+$InvLen)."\tLength:$InvLen\n$inv\n";
	}
}

#make Duplication variation and output the position where it come from and has been duplicate and locate to where
#usage: make_Duplication(\$string,\@num);
#############################################
sub make_Duplication
{
	my $seq_dup=shift;
	my $num_dup=shift;
	my $hel_dup=shift;
	my $duplication;
	my $DupPos1;
	my $DupPos2;
	my $DupLen;
	
	my ($num,$start,$end,$size)=@$num_dup;
	
	$$seq_dup =~ s/\s//g;
	for (my $xx=0;$xx<$num;$xx++)
	{
		$DupLen = $size;
		if ($xx>1)
		{
			$DupLen = (defined $Size) ? ( ($Size=~/\~(\d+)/) && (defined $InsertSizeSD) ) ? int ($1+rand($InsertSizeSD)+0.5) : $size : 100+int(rand(1000));
		}
		$DupPos1 = $start*($xx+1);
		$DupPos1 = ( int ( (rand (length($$seq_dup) )/$num) ) ) if (( $DupPos1<0) || ($DupPos1>length($$seq_dup)-$DupLen-1) );
		my $td_rand=int rand(2);
		if ($td_rand==1)
		{
			$DupPos2 = int (rand (length($$seq_dup)-($size+1000)*$num));
		}
		else
		{
			my $ud_rand=int rand(2);
			if ($ud_rand==1)
			{
			$DupPos2 = $DupPos1 + $size + int(rand(100*$xx));
			}
			else
			{
			$DupPos2 = $DupPos1 - $size - int(rand(100*$xx));
			}
		}
		$DupPos2 = ( int ( (rand (length($$seq_dup) ) )/($num+1)) ) if (($DupPos2<0) || ($DupPos2==$DupPos1) || ($DupPos2>length($$seq_dup)-$DupLen-1));
		$duplication = substr($$seq_dup,$DupPos1,$DupLen);
		substr($$seq_dup,$DupPos2,0,$duplication);
		print "Duplcation_$num $hel_dup ".($xx+1).": From: $DupPos1\t".($DupPos1+$DupLen-1)."\t";
		print "To: ".($DupPos2+1)."\t".($DupPos2+$DupLen)."\tLength:$DupLen\n$duplication\n";
	}
}

#make Translocation variation and output the position where it was from and has been translocated to where
#it equal to Deletion+Invsertion
#usage: make_Translocation(\$string,\@num);
#############################################
sub make_Translocation
{
	my $seq_t=shift;
	my $num_t=shift;
	my $hel_t=shift;
	my $translocation;
	my $TranPos1;
	my $TranPos2;
	my $TranLen;
	
	my ($num,$start,$end,$size)=@$num_t;
	
	$$seq_t =~ s/\s//g;
	for (my $yy=0;$yy<$num;$yy++)
	{
		$TranLen = $size;
		if ($yy>1)
		{
			$TranLen = (defined $Size) ? ( ($Size=~/\~(\d+)/) && (defined $InsertSizeSD) ) ? int ($1+rand($InsertSizeSD)+0.5) : $size : 100+int(rand(1000));
		}
		$TranPos1 = $start*($yy+1);
		$TranPos1 = ( int ( (rand (length($$seq_t) )/$num)) ) if (($TranPos1<0) || ($TranPos1>length($$seq_t)-$TranLen-1));
		$TranPos2 = int (rand (length($$seq_t)-($TranLen+1000)*$num));
		$TranPos2 = ( int ( (rand (length($$seq_t) ) )/($num+1)) ) if ( ($TranPos2<0) || ($TranPos2>length($$seq_t)-$TranLen-1) || ($TranPos2==$TranPos1) || ($TranPos1>length($$seq_t)-$TranLen-1));;
		$translocation = substr($$seq_t,$TranPos1,$TranLen);
		substr($$seq_t,$TranPos1,$TranLen,"");		#deletion place
		substr($$seq_t,$TranPos2,0,$translocation);	#insertion place
		print "Translocation_$num $hel_t ".($yy+1).": From: $TranPos1\t".($TranPos1+$TranLen-1)."\t";
		print "To: ".($TranPos2+1)."\t".($TranPos2+$TranLen)."\tLength:$TranLen\n$translocation\n";
	}
}

######################################## Bottom ########################################