#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my $usage = qq(
Transformers -- biological data format transforming tools.
Author : BENM <BinxiaoFeng\@gmail.com>
Version: 0.2.10 alpha 2013-03-11
Birthday : 2009-08-09 0.1.0
Usage: Transformers <command> [Option]
Command:
##Sequence transformer
         fa2std		Convert FASTA to the standard(phred64-quals) FASTQ
         fq2faq		Convert various FASTQ-like format to FASTA and spliting out QUAL
         fa2cs		Convert FASTA to SOLiD/ABi CSFASTA
         cs2fa		Convert SOLiD/ABi CSFASTA to FASTA
         seqal2fq	Combine Seq and Qual files to FASTQ file
         qseq2fq	Covever Solexa/Illumina qseq format to FASTQ format
         
         sol2std	Convert Solexa/Illumina (ASCII offset value of 64) format to the standard format
         fqint2std	Convert FASTQ-int format (3730,454,SOLiD) to the standard/Sanger(phred64-quals) format
         scarf2std	Convert SCARF format to the standard/Sanger(phred64-quals) format

         std2sol	Convert standard/Sanger(phred64-quals) FASTQ format to Solexa format
         std2solfq	Convert standard/Sanger(phred64-quals) FASTQ format to Solexa FASTQ format
         std2fqint	Convert standard/Sanger(phred64-quals) FASTQ format to FASTQ-int format (3730,454,SOLiD)
         std2scarf	Convert standard/Sanger(phred64-quals) FASTQ format to SCARF format

         solid2fq	Convert SOLiD/ABi color space FASTQ like data to FASTQ file
         fq2solid	Convert FASTQ file to SOLiD/ABi color space FASTQ like data, relative to solid2fq

         sam2fq		Convert Sequence Alignment/Map(SAM) format to FASTQ format
         sam2pe		Convert Sequence Alignment/Map(SAM) format to PE FASTQ format

         seq4newbler	Change sequencing reads for NEWBLER assembling
         seq4phrap	Change sequencing reads for PHRAP assembling
         gs2pe		Covert 454/Roche Genome Sequencer FLX System Paired End reads to standard PE

         shuffle	shuffle two PE sequences to one PE sequence for VELVET
         distribute	distribute one shuffled PE sequence to two distributed PE sequences

         mergepe	merge PE overlap

         cds2aa		translate nucleic acid to amino acid

         sort2pe	sort Solexa/Illumina reads to paired reads

##Assembly/Alignment format transformer:
psl2m8, m82psl, maf2cns, maqview, maq2ace, maq2sam, sam2maq, soap2sam, sam2soap, ace2afg, afg2bank, maf2sam, sam2bam, sam2bed

##Var annotated file format transformer:
sam2vcf, bam2vcf, cns2vcf, vcf2gff, vcf2bed, vcf2fq

#Annotaion file format transformer:
gtf2gff, wig2gff, bed2gff, psl2gff; gff2gtf, wig2gtf, bed2gtf, psl2gtf;
gff2wig, gtf2wig, bed2wig, psl2wig; gff2bed, gtf2bed, wig2bed, psl2bed;
gff2psl, gtf2psl, bed2psl, wig2psl

         instruction	Explanation to different format
         example	Show examples of various formats

Option:  -S <str>		seq file
         -Q <str>		qual file
         -FQ <str>		FASTQ file
         -p <str>		output file's prefix name
         -f <str>		output format for std2fqint command [Solexa|454|3730|SOLiD]
         -qf <str>		ASCII quality format [SXIJL] follow FASTQ quality instruction, default: L
         -t <str>		set typeset
         -header <str>		SOLiD header base, for colorspace coverting, default: G
         -RC			reversed and complment sequence
         -eh			SOLiD header base is existing in bases reads, for command "fa2cs", "std2fqint" && "fq2solid"
         -mp <str>		input mate pairs(paired ends) sequence according 'shuffle' process
         -linker <str>		input Titanium linkers sequence (Circularization Adaptor sequence)
         -pelib			pairlab trnaforming
         -pe_strand <str>	set two strand of paired-end(or mate-pair) [FR|RF|FF|RR], default: FR
         -author		Author & Copyright information
         -version		show update info
         -help			Help Information

Example:
 (1) Transformers fqint2std -S BARB_20071114_2_YorubanMP-BC3_F3.csfasta -Q BARB_20071114_2_YorubanMP-BC3_F3_QV.qual -p BARB_20071114_2_YorubanMP-BC3_F3
 (2) Transformers std2fqint -FQ BARB_20071114_2_YorubanMP-BC3_F3.fq -header T -p BARB_20071114_2_YorubanMP-BC3_F3 -f SOLiD
\n);

my $update_log=qq(
Transformers -- NGS sequence data format transforming tools.
Latest Version: 0.2.5 alpha 
Update1: 2009-08-12 0.1.1 Add TRASFORMERS LOGO and amend some bugs of \'fqint2std\' & \'std2scarf\'
Update2: 2009-08-14 0.1.2 Add \'solid2fq\' and \'fq2solid\' function for color space FASTQ like data transforming
Update3: 2009-08-24 0.1.3 Add \'-RC\' option
Update4: 2009-10-06 0.1.4 Modify SOLiD color space decode & endecode parts
Update5: 2009-12-20 0.1.5 Add \'seq4newbler\' and \'seq4phrap\' functions
Update6: 2010-03-15 0.1.6 Change \'fq2fa\' into \'fq2faq\'
Update7: 2010-05-05 0.1.7 Add std format for \'seq4phrap\'
Update8: 2010-06-19 0.1.8 Add \'shuffle\' and \'distribute\' function
Update9: 2010-06-19 0.1.9 Add \'sam2fq\' and \'sam2pe\'
Update10: 2010-07-14 0.2.0 Revised some bugs of \'sam2pe\'
Update11: 2010-08-05 0.2.1 Add \'-pe_strand\' option for \'seq4newbler\' and \'seq4phrap\'
Update12: 2010-09-02 0.2.2 Add \'-mp\' option for \'seq4newbler\'
Update13: 2010-09-03 0.2.3 Add \'qseq2fq\' function
Update14: 2010-09-08 0.2.4 Modify \'fq2faq\' function
Update15: 2010-09-25 0.2.5 Add \'cds2aa\' function
Update16: 2010-09-26 0.2.6 Add \'gs2pe\' function
Update17: 2011-05-22 0.2.7 Add \'sort2pe\' function
Update18: 2011-06-22 0.2.8 Revised \'sort2pe\' function
Update19: 2012-05-04 0.2.9 Revised \'sort2pe\' function
Update20: 2012 0.3.0 Add \'sam2vcf\', \'bam2vcf\', \'cns2vcf\', \'vcf2gff\', \'vcf2bed\' (on going...)
Update21: 2013-03-11 0.2.10 Add \'merge\'
);

my ($Seq,$Qual,$Fastq,$Prefix,$Format,$Header,$Existh,$Chr,$PElib,$PE_strand,$RC,$MatePair,$Linker,$Noq,$Author,$Version,$Help);
my ($Typeset,$Transet,$Outset,$Qualset,$Lenset,$Seed,$Diffbase,$QFormat);
my %opts;
GetOptions
(
	\%opts,
	"S:s"=>\$Seq,
	"Q:s"=>\$Qual,
	"FQ:s"=>\$Fastq,
	"p:s"=>\$Prefix,
	"f:s"=>\$Format,
	"qf:s"=>\$QFormat,
	"header:s"=>\$Header,
	"RC"=>\$RC,
	"eh"=>\$Existh,
	"chr"=>\$Chr,
	"pelib"=>\$PElib,
	"pe_strand:s"=>\$PE_strand,
	"mp:s"=>\$MatePair,
	"linker:s"=>\$Linker,
	"noq"=>\$Noq,
	"author"=>\$Author,
	"version"=>\$Version,
	"t:s"=>\$Typeset,
	"translate:s"=>\$Transet,
	"o:i"=>\$Outset,
	"qual:i"=>\$Qualset,
	"len:i"=>\$Lenset,
	"seed:i"=>\$Seed,
	"diff:i"=>\$Diffbase,
	"help"=>\$Help
);
Transformers("BENM") if (defined $Author);
die ($update_log) if (defined $Version);
die($usage) if ((@ARGV==0)||($Help));

$PE_strand ||= "FR";

$QFormat ||= "L";
my %qualth=('S'=>33,'X'=>59,'I'=>64,'J'=>66,'L'=>33);

if (!defined $Prefix)
{
	#warn("-- The default prefix name of output file is set to 'out'. Use '-p' at the command line to change the default.\n");
	$Prefix ||= "out";
}

# Solexa->Sanger quality conversion table
my @conv_table;
for (-64..64) {
	$conv_table[$_+64] = chr(int($qualth{$QFormat} + 10*log(1+10**($_/10.0))/log(10)+.499));
}
# Solexa->454 quality conversion table
my @conv_newbler;
for (-64..64) {
	$conv_newbler[$_+64] = int(10*log(1+10**($_/10.0))/log(10)+.499);
}

# SOLiD color code
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

#454 Paired End Titanium Linkers
my $titanium_f;
my $titanium_r;
if (defined $Linker)
{
	if (-f $Linker)
	{
		open (LK,$Linker) || die $!;
		my ($f,$r)=(0,0);
		while(<LK>)
		{
			if (/^\>/)
			{
				if ((/f$/i)||(/5/)||(/five/))
				{
					$f=1;
				}
				elsif ((/r$/i)||(/3/)||(/three/))
				{
					$r=1;
				}
			}
			else
			{
				chomp;
				$titanium_f.=$_ if ($f==1);
				$titanium_r.=$_ if ($r==1);
			}
		}
		close LK;
	}
	else
	{
		$titanium_f=$Linker;
		$titanium_r=reverse(uc($Linker));
		$titanium_r=~tr/ACTG/TGCA/;
	}
}
if ((!defined $titanium_f)||(!defined $titanium_r))
{
	if ((defined $titanium_f)&&(!defined $titanium_r))
	{
		$titanium_r=reverse(uc($titanium_f));
		$titanium_r=~tr/ACTG/TGCA/;
	}
	elsif ((!defined $titanium_f)&&(defined $titanium_r))
	{
		$titanium_f=reverse(uc($titanium_r));
		$titanium_f=~tr/ACTG/TGCA/;
	}
	else
	{
		$titanium_r ||= 'TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG';
		$titanium_f ||= 'CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA';
	}
}

## amino acid translate table
my %CODE = (
		"nuclear" =>
		{
			'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
			'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
			'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Acid
			'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Acid
			'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanine
			'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
			'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
			'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',                                             # Isoleucine
			'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
			'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
			'ATG' => 'M',                                                                         # Methionine
			'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
			'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
			'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
			'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',   # Arginine
			'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',   # Serine
			'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
			'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
			'TGG' => 'W',                                                                         # Tryptophan
			'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
			'TAA' => 'U', 'TAG' => 'U', 'TGA' => 'U'                                              # Stop
		}
	);

# parsing command line
my $cmd = shift;
my %cmd_hash = (scarf2std=>\&scarf2std, std2scarf=>\&std2scarf, fqint2std=>\&fqint2std, std2fqint=>\&std2fqint,qseq2fq=>\&qseq2fq,
	sol2std=>\&sol2std, std2sol=>\&std2sol, std2solfq=>\&std2solfq, fa2std=>\&fa2std, fq2faq=>\&fq2faq, fa2cs=>\&fa2cs, cs2fa=>\&cs2fa,
	seqal2fq=>\&seqal2fq, solid2fq=>\&solid2fq, fq2solid=>\&fq2solid, sam2fq=>\&sam2fq, sam2pe=>\&sam2pe, sam2vcf=>\&sam2vcf,
	cns2vcf=>\&cns2vcf, vcf2gff=>\&vcf2gff, vcf2bed=>\&vcf2bed,
	cds2aa=>\&cds2aa, gs2pe=>\&gs2pe, example=>\&example, instruction=>\&instruction,
	seq4newbler=>\&seq4newbler,seq4phrap=>\&seq4phrap,shuffle=>\&shuffle,distribute=>\&distribute,sort2pe=>\&sort2pe,mergepe=>\&mergepe,
	bed2gff=>\&bed2gff,
	);
if (defined($cmd_hash{$cmd})) {
  &{$cmd_hash{$cmd}};
} else {
  Transformers(0);
  die("** Unrecognized command $cmd\n");
}

#######################################################################################
#---------------------------------- Main  Function ----------------------------------#
######################################################################################
sub seqal2fq
{
	Transformers(1);
	die "No inputfile!\nExample: Transformers seqal2fq -S seqence.txt -Q qual.txt -p out\n" if ((!defined $Seq)||(!defined $Qual));
	open (IN1,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
	open (IN2,$Qual) || die "Fail to open Quality file:$Qual for reading\n";
	open (OUT,">$Prefix.fastq") || die "Fail to create FASTQ file: $Prefix.fastq for writing\n";
	while (<IN1>)
	{
		my ($name,$seq,$qual)=("","","");
		next if ($_=~/^#/);
		$name = $1 if ($_=~/^\>(.*)/);
		$seq = <IN1>;
		chomp $seq;
		$_ = <IN2>; if ($_=~/^\#/) {<IN2>;} $qual=<IN2>;
		chomp $qual;
		if (defined $RC)
		{
			RC(\$seq);
			RC(\$qual);
		}
		print OUT "\@$name\n$seq\n+$name\n$qual\n";
	}
	close IN1;
	close IN2;
	close OUT;
}

sub qseq2fq
{
	Transformers(1);
	die "No inputfile!\nExample: Transformers qseq2fq -S file.qseq -p out\n" if ((!defined $Seq)||(!defined $Prefix));
	open (OUT,">$Prefix.fastq");
	open (IN,$Seq) || die $!;
	while(<IN>)
	{
		my @t=split;
		my $name="@".(join ":",@t[0,2..5])."#$t[6]/$t[7]";
		my $seq=$t[8];
		my $qual=$t[9];
		$seq=~tr/ACGTN/N/c;
		print OUT "$name\n$seq\n+\n$qual\n";
	}
	close IN;
	close OUT;
}

sub fa2std
{
	Transformers(2);
	my $q = ((defined $Qualset)&&($Qualset>0)) ? $Qualset : 25;
	die "No inputfile!\nExample: Transformers fa2std -S file.fasta -p out\n" if ((!defined $Seq)||(!defined $Prefix));
	#warn("-- The default quality is set to $Qualset. Use '-qual' at the command line to change the default.\n");
	my $output=$Prefix.".fq";
	open (IN,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
	open (OUT,">$output") || die "Fail to create FASTQ file:$output for writing\n";
	while (<IN>)
	{
		if (/^>(.*)/)
		{
			print OUT "\@$1\n";
			$_ = <IN>;
			s/\s+$//;
			my $seq = $_;
			RC(\$seq) if (defined $RC);
			print OUT "$seq\n+\n", $q x length($seq), "\n";
		}
	}
	close IN;
	close OUT;
}

sub fq2faq
{
	Transformers(3);
	die "No inputfile!\nExample: Transformers fq2faq -FQ file.fastq -p out [-RC] [-noq]\n" if (!defined $Fastq);
	my $outfa = $Prefix.".fa";
	my $outqual = $Prefix.".qual";
	open (IN,$Fastq) || die "Fail to open FASTQ file:$Fastq for reading\n";
	open (FA,">$outfa") || die "Fail to create FASTA file:$outfa for writing\n";
	open (QUAL,">$outqual") || die "Fail to create FASTA file:$outqual for writing\n" unless (defined $Noq);
	while (<IN>)
	{
		if (/^@(.*)/)
		{
			print FA ">$1\n";
			print QUAL ">$1\n" unless (defined $Noq);
			$_ = <IN>;
			s/\s+$//;
			my $seq=$_;
			RC(\$seq) if (defined $RC);
			print FA "$seq\n";
			<IN>;
			$_ = <IN>;
			next if (defined $Noq);
			chomp;
			my $qual=$_;
			if ($qual=~/\d+\s+\d+/)
			{
				$qual=join ' ',(reverse(split /\s+/,$_));
			}
			else
			{
				$qual=reverse($qual);
			}
			print QUAL $_ unless (defined $Noq);
		}
	}
	close IN;
	close FA;
	close QUAL unless (defined $Noq);
}

sub fa2cs
{
	Transformers(4);
	die "No inputfile!\nExample: Transformers fa2cs -S seqence.fasta -p out\n" if (!defined $Seq);
	my $output=$Prefix.".csfasta";
	open (IN,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
	open (OUT,">$output") || die "Fail to create CSFASTA file:$output for writing\n";
	print OUT "# Title: $Prefix\_";
	while (<IN>)
	{
		if (/^>/)
		{
			print OUT $_;
			$_ = <IN>;
			s/\s+$//;
			$_ =~ s/N//ig;
			$_ =~ s/X//ig;
			my $seq = $_;
			RC(\$seq) if (defined $RC);
			$seq=base2col($seq) if ($seq ne "");
			print OUT "$seq\n";
		}
	}
	close IN;
	close OUT;
}

sub cs2fa
{
	Transformers(5);
	die "No inputfile!\nExample: Transformers cs2fa -S seqence.csfasta -p out\n" if (!defined $Seq);
	my $output=$Prefix.".fasta";
	open (IN,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
	open (OUT,">$output") || die "Fail to create FASTA file:$output for writing\n";
	my $title="";
	while (<IN>)
	{
		if (/Title\:\s+(\S+)/)
		{
			$title=$1;
		}
		elsif (/^>(\S+)/)
		{
			print OUT ">$title$1\n";
			$_ = <IN>;
			s/\s+$//;
			my $seq=col2base($_);
			RC(\$seq) if (defined $RC);
			print OUT "$seq\n";
		}
	}
	close IN;
	close OUT;
}

sub scarf2std
{
	Transformers(6);
	die "No inputfile!\nExample: Transformers scarf2std -FQ seqence.scarf.txt -p out\n" if (!defined $Fastq);
	my $output=$Prefix.".fq";
	open (IN,$Fastq) || die "Fail to open SCARF file:$Fastq for reading\n";
	open (OUT,">$output") || die "Fail to create FASTQ file:$output for writing\n";
	while (<IN>)
	{
		my @t = split(':', $_);
		my ($name,$seq,$qual)=("","","");
		my @q;
		if (($t[1]=~/^[ACGTN]\w+[ACGTN]$/i)&&(($t[2]=~/\d+\s+\d+/)||($_=~/^\d+$/)))
		{
			$name=$t[0];
			$seq=$t[1];
			@q = split (/\s/, $t[2]);
		}
		else
		{
			$name = join('_', @t[0..4]);
			$seq = $t[5];
			@q = split(/\s/, $t[6]);
		}
		$qual .= $conv_table[$_+64] for (@q);
		if (defined $RC)
		{
			RC(\$seq);
			RC(\$qual);
		}
		print OUT "\@$name\n$seq\n+\n$qual\n";
	}
	close IN;
	close OUT;
}

sub std2scarf
{
	Transformers(7);
	die "No inputfile!\nExample: Transformers std2scarf -FQ seqence.fq -p out\n" if (!defined $Fastq);
	my $output = $Prefix.".scarf.txt";
	open (IN,$Fastq) || die "Fail to open FASTA file:$Fastq for reading\n";
	open (OUT,">$output") || die "Fail to create SCARF file:$output for writing\n";
	while(<IN>)
	{
		if ($_=~/^\@(.*)/)
		{
			my ($name,$seq,$qual)=("","","");
			$name=$1;
			$seq = <IN>;
			chomp $seq;
			<IN>;
			$_ = <IN>;
			s/\s+$//;
			my @t = split "",$_;
			$qual .= sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (@t);
			if (defined $RC)
			{
				RC(\$seq);
				RC(\$qual);
			}
			print OUT "$name:$seq:$qual\n";
		}
	}
	close OUT;
}

sub fqint2std
{
	Transformers(8);
	die "No inputfile!\nExample: Transformers fqint2std -FQ seqence.fastq -p out\nOr: Transformers fqint2std -S seqence.txt -Q qual.txt -p out\n" if ((!defined $Fastq)&&((!defined $Seq)||(!defined $Qual)));
	my $output=$Prefix.".fq";
	open (OUT,">$output") || die "Fail to create FASTA file:$output for writing\n";
	if (defined $Fastq)
	{
		open (IN,$Fastq) || die "Fail to open FASTQ file:$Fastq for reading\n";
		while (<IN>)
		{
			if (/^@/)
			{
				print OUT $_;
				$_ = <IN>;
				s/\s+$//;
				my $seq = ($_=~/^[ACGT]\d+/) ? col2base($_) : $_;
				RC(\$seq) if (defined $RC);
				print OUT "$seq\n";
				$_ = <IN>; $_ = <IN>;
				my @t = split;
				my $qual = '';
				$qual .= $conv_table[$_+64] for (@t);
				#my $qual=$_;
				#$qual =~ s/\s*(\d+)\s*/chr($1+$qualth{$QFormat})/eg;
				RC(\$qual) if (defined $RC);
				print OUT "+\n$qual\n";
			}
		}
		close IN;
	}
	elsif ((defined $Seq)&&(defined $Qual))
	{
		my ($Title,$Name,$seq,$qual)=("","","","");
		open (IN1,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
		open (IN2,$Qual) || die "Fail to open Quality file:$Qual for reading\n";
		while(<IN2>)
		{
			if ($_=~/^\#\s+Title\:\s+(\S+)/)
			{
				$Title=$1;
				my $in1=<IN1>;
				while($in1 !~ /^\#\s+Title\:\s+(\S+)/)
				{
					$in1=<IN1>;
				}
			}
			else
			{
				if ($_=~/^#\s\w+/)
				{
					<IN2>;<IN1>;
				}
				elsif (/^[\@\>](\S+)/)
				{
					chomp $seq;
					chomp $qual;
					substr($qual,0,1,"") if (length($qual)>length($seq));
					warn "The length of quality is not match with seqence!\nSeq: $seq\nQual: $qual\n" if (length($qual)!= length($seq));
					print OUT "$Name\n$seq\n+\n$qual\n" if ($seq ne "")&&($qual ne "");
					($seq,$qual)=("","");
					$Name="@".$Title.$1;
					$_ = <IN1>;
					if (/^[\@\>]/)
					{
						$_ = <IN1>;
					}
					my $tmp_seq;
					while (($_ ne "")&&($_ !~ /^[\@\>]/)&&(!eof))
					{
						s/\s+$//;
						$tmp_seq.=$_;
						$_=<IN1>;
						chomp $_;
					}
					$tmp_seq .= $_ if ((eof)&&($_ !~ /^[\@\>]/));
					chomp $tmp_seq;
					$seq = (($tmp_seq=~/\d+/)) ? col2base($tmp_seq) : $tmp_seq;
					$_ = <IN2>;
					s/^\s+//;
					s/\s+$//;
					if (($_=~/\d+\s+\d+/)||($_=~/^\d+$/))
					{
						$qual .= " " if ($qual=~/\d+$/);
						my @t = split /\s+/,$_;
						$qual .= $conv_table[$_+64] for (@t);
					}
					else
					{
						$qual .= $_;
					}
				}
				else
				{
					s/^\s+//;
					s/\s+$//;
					if (($_=~/\d+\s+\d+/)||($_=~/^\d+$/))
					{
						$qual .= " " if ($qual=~/\d+$/);
						my @t = split /\s+/,$_;
						$qual .= $conv_table[$_+64] for (@t);
					}
					else
					{
						$qual .= $_;
					}
				}
			}
		}
		
		substr($qual,0,1,"") if (length($qual)>length($seq));
		if (defined $RC)
		{
			RC(\$seq);
			RC(\$qual);
		}
		print OUT "$Name\n$seq\n+\n$qual\n" if ($seq ne "")&&($qual ne "");
		close IN1;
		close IN2;
	}
	close OUT;
}

sub std2fqint
{
	Transformers(9);
	die "No inputfile!\nExample: Transformers std2fqint -FQ seqence.fq -p out -f outputformat\n" if (!defined $Fastq);
	if (!defined $Format)
	{
		warn("-- The default output format is set by -f. Use '-f' at the command line to change the default.\n");
		$Format = "SOLiD";
	}
	my ($seq_file,$qual_file) = ($Format =~ /SOLiD/i) ? ("$Prefix.csfasta","$Prefix\_QV.qual") :
	($Format =~ /3730/) ? ("$Prefix.seq","$Prefix.seq.qual") : ("$Prefix.fna","$Prefix.qual");
	open (IN,$Fastq) || die "Fail to open FASTA file:$Fastq for reading\n";
	open (OUT1,">$seq_file") || die "Fail to create Seq file:$seq_file for writing\n";
	open (OUT2,">$qual_file") || die "Fail to create Quality file:$qual_file for writing\n";
	while (<IN>)
	{
		if ($_=~/^\@(.*)/)
		{
			my ($name,$seq,$qual)=("","","");
			$name=$1;
			$seq = <IN>;
			chomp $seq;
			RC(\$seq) if (defined $RC);
			$seq = base2col($seq) if ($Format=~/SOLiD/i);
			<IN>;
			$_ = <IN>;
			s/\s+$//;
			my @t = split "",$_;
			$qual .= sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (@t);
			RC(\$qual) if (defined $RC);
			my $tmp = ($Format=~/SOLiD/i) ? $conv_table[int(rand(40))+64] : "";
			print OUT1 ">$name\n$seq\n";
			print OUT2 ">$name\n$tmp$qual\n";
		}
	}
	close IN;
	close OUT1;
	close OUT2;
}

sub sol2std
{
	Transformers(10);
	die "No inputfile!\nExample: Transformers sol2std -FQ seqence.fastq -p out\nOr: Transformers sol2std -S seqence.txt -Q qual.txt -p out\n" if ((!defined $Fastq)&&((!defined $Seq)||(!defined $Qual)));
	my $output=$Prefix.".fq";
	open (OUT,">$output") || die "Fail to create FASTA file:$output for writing\n";
	if (defined $Fastq)
	{
		open (IN,$Fastq) || die "Fail to open FASTQ file:$Fastq for reading\n";
		while (<IN>)
		{
			if (/^@/) {
				print OUT $_;
				$_ = <IN>;
				s/\s+$//;
				my $seq = $_;
				RC(\$seq) if (defined $RC);
				print OUT "$seq\n";
				$_ = <IN>; $_ = <IN>;
				s/\s+$//;
				my @t = split('', $_);
				my $qual = '';
				$qual .= $conv_table[ord($_)] for (@t);
				RC(\$qual) if (defined $RC);
				print OUT "+\n$qual\n";
			}
		}
		close IN;
	}
	elsif ((defined $Seq)&&(defined $Qual))
	{
		my ($Name,$seq,$qual)=("","","","");
		open (IN1,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
		open (IN2,$Qual) || die "Fail to open Quality file:$Qual for reading\n";
		while(<IN1>)
		{
			if ($_=~/[\@\>](.*)/)
			{
				$Name="@".$1;
				print OUT "$Name\n";
				$_ = <IN1>;
				s/\s+$//;
				$seq = $_;
				RC(\$seq) if (defined $RC);
				print OUT "$seq\n";
				$_ = <IN2>; $_ = <IN2>;
				s/\s+$//;
				my @t = split('', $_);
				$qual .= $conv_table[ord($_)] for (@t);
				RC(\$qual) if (defined $RC);
				print OUT "+\n$qual\n";
				$qual = "";
			}
		}
		close IN1;
		close IN2;
	}
	close OUT;
}

sub std2sol
{
	Transformers(11);
	die "No inputfile!\nExample: Transformers std2sol -FQ seqence.fq -p Solexa\n" if (!defined $Fastq);
	my ($seq_file,$qual_file)=("$Prefix\_sequence.txt","$Prefix\_qcal.txt");
	open (IN,$Fastq) || die "Fail to open FASTA file:$Fastq for reading\n";
	open (OUT1,">$seq_file") || die "Fail to create Seq file:$seq_file for writing\n";
	open (OUT2,">$qual_file") || die "Fail to create Quality file:$qual_file for writing\n";
	while(<IN>)
	{
		if ($_=~/^\@(.*)/)
		{
			my ($name,$seq,$qual)=("","","");
			$name=$1;
			$seq = <IN>;
			chomp $seq;
			RC(\$seq) if (defined $RC);
			<IN>;
			$_ = <IN>;
			s/\s+$//;
			my @t = split "",$_;
			$qual .= chr(ord($_)-$qualth{$QFormat}+64) for (@t);
			RC(\$qual) if (defined $RC);
			print OUT1 ">$name\n$seq\n";
			print OUT2 ">$name\n$qual\n";
		}
	}
	close IN;
	close OUT1;
	close OUT2;
}

sub std2solfq
{
	Transformers(12);
	die "No inputfile!\nExample: Transformers std2solfq -FQ seqence.fq -p Illumina:SOlexa\n" if (!defined $Fastq);
	my $solexa_fastq=$Prefix.".solfq";
	open (IN,$Fastq) || die "Fail to open FASTQ file:$Fastq for reading\n";
	open (OUT,">$solexa_fastq") || die "Fail to create Seq file:$solexa_fastq for writing\n";
	while(<IN>)
	{
		if ($_=~/^\@(.*)/)
		{
			my ($name,$seq,$qual)=("","","");
			$name=$1;
			$seq = <IN>;
			chomp $seq;
			RC(\$seq) if (defined $RC);
			<IN>;
			$_ = <IN>;
			s/\s+$//;
			my @t = split "",$_;
			$qual .= chr(ord($_)-$qualth{$QFormat}+64) for (@t);
			RC(\$qual) if (defined $RC);
			print OUT "\@$name\n$seq\n\+$name\n$qual\n";
		}
	}
	close IN;
	close OUT;
}

sub col2base
{
	my $col = shift;
	my @colors = split '',$col;
	my $string = $colors[0];
	if($string !~ /[acgtn]/i){
		warn "$col has no header base\n";
		return 0;
	}
	my $last_base = $string;
	my $current_base = '';
	for(my $i=1;$i<@colors;$i++)
	{
		if (($last_base=~/N/i)&&(($colors[$i] eq ".")||($colors[$i]==5)))
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
	substr($string,0,1,"");
	return $string;
}

sub base2col
{
	my $base = shift;
	my @bases = split '',uc($base);
	my $j=0;
	my $col_code = (defined $Header) ? $Header : "G";
	if (defined $Existh)
	{
		$col_code=$bases[0];
		$j = 1;
	}
	my $current_base = $col_code;
	my $code = '';
	for(my $i=$j;$i<@bases;$i++)
	{
		$code = (exists $colcode{"$current_base$bases[$i]"}) ? $colcode{"$current_base$bases[$i]"} : $othercode[int(rand(@othercode))];
		$col_code .= $code;
		$current_base = $bases[$i];
	}
	return $col_code;
}

sub solid2fq
{
	Transformers(13);
	die "No inputfile!\nExample: Transformers solid2fq -FQ solid.fq -p ABi:SOLiD\n" if (!defined $Fastq);
	my $solid_fastq=$Prefix.".fq";
#NCBIs SOLiD color space fastq data
#\@SRR015253.1 LIZ_20071025_2_GrimmondsMES_SS7747_13_23_38_F3 length=35
#T1011122220100230032132.2111111002.1
#+SRR015253.1 LIZ_20071025_2_GrimmondsMES_SS7747_13_23_38_F3 length=35
#!)+%.*%*+2'0%%%-%+%*5'%!%9+'%+<+0%!%
#@SRR015253.2 LIZ_20071025_2_GrimmondsMES_SS7747_13_23_119_F3 length=35
#T0101233211103200232333.2111211002.1
#+SRR015253.2 LIZ_20071025_2_GrimmondsMES_SS7747_13_23_119_F3 length=35
#!,.+'+')'390%%%%%%%'%%%!-<++++<99%!%
#@SRR015253.3 LIZ_20071025_2_GrimmondsMES_SS7747_13_23_146_F3 length=35
#T0312202213101213131111.1110131102.1
#+SRR015253.3 LIZ_20071025_2_GrimmondsMES_SS7747_13_23_146_F3 length=35
#!93<*/18+%:9%+075*%:;+6!3<26%/<%-%!%
	open (OUT,">$solid_fastq") || die "Fail to create Seq file:$solid_fastq for writing\n";
	open (IN,$Fastq) || die "Fail to open FASTQ file:$Fastq for reading\n"; 
	while (<IN>)
	{
		if (/^\@(.*)/) #(/(\d+)_(\d+)_(\d+)_([FR]3)/)
		{
			my $name=$1;
			if (/(\d+)_(\d+)_(\d+)_([FR]3)/)
			{
				my ($a,$b,$c)=($1,$2,$3);
				my $i = ($4 =~ /F3/) ? 1 : 2;
				$name = "$Prefix\:$a\:$b\:$c/$i";
			}
			$_ = <IN>;
			#$_ = substr(<IN>, 1);
			s/\s+$//;
			#tr/0123./ACGTN/;
			my $seq = col2base($_);
			RC(\$seq) if (defined $RC);
			<IN>; $_ = <IN>;
			s/\s+$//;
			if ($_=~/^\d+\s*/)
			{
				s/^(\d+)\s*//;
				s/(\d+)\s*/chr($1+$qualth{$QFormat})/eg;
			}
			else
			{
				substr($_,0,1,'') if (length($_)>length($seq));
			}
			my $qual=$_;
			RC(\$qual) if (defined $RC);
			print OUT qq/\@$name\n$seq\n+\n$qual\n/;
		}
	}
	close IN;
	close OUT;
}

sub fq2solid
{
	Transformers(14);
	die "No inputfile!\nExample: Transformers fq2solid -FQ FASTQ.fq -p ABi:SOLiD [-header G | -eh]\n" if (!defined $Fastq);
	if (!defined $Existh)
	{
		warn("-- The default header base is set to 'G'. Use '-header' at the command line to change the default.\n") if (!defined $Header);
		$Header ||= "G";
	}
	my $fastq=$Prefix.".solidfq";
	open (OUT,">$fastq") || die "Fail to create Seq file:$fastq for writing\n";
	open (IN,$Fastq) || die "Fail to open FASTQ file:$Fastq for reading\n"; 
	while (<IN>)
	{
		my $name = $_;
		$name =~ s/\:/\_/g;
		$name =~ s/\/1/\_F3/;
		$name =~ s/\/2/\_R3/;
		$_ = <IN>;
		s/\s+$//;
		my $seq = $_;
		$seq = base2col($seq);
		#$seq =~ tr/ACGTNacgtn/0123.0123./;
		#$seq = ((defined $Header)&&(!defined $Existh)) ? $Header.$seq : $seq;
		<IN>; $_ = <IN>;
		my $qual = (length($_)<length($seq)) ? "H".$_ : $_;
		if (defined $RC)
		{
			RC(\$seq);
			RC(\$qual);
		}
		print OUT qq/$name$seq\n+\n$qual/;
	}
	close IN;
	close OUT;
}

## standart Solexa paired-end read name format: @HWI-EAS6_2_FC12994_PE:8:1:678:777
## SOLiD Name format: >3_16_150_F3
## 454 Name format: @SRR000921.7  E7OAHBI01A9IFH  length=173
## 3730 Name format : read eeq03a02.b1 is univ fwd template: eeq03a02 library: eeq03
##   or : mgsaea0_000103.z1.scf CHROMAT_FILE: mgsaea0_000103.z1.scf PHD_FILE: mgsaea0_000103.z1.scf.phd.1 CHEM: unknown DYE: unknown TIME: Wed Dec 29 09:56:26 2004
sub seq4newbler
{
	Transformers(15);
	if ((!defined $Fastq)&&(!defined $Seq)&&(!defined $MatePair))
	{
		die "No inputfile!\nExample: Transformers seq4newbler -FQ seqence.fq [-S seq.fa -Q seq.qual | -mp mate_pair.prefix] -p out.fix -f inputformat [-pelib] [-pe_strand FR]\n";
	}
	if (!defined $Format)
	{
		warn("-- The default input format is set by -f. Use '-f' at the command line to change the default.\n");
		$Format = "3730";
	}
	if (($Format=~/3730/)||($Format=~/sanger/i))
	{
		if (defined $Fastq)
		{
			sangerfq4newbler($Fastq,$Prefix);
		}
		elsif ((defined $Seq)&&(defined $Qual))
		{
			sangersq4newbler($Seq,$Qual,$Prefix);
		}
		elsif (defined $Seq)
		{
			sangerseq4newbler($Seq,$Prefix);
		}
	}
	else
	{
		solexa4newbler();
	}
}

sub solexa4newbler
{
	my %mark;
	$mark{1}=($PE_strand=~/^F/i)?"F":"R";
	$mark{2}=($PE_strand=~/R$/i)?"R":"F";
	my %LR=(1,"left",2,"right");
	if (defined $RC)
	{
		$mark{1}=($mark{1} eq "F")?"R":"F";
		$mark{2}=($mark{2} eq "R")?"F":"R";
	}
	if (defined $MatePair)
	{
		open (IN,$MatePair) || die "Can't open matepair(paired end) $MatePair for reading\n";
		open (OUT1,">$Prefix.fna") || die "Can't write to $Prefix.fna\n";
		open (OUT2,">$Prefix.qual") || die "Can't write to $Prefix.qual\n";
		#open (OUT3,">$Prefix.xml") || die "Can't write to $Prefix.xml\n";
		#print OUT3 "\<\?xml version=\"1.0\"\?\>\n\<trace_volume\>\n";
		while(<IN>)
		{
			chomp;
			my $name1=$_;
			$_=<IN>;
			chomp;
			my $seq1=$_;
			<IN>;
			$_=<IN>;
			chomp;
			my $qual1=$_;
			$_=<IN>;
			chomp;
			my $name2=$_;
			$_=<IN>;
			chomp;
			my $seq2=$_;
			<IN>;
			$_=<IN>;
			chomp;
			my $qual2=$_;
			my ($print_name,$print_seq,$print_qual)=print_mp($name1,$seq1,$qual1,$name2,$seq2,$qual2);
			print OUT1 ">$print_name\n$print_seq";
			print OUT2 ">$print_name\n$print_qual";
			#print OUT3 $print_xml;
		}
		close IN;
		close OUT1;
		close OUT2;
		#print OUT3 "\<\/trace_volume\>\n";
		#close OUT3;
	}
	elsif (defined $Fastq)
	{
		open (IN,$Fastq) || die "Can't open fastq $Fastq for reading\n";
		open (OUT1,">$Prefix.fna") || die "Can't write to $Prefix.fna\n";
		open (OUT2,">$Prefix.qual") || die "Can't write to $Prefix.qual\n";
		while(<IN>)
		{
			s/\s*$//;
			if (/^\@(.*)\/([12])/)
			{
				my $name=$1;
				my ($head,$tail);
				$tail=$2;
				if ($name=~/(\S+)\:(\d+)\:(\d+)\:(\d+)\:(\d+)/)
				{
					$head="$1$2$3$4$5";
				}
				$head =~ s/\:/\_/g;
				my $lib=(split /\_/,$head)[0];
				my $seq=<IN>;
				chomp $seq;
				my $len=length($seq);
				RC(\$seq) if (defined $RC);
				print OUT1 (">$head\_$LR{$tail} template=$head dir=$mark{$tail} $lib length=$len ");
				print OUT1 "library=pairlab" if (defined $PElib);
				print OUT1 "\n$seq\n";
				<IN>;
				$_=<IN>;
				s/\s+$//;
				my @t=split "",$_;
				my $qual="";
				if ($Format =~ /solexa/i)
				{
					$qual .= sprintf ("%02d ",$conv_newbler[ord($_)]) for (@t);
				}
				else
				{
					$qual .= sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (@t);
				}
				chomp $qual;
				RC(\$seq) if (defined $RC);
				print OUT2 (">$head\_$LR{$tail} template=$head dir=$mark{$tail} $lib length=$len ");
				print OUT2 "library=pairlab" if (defined $PElib);
				print OUT2 "\n$qual\n";
			}
		}
		close IN;
		close OUT1;
		close OUT2;
	}
	elsif ((defined $Seq)&&(defined $Qual))
	{
		my ($Name,$seq,$qual)=("","","");
		open (IN1,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
		open (IN2,$Qual) || die "Fail to open Quality file:$Qual for reading\n";
		open (OUT1,">$Prefix.fna") || die "Can't write to $Prefix.fna\n";
		open (OUT2,">$Prefix.qual") || die "Can't write to $Prefix.qual\n";
		while(<IN1>)
		{
			if ($_=~/\>(.*)/)
			{
				my $name=$1;
				my ($head,$tail);
				$tail=$2;
				if ($name=~/(\S+)\:(\d+)\:(\d+)\:(\d+)\:(\d+)/)
				{
					$head="$1$2$3$4$5";
				}
				$head =~ s/\:/\_/g;
				my $lib=(split /\_/,$head)[0];
				$_ = <IN1>;
				s/\s+$//;
				$seq = $_;
				my $len=length($seq);
				RC(\$seq) if (defined $RC);
				$_ = <IN2>; $_ = <IN2>;
				s/\s+$//;
				my @t = split('', $_);
				if ($Format =~ /solexa/i)
				{
					$qual .= sprintf ("%02d ",$conv_newbler[ord($_)]) for (@t);
				}
				else
				{
					$qual .= sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (@t);
				}
				$Name = "$head\_$LR{$tail} template=$head dir=$mark{$tail} $lib length=$len ";
				$Name .= "library=pairlab" if (defined $PElib);
				RC(\$qual) if (defined $RC);
				print OUT1 "$Name\n$seq\n";
				print OUT2 "$Name\n$qual\n";
				$qual = "";
			}
		}
		close IN1;
		close IN2;
		close OUT1;
		close OUT2;
	}
}

sub print_mp
{
	my ($name1,$seq1,$qual1,$name2,$seq2,$qual2)=@_;
	my ($name_o,$seq_o,$qual_o)=("","","");
	if ($name1=~/^\@(\S+)\/([12])$/)
	{
		my ($a,$b)=($1,$2);
		my ($c,$d)=($1,$2) if ($name2=~/^\@(\S+)\/([12])$/);
		if ($a eq $c)
		{
			if ($b==2)
			{
				my $tmp_seq=$seq2;
				my $tmp_qual=$qual2;
				$seq2=$seq1;
				$qual2=$qual1;
				$seq1=$tmp_seq;
				$qual1=$tmp_qual;
			}
			if ($a=~/(\S+)\:(\d+)\:(\d+)\:(\d+)\:(\d+)/)
			{
				$name_o="$1$2$3$4$5 rank=".sprintf("%03d",$2).sprintf("%03d",$3)." x=$4.0 y=$5.0 length=";
			}
		}
		else
		{
			warn "$name1 not eq $name2\n";
		}
	}
	my $linker;
	if ($PE_strand eq "FR")
	{
		RC(\$seq2);
		RC(\$qual2);
		$linker=$titanium_f;
	}
	elsif ($PE_strand eq "RF")
	{
		RC(\$seq2);
		RC(\$qual2);
		my $tmp_seq=$seq2;
		my $tmp_qual=$qual2;
		$seq2=$seq1;
		$qual2=$qual1;
		$seq1=$tmp_seq;
		$qual1=$tmp_qual;
		$linker=$titanium_r;
	}
	elsif ($PE_strand eq "FF")
	{
		$linker=$titanium_f;
	}
	elsif ($PE_strand eq "RR")
	{
		my $tmp_seq=$seq2;
		my $tmp_qual=$qual2;
		$seq2=$seq1;
		$qual2=$qual1;
		$seq1=$tmp_seq;
		$qual1=$tmp_qual;
		$linker=$titanium_r;
	}
	$seq_o=$seq1.$linker.$seq2;
	if (($qual1=~/\d{2}\s\d{2}/)&&($qual2=~/\d{2}\s\d{2}/))
	{
		$qual_o=$qual1." ".("40 "x(length($linker)))." ".$qual2;
	}
	else
	{
		if ($Format =~ /solexa/i)
		{
			$qual_o.=sprintf ("%02d ",$conv_newbler[ord($_)]) for (split '', $qual1);
			$qual_o.=" ".("40 "x(length($linker)))." ";
			$qual_o.=sprintf ("%02d ",$conv_newbler[ord($_)]) for (split '', $qual2);
		}
		else
		{
			$qual_o.=sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (split '', $qual1);
			$qual_o.=" ".("40 "x(length($linker)))." ";
			$qual_o.=sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (split '', $qual2);
		}
	}
	$name_o.=length($seq_o);
	Display_str(\$seq_o);
	Display_str(\$qual_o);
	return ($name_o,$seq_o,$qual_o);
}

sub Display_str
{
	my $str=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;
	$$str =~ s/\s$//;
	if ($$str=~/\d+\s\d+/)
	{
		my @t=split /\s+/,$$str;
		for (my $i=0;$i<@t;$i+=$num_line)
		{
			if ($i+59<@t-1)
			{
				$disp .= (join " ",@t[$i..($i+59)])."\n";
			}
			else
			{
				$disp .= (join " ",@t[$i..(@t-1)])."\n";
			}
		}
	}
	else
	{
		for (my $i=0; $i<length($$str); $i+=$num_line)
		{
			$disp .= substr($$str,$i,$num_line)."\n";
		}
	}
	$$str = ($disp) ?  $disp : "\n";
}

sub sangerfq4newbler
{
	my ($input,$outprefix)=@_;
	open (IN,$input) || die "Can't open fastq $input for reading\n";
	open (OUT1,">$outprefix.fna") || die "Can't write to $outprefix.fna\n";
	open (OUT2,">$outprefix.qual") || die "Can't write to $outprefix.qual\n";
	while(<IN>)
	{
		if (/^\@(.*)/)
		{
			my $name=$1;
			my $mark;
			my $head = $1 if ($name=~/^(\S+)/);
			$head=~s/\.[zy]1\.scf//g;
			$head=~s/[FR]//;
			$head=~s/[\-\_]$//;
			if (($name=~/y/)||($name=~/fwd/)||($name=~/F/))
			{
				$mark=(defined $RC)?"R":"F";
			}
			elsif (($name=~/z/)||($name=~/rev/)||($name=~/R/))
			{
				$mark=(defined $RC)?"F":"R";
			}
			elsif ($head=~/(\S+)\.ab1$/)
			{
				$head=$1;
				$mark="ab1";
			}
			else
			{
				warn "unrecognized name\n";
			}
			my $lib=(split /\_/,$head)[0];
			my $seq=<IN>;
			chomp $seq;
			my $len=length($seq);
			RC(\$seq) if (defined $RC);
			print OUT1 (">$head$mark template=$head dir=$mark $lib length=$len ");
			print OUT1 "library=pairlab" if (defined $PElib);
			print OUT1 "\n$seq\n";
			<IN>;
			$_=<IN>;
			s/\s+$//;
			my $qual="";
			if ($_=~/^\d+\s+\d+/)
			{
				$qual = $_;
			}
			else
			{
				my @t = split "",$_;
				$qual .= sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (@t);
			}
			RC(\$qual) if (defined $RC);
			print OUT2 (">$head$mark template=$head dir=$mark $lib length=$len ");
			print OUT2 "library=pairlab" if (defined $PElib);
			print OUT2 "\n$qual\n";
		}
	}
	close IN;
	close OUT1;
	close OUT2;
}

sub sangersq4newbler
{
	my ($seq_file,$qual_file,$outprefix)=@_;
	my ($Name,$head,$mark,$seq,$qual)=("","","","","");
	open (IN1,$seq_file) || die "Fail to open Seq file:$seq_file for reading\n";
	open (IN2,$qual_file) || die "Fail to open Quality file:$qual_file for reading\n";
	open (OUT1,">$outprefix.fna") || die "Can't write to $outprefix.fna\n";
	open (OUT2,">$outprefix.qual") || die "Can't write to $outprefix.qual\n";
	while(<IN2>)
	{
		if (/^[\@\>](.*)/)
		{
			my $name=$1;
			chomp $seq;
			chomp $qual;
			#substr($qual,0,1,"") if (length($qual)>length($seq));
			#warn "The length of quality is not match with seqence!\nSeq: $seq\nQual: $qual\n" if (length($qual)!= length($seq));
			if (($seq ne "")&&($qual ne ""))
			{
				my $len=length($seq);
				$Name.="length=$len ";
				$Name.="library=pairlab" if (defined $PElib);
				print OUT1 "$Name\n$seq\n";
				print OUT2 "$Name\n$qual\n";
				$Name="";
			}
			($seq,$qual)=("","");
			my $head = $1 if ($name=~/^(\S+)/);
			$head=~s/\.[zy]1\.scf//g;
			$head=~s/[FR]//;
			$head=~s/[\-\_]$//;
			if (($name=~/y/)||($name=~/fwd/)||($name=~/F/))
			{
				$mark=(defined $RC)?"R":"F";
			}
			elsif (($name=~/z/)||($name=~/rev/)||($name=~/R/))
			{
				$mark=(defined $RC)?"F":"R";
			}
			elsif ($head=~/(\S+)\.ab1$/)
			{
				$head=$1;
				$mark="ab1";
			}
			else
			{
				warn "unrecognized name\n";
			}
			my $lib=(split /\_/,$head)[0];
			$Name=">$head.$mark template=$head dir=$mark $lib ";
			$_ = <IN1>;
			if (/^[\@\>](.*)/)
			{
				die "sequence name not equal in seq file and qual file: $name ne $1" if ($name ne $1);
				$_ = <IN1>;
			}
			my $tmp_seq;
			while (($_ ne "")&&($_ !~ /^[\@\>]/)&&(!eof))
			{
				s/\s+$//g;
				$tmp_seq.=$_;
				$_=<IN1>;
				chomp $_;
			}
			$tmp_seq .= $_ if ((eof)&&($_ !~ /^[\@\>]/));
			chomp $tmp_seq;
			$seq = (($tmp_seq=~/\d+/)) ? col2base($tmp_seq) : $tmp_seq;
			$_ = <IN2>;
			s/\s+$//;
			s/^\s+//;
			if (($_=~/\d+\s+\d+/)||($_=~/^\d+$/))
			{
				$qual .= " " if ($qual=~/\d+$/);
				$qual .= $_;
			}
			else
			{
				my @t = split "",$_;
				$qual .= sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (@t);
			}
		}
		else
		{
			s/\s+$//;
			s/^\s+//;
			if (($_=~/\d+\s+\d+/)||($_=~/^\d+$/))
			{
				$qual .= " " if ($qual=~/\d+$/);
				$qual .= $_;
			}
			else
			{
				my @t = split "",$_;
				$qual .= sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (@t);
			}
		}
	}
	if (defined $RC)
	{
		RC(\$seq);
		RC(\$qual);
	}
	substr($qual,0,1,"") if (length($qual)>length($seq));
	my $seqlen=length($seq);
	$Name.="length=$seqlen ";
	$Name.="library=pairlab" if (defined $PElib);
	print OUT1 "$Name\n$seq\n";
	print OUT2 "$Name\n$qual\n";
	close IN1;
	close IN2;
}

sub sangerseq4newbler
{
	my ($seq_file,$outprefix)=@_;
	my ($Name,$head,$mark,$seq,$qual,$len)=("","","","","",0);
	open (IN,$seq_file) || die "Fail to open Seq file:$seq_file for reading\n";
	open (OUT,">$outprefix.fna") || die "Can't write to $outprefix.fna\n";
	while(<IN>)
	{
		if (/^\>(.*)/)
		{
			$len=length($seq);
			if ($len>0)
			{
				$Name .="length=$len ";
				$Name .="library=pairlab" if (defined $PElib);
				print OUT "$Name\n";
				RC(\$seq) if (defined $RC);
				print OUT "$seq\n";
			}
			$seq="";
			my $name=$1;
			my $head = $1 if ($name=~/^(\S+)/);
			$head=~s/\.[zy]1\.scf//g;
			$head=~s/[FR]//;
			$head=~s/[\-\_]$//;
			if (($name=~/y/)||($name=~/fwd/)||($name=~/F/))
			{
				$mark=(defined $RC)?"R":"F";
			}
			elsif (($name=~/z/)||($name=~/rev/)||($name=~/R/))
			{
				$mark=(defined $RC)?"F":"R";
			}
			elsif ($head=~/(\S+)\.ab1$/)
			{
				$head=$1;
				$mark="ab1"
			}
			else
			{
				warn "unrecognized name\n";
			}
			my $lib=(split /\_/,$head)[0];
			$Name = ">$head.$mark template=$head dir=$mark $lib ";
		}
		else
		{
			s/^\s+//;
			s/\s*$//;
			$seq.=$_;
		}
	}
	$len=length($seq);
	if ($len>0)
	{
		$Name .="length=$len ";
		$Name .="library=pairlab" if (defined $PElib);
		print OUT "$Name\n";
		RC(\$seq) if (defined $RC);
		print OUT "$seq\n";
		$seq="";
	}
	close OUT;
	close IN;
}

sub gs2pe
{
	Transformers(26);
	die "No inputfile!\nExample: Transformers gs2pe -S seq.fna -Q seq.qual -p output_prefix [-linker linker.fa] [-pe_strand output_strand,default:$PE_strand]\n" if ((!defined $Seq)||(!defined $Qual)||(!defined $Prefix));
	$Format ||= "454";
	my ($name,$head,$seq,$qual)=("","","","");
	open (IN1,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
	open (IN2,$Qual) || die "Fail to open Quality file:$Qual for reading\n";
	open (PE1,">$Prefix\_PE1.fastq") || die "Can't write to $Prefix\_PE1.fastq\n";
	open (PE2,">$Prefix\_PE2.fastq") || die "Can't write to $Prefix\_PE2.fastq\n";
	open (OT,">$Prefix\_unidentified.fastq") || die "Can't write to $Prefix\_unidentified.fastq\n";
	while(<IN2>)
	{
		if (/^\>(.*)/)
		{
			if (($seq ne "")&&($qual ne ""))
			{
				my ($pe1,$pe2)=check_454PE($name,$seq,$qual);
				if ((defined $pe1)&&(defined $pe2))
				{
					print PE1 $pe1;
					print PE2 $pe2;
				}
				else
				{
					print OT $pe1;
				}
			}
			$name=$1;
			($seq,$qual)=("","");
			$_ = <IN1>;
			if (/^\>(.*)/)
			{
				die "Error: sequence name not equal in seq file and qual file: $name ne $1" if ($name ne $1);
				$_ = <IN1>;
			}
			while (($_ ne "")&&($_ !~ /^\>/)&&(!eof))
			{
				s/\s+$//;
				s/^\s+//;
				$seq.=$_;
				$_=<IN1>;
				chomp $_;
			}
			$seq .= $_ if ((eof)&&($_ !~ /^[\@\>]/));
			$_ = <IN2>;
			s/\s+$//;
			s/^\s+//;
			if (($_=~/\d+\s+\d+/)||($_=~/^\d+$/))
			{
				$qual .= " " if ($qual=~/\d+$/);
				$qual .= $_;
			}
		}
		else
		{
			s/\s+$//;
			s/^\s+//;
			if (($_=~/\d+\s+\d+/)||($_=~/^\d+$/))
			{
				$qual .= " " if ($qual=~/\d+$/);
				$qual .= $_;
			}
		}
	}
	if (($seq ne "")&&($qual ne ""))
	{
		my ($pe1,$pe2)=check_454PE($name,$seq,$qual);
		if ((defined $pe1)&&(defined $pe2))
		{
			print PE1 $pe1;
			print PE2 $pe2;
		}
		else
		{
			print OT $pe1;
		}
	}
	close IN1;
	close IN2;
	close PE1;
	close PE2;
	close OT;
}

## 3730 Name format : read eeq03a02.b1 is univ fwd template: eeq03a02 library: eeq03
sub seq4phrap
{
	Transformers(16);
	if ((!defined $Fastq)&&(!defined $Seq))
	{
		die "No inputfile!\nExample: Transformers seq4phrap -FQ seqence.fq [-S seq.fa -Q seq.qual] -p out.prefix -f inputformat  [-pe_strand FR]\n";
	}
	if (!defined $Format)
	{
		warn("-- The default input format is set by -f. Use '-f' at the command line to change the default.\n");
		$Format = "454";
	}
	if (($Format=~/454/)||($Format=~/roche/i))
	{
		if (defined $Fastq)
		{
			open (IN,$Fastq) || die "Can't open fastq $Fastq for reading\n";
			open (OUT1,">$Prefix.seq") || die "Can't write to $Prefix.seq\n";
			open (OUT2,">$Prefix.seq.qual") || die "Can't write to $Prefix.seq.qual\n";
			while(<IN>)
			{
				if (/^\@(.*)/)
				{
					my $name=$1;
					my ($mark1,$mark2);
					my $head = $1 if ($name=~/^(\S+)/);
					if (($head=~/F$/)||($name=~/die\=F/))
					{
						$mark1=(defined $RC)?"z1":"y1";
						$mark2=(defined $RC)?"rev":"fwd";
					}
					elsif (($head=~/R$/)||($name=~/die\=R/))
					{
						$mark1=(defined $RC)?"y1":"z1";
						$mark2=(defined $RC)?"fwd":"rev";
					}
					else
					{
						warn "unrecognized name\n";
					}
					my $lib=(split /\_/,$head)[0];
					my $seq=<IN>;
					chomp $seq;
					my $len=length($seq);
					RC(\$seq) if (defined $RC);
					print OUT1 (">$head.$mark1.scf is univ $mark2 template: $head library: $lib length: $len\n");
					print OUT1 "$seq\n";
					<IN>;
					$_=<IN>;
					s/\s+$//;
					my $qual="";
					if ($_=~/^\d+\s+\d+/)
					{
						$qual = $_;
					}
					else
					{
						my @t = split "",$_;
						$qual .= sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (@t);
					}
					RC(\$qual) if (defined $RC);
					print OUT2 (">$head.$mark1.scf is univ $mark2 template: $head library: $lib length: $len\n");
					print OUT2 "$qual\n";
				}
			}
			close IN;
			close OUT1;
			close OUT2;
		}
		elsif ((defined $Seq)&&(defined $Qual))
		{
			my ($Name,$head,$seq,$qual)=("","","","");
			open (IN1,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
			open (IN2,$Qual) || die "Fail to open Quality file:$Qual for reading\n";
			open (OUT1,">$Prefix.seq") || die "Can't write to $Prefix.seq\n";
			open (OUT2,">$Prefix.seq.qual") || die "Can't write to $Prefix.seq.qual\n";
			while(<IN2>)
			{
				if (/^\>(.*)/)
				{
					my $name=$1;
					chomp $seq;
					chomp $qual;
					#substr($qual,0,1,"") if (length($qual)>length($seq));
					#warn "The length of quality is not match with seqence!\nSeq: $seq\nQual: $qual\n" if (length($qual)!= length($seq));
					if (($seq ne "")&&($qual ne ""))
					{
						my $len=length($seq);
						$Name.="length: $len ";
						if (defined $RC)
						{
							RC(\$seq);
							RC(\$qual);
						}
						print OUT1 "$Name\n$seq\n";
						print OUT2 "$Name\n$qual\n";
						$Name="";
					}
					($seq,$qual)=("","");
					my $head = $1 if ($name=~/^(\S+)/);
					my ($mark1,$mark2);
					if (($head=~/F$/)||($name=~/die\=F/))
					{
						$mark1=(defined $RC)?"z1":"y1";
						$mark2=(defined $RC)?"rev":"fwd";
					}
					elsif (($head=~/R$/)||($name=~/die\=R/))
					{
						$mark1=(defined $RC)?"y1":"z1";
						$mark2=(defined $RC)?"fwd":"rev";
					}
					else
					{
						warn "unrecognized name\n";
					}
					my $lib=(split /\_/,$head)[0];
					$Name=">$head.$mark1.scf is univ $mark2 template: $head library: $lib";
					$_ = <IN1>;
					if (/^[\@\>](.*)/)
					{
						die "sequence name not equal in seq file and qual file: $name ne $1" if ($name ne $1);
						$_ = <IN1>;
					}
					my $tmp_seq;
					while (($_ ne "")&&($_ !~ /^[\@\>]/)&&(!eof))
					{
						s/\s+$//;
						s/^\s+//;
						$tmp_seq.=$_;
						$_=<IN1>;
						chomp $_;
					}
					$tmp_seq .= $_ if ((eof)&&($_ !~ /^[\@\>]/));
					chomp $tmp_seq;
					$seq = (($tmp_seq=~/\d+/)) ? col2base($tmp_seq) : $tmp_seq;
					$_ = <IN2>;
					s/\s+$//;
					s/^\s+//;
					if (($_=~/\d+\s+\d+/)||($_=~/^\d+$/))
					{
						$qual .= " " if ($qual=~/\d+$/);
						$qual .= $_;
					}
					else
					{
						my @t = split "",$_;
						$qual .= sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (@t);
					}
				}
				else
				{
					s/\s+$//;
					s/^\s+//;
					if (($_=~/\d+\s+\d+/)||($_=~/^\d+$/))
					{
						$qual .= " " if ($qual=~/\d+$/);
						$qual .= $_;
					}
					else
					{
						my @t = split "",$_;
						$qual .= sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (@t);
					}
				}
			}
			if (defined $RC)
			{
				RC(\$seq);
				RC(\$qual);
			}
			substr($qual,0,1,"") if (length($qual)>length($seq));
			my $seqlen=length($seq);
			$Name.="length: $seqlen";
			print OUT1 "$Name\n$seq\n";
			print OUT2 "$Name\n$qual\n";
			close IN1;
			close IN2;
		}
		elsif (defined $Seq)
		{
			my ($Name,$head,$seq,$qual,$len)=("","","","",0);
			open (IN,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
			open (OUT,">$Prefix.seq") || die "Can't write to $Prefix.seq\n";
			while(<IN>)
			{
				if (/^\>(.*)/)
				{
					$len=length($seq);
					if ($len>0)
					{
						$Name.="length: $len\n";
						print OUT $Name;
						RC(\$seq) if (defined $RC);
						print OUT "$seq\n";
						$seq="";
					}
					my $name=$1;
					my $head = $1 if ($name=~/^(\S+)/);
					my ($mark1,$mark2);
					if (($head=~/F$/)||($name=~/die\=F/))
					{
						$mark1=(defined $RC)?"z1":"y1";
						$mark2=(defined $RC)?"rev":"fwd";
					}
					elsif (($head=~/R$/)||($name=~/die\=R/))
					{
						$mark1=(defined $RC)?"y1":"z1";
						$mark2=(defined $RC)?"fwd":"rev";
					}
					else
					{
						warn "unrecognized name\n";
					}
					my $lib=(split /\_/,$head)[0];
					$Name = ">$head.$mark1.scf is univ $mark2 template: $head library: $lib ";
				}
				else
				{
					s/^\s+//;
					s/\s+$//;
					$seq.=$_;
				}
			}
			$len=length($seq);
			if ($len>0)
			{
				$Name.="length: $len\n";
				print OUT $Name;
				RC(\$seq) if (defined $RC);
				print OUT "$seq\n";
			}
			close OUT;
			close IN;
		}
	}
	elsif (($Format=~/solexa/i)||($Format=~/illumina/i)||($Format=~/std/i)||($Format=~/sanger/i))
	{
		my %mark;
		$mark{1}{1}=($PE_strand=~/^F/i)?"y1":"z1";
		$mark{1}{2}=($PE_strand=~/R$/i)?"z1":"y1";
		$mark{2}{1}=($PE_strand=~/^F/i)?"fwd":"rev";
		$mark{2}{2}=($PE_strand=~/R$/i)?"rev":"fwd";
		if (defined $RC)
		{
			$mark{1}{1}=($mark{1}{1} eq "y1")?"z1":"y1";
			$mark{1}{2}=($mark{1}{2} eq "z1")?"y1":"z1";
			$mark{2}{1}=($mark{2}{1} eq "fwd")?"rev":"fwd";
			$mark{2}{2}=($mark{2}{2} eq "rev")?"fwd":"rev";
		}
		if (defined $Fastq)
		{
			open (IN,$Fastq) || die "Can't open fastq $Fastq for reading\n";
			open (OUT1,">$Prefix.seq") || die "Can't write to $Prefix.seq\n";
			open (OUT2,">$Prefix.seq.qual") || die "Can't write to $Prefix.seq.qual\n";
			while(<IN>)
			{
				if (/^\@(.*)/)
				{
					my $name=$1;
					my ($head,$tail)=($1,$2) if ($name=~/(\S+)\/(\d{1})/);
					$head =~ s/\:/\_/g;
					my $lib=(split /\_/,$head)[0];
					my $seq=<IN>;
					chomp $seq;
					my $len=length($seq);
					RC(\$seq) if (defined $RC);
					print OUT1 (">$head.$mark{1}{$tail}.scf is univ $mark{2}{$tail} template: $head library: $lib length: $len\n");
					print OUT1 "$seq\n";
					<IN>;
					$_=<IN>;
					s/\s+$//;
					my @t=split "",$_;
					my $qual="";
					if ($Format=~/solexa/i)
					{
						$qual .= sprintf ("%02d ",$conv_newbler[ord($_)]) for (@t);
					}
					else
					{
						$qual .= sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (@t);
					}
					chomp $qual;
					RC(\$qual) if (defined $RC);
					print OUT2 (">$head.$mark{1}{$tail}.scf is univ $mark{2}{$tail} template: $head library: $lib length: $len\n");
					print OUT2 "$qual\n";
				}
			}
			close IN;
			close OUT1;
			close OUT2;
		}
		elsif ((defined $Seq)&&(defined $Qual))
		{
			my ($Name,$seq,$qual)=("","","");
			open (IN1,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
			open (IN2,$Qual) || die "Fail to open Quality file:$Qual for reading\n";
			open (OUT1,">$Prefix.seq") || die "Can't write to $Prefix.seq\n";
			open (OUT2,">$Prefix.seq.qual") || die "Can't write to $Prefix.seq.qual\n";
			while(<IN1>)
			{
				if ($_=~/[\@\>](.*)/)
				{
					my $name=$1;
					my ($head,$tail)=($1,$2) if ($name=~/(\S+)\/(\d{1})/);
					$head =~ s/\:/\_/g;
					my $lib=(split /\_/,$head)[0];
					$_ = <IN1>;
					s/\s+$//;
					$seq = $_;
					my $len=length($seq);
					RC(\$seq) if (defined $RC);
					$_ = <IN2>; $_ = <IN2>;
					s/\s+$//;
					my @t = split('', $_);
					if ($Format=~/solexa/i)
					{
						$qual .= sprintf ("%02d ",$conv_newbler[ord($_)]) for (@t);
					}
					else
					{
						$qual .= sprintf ("%02d ",(ord($_)-$qualth{$QFormat})) for (@t);
					}
					$Name = ">$head.$mark{1}{$tail}scf is univ $mark{2}{$tail} template: $head library: $lib length: $len";
					RC(\$qual) if (defined $RC);
					print OUT1 "$Name\n$seq\n";
					print OUT2 "$Name\n$qual\n";
					$qual = "";
				}
			}
			close IN1;
			close IN2;
			close OUT1;
			close OUT2;
		}
		elsif ((defined $Seq)&&(!defined $Qual))
		{
			my ($Name,$seq)=("","");
			my ($lib,$head,$tail,$len);
			open (IN,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
			open (OUT,">$Prefix.seq") || die "Can't write to $Prefix.seq\n";
			while(<IN>)
			{
				if ($_=~/[\@\>](.*)/)
				{
					my $name=$1;
					if ($seq ne "")
					{
						$len=length($seq);
						RC(\$seq) if (defined $RC);
						$Name = ">$head.$mark{1}{$tail}scf is univ $mark{2}{$tail} template: $head library: $lib length: $len";
						print OUT "$Name\n$seq\n";
					}
					($head,$tail)=($1,$2) if ($name=~/(\S+)\/(\d{1})/);
					$head =~ s/\:/\_/g;
					$lib=(split /\_/,$head)[0];
					$seq="";
				}
				else
				{
					s/\s*//g;
					$seq.=$_ if ($_=~/^[ACGTN]+$/i);
				}
			}
			if ($seq ne "")
			{
				$len=length($seq);
				RC(\$seq) if (defined $RC);
				$Name = ">$head.$mark{1}{$tail}scf is univ $mark{2}{$tail} template: $head library: $lib length: $len";
				print OUT "$Name\n$seq\n";
			}
			close IN;
			close OUT;
		}
	}
}

sub shuffle
{
	Transformers(17);
	die "Usage:\nTransformers shuffle <IN1:PE1.fa|PE1.fq> <IN2:PE2.fa|PE2.fq> <OUT:out.fa|out.fq>\n" if @ARGV!=3;
	my ($PE1,$PE2,$SEQ)=@ARGV[0..2];
	my $aq1=check_fafq($PE1);
	my $aq2=check_fafq($PE2);
	die "\nInput sequences are not in the same format!\n" if ($aq1!=$aq2);
	open (IN1,"<$PE1") || die $!;
	open (IN2,"<$PE2") || die $!;
	open (OUT,">$SEQ") || die $!;
	my $out="";
	while(<IN1>)
	{
		if ($aq1==1)
		{
			print OUT $_;
			$_=<IN1>;
			my $seq1=$_;
			chomp $seq1;
			RC(\$seq1) if (defined $RC);
			print OUT "$seq1\n";

			$_=<IN2>;
			print OUT $_;
			$_=<IN2>;
			my $seq2=$_;
			chomp $seq2;
			RC(\$seq2) if (defined $RC);
			print OUT "$seq2\n";
		}
		else
		{
			print OUT $_;
			$_=<IN1>;
			my $seq1=$_;
			chomp $seq1;
			RC(\$seq1) if (defined $RC);
			print OUT "$seq1\n";
			$_=<IN1>;
			print OUT $_;
			$_=<IN1>;
			print OUT $_;

			$_=<IN2>;
			print OUT $_;
			$_=<IN2>;
			my $seq2=$_;
			chomp $seq2;
			RC(\$seq2) if (defined $RC);
			print OUT "$seq2\n";
			$_=<IN2>;
			print OUT $_;
			$_=<IN2>;
			print OUT $_;
		}
	}
	close IN1;
	close IN2;
	close OUT;
}

sub distribute
{
	Transformers(18);
	die "Usage:\nTransformers distribute <IN:PE.fa|PE.fq> <OUT1:out1.fa|out1.fq> <OUT2:out2.fa|out2.fq>\n" if @ARGV!=3;
	my ($SEQ,$PE1,$PE2)=@ARGV[0..2];
	my $aq=check_fafq($SEQ);
	open (IN,"<$SEQ") || die $!;
	open (OUT1,">$PE1") || die $!;
	open (OUT2,">$PE2") || die $!;
	while(<IN>)
	{
		if ($aq==1)
		{
			print OUT1 $_;
			$_=<IN>;
			my $seq1=$_;
			chomp $seq1;
			RC(\$seq1) if (defined $RC);
			print OUT1 "$seq1\n";
			
			$_=<IN>;
			print OUT2 $_;
			$_=<IN>;
			my $seq2=$_;
			chomp $seq2;
			RC(\$seq2) if (defined $RC);
			print OUT2 "$seq2\n";
		}
		else
		{
			print OUT1 $_;
			$_=<IN>;
			my $seq1=$_;
			chomp $seq1;
			RC(\$seq1) if (defined $RC);
			print OUT1 "$seq1\n";
			$_=<IN>;
			print OUT1 $_;
			$_=<IN>;
			print OUT1 $_;
			
			$_=<IN>;
			print OUT2 $_;
			$_=<IN>;
			my $seq2=$_;
			chomp $seq2;
			RC(\$seq2) if (defined $RC);
			print OUT2 "$seq2\n";
			$_=<IN>;
			print OUT2 $_;
			$_=<IN>;
			print OUT2 $_;
		}
	}
	close IN;
	close OUT1;
	close OUT2;
}

sub mergepe
{
	Transformers(17);
	die "Usage:\nTransformers shuffle <IN1:PE1.fa|PE1.fq> <IN2:PE2.fa|PE2.fq> <-p:output_prefix> [-seed 11] [-PE_strand FR]\n" if (@ARGV!=2 || !defined $Prefix);
	$Seed ||= 11;
	my ($Reads1,$Reads2)=@ARGV;
	my $aq1=check_fafq($Reads1);
	my $aq2=check_fafq($Reads2);
	my $suffix=($aq1==1)?"fasta":"fastq";
	die "$Reads1 is not the same format as $Reads2\n" if ($aq1 != $aq2);
	my ($Name1,$Name2,$name1,$name2,$seq1,$seq2,$qual1,$qual2)=();
	open (R1,$Reads1) || die $!;
	open (R2,$Reads2) || die $!;
	open (PE1,">$Prefix\_R1.$suffix") || die $!;
	open (PE2,">$Prefix\_R2.$suffix") || die $!;
	open (MO,">$Prefix\_merged.$suffix") || die $!;
	while(<R1>)
	{
		if ($aq1==1)
		{
			$name1=$1;
			$Name1=$_;
			$_=<R1>;chomp;
			$seq1=$_;
			$name2=$1 if (/^\@([^\/\s]+)/);
			$Name2=$_;
			$_=<R2>;chomp;
			$seq2=$_;
			die "Reads name error:\n$name1\n$name2\n" if ($name1 ne $name2);
			if ($PE_strand =~ /F$/)
			{
				$seq2=reverse($seq2);
				$seq2=~tr/ACGTacgt/TGCAtgca/;
			}
			if ($PE_strand =~ /^R/)
			{
				$seq1=reverse($seq1);
				$seq1=~tr/ACGTacgt/TGCAtgca/;
			}
			my ($mo,$seq)=MergeSeq($seq1,$seq2);
			if (defined $mo && defined $seq)
			{
				chomp $Name1;
				print MO "$Name1 $mo\n$seq\n";
			}
			else
			{
				print PE1 "$Name1$seq1\n";
				print PE2 "$Name2$seq2\n";
			}
		}
		else
		{
			$name1=$1;
			$Name1=$_;
			$_=<R1>;chomp;
			$seq1=$_;
			$_=<R1>;$_=<R1>;chomp;
			$qual1=$_;
			$_=<R2>;
			$name2=$1 if (/^\@([^\/\s]+)/);
			$Name2=$_;
			$_=<R2>;chomp;
			$seq2=$_;
			$_=<R2>;$_=<R2>;chomp;
			$qual2=$_;
			die "Reads name error:\n$name1\n$name2\n" if ($name1 ne $name2);
			if ($PE_strand =~ /F$/)
			{
				$seq2=reverse($seq2);
				$seq2=~tr/ACGTacgt/TGCAtgca/;
				$qual2=reverse($qual2);
			}
			if ($PE_strand =~ /^R/)
			{
				$seq1=reverse($seq1);
				$seq1=~tr/ACGTacgt/TGCAtgca/;
				$qual1=reverse($qual1);
			}
			my ($mo,$seq,$qual)=MergeSeq($seq1,$seq2,$qual1,$qual2);
			if (defined $mo && defined $seq)
			{
				chomp $Name1;
				print MO "$Name1 $mo\n$seq\n\+\n$qual\n";
			}
			else
			{
				print PE1 "$Name1$seq1\n\+\n$qual2\n";
				print PE2 "$Name2$seq2\n\+\n$qual2\n";
			}
		}
	}
	close R1;
	close R2;
	close PE1;
	close PE2;
	close MO;
}

sub MergeSeq
{
	my ($seq1,$seq2,$qual1,$qual2)=@_;
	my ($mo,$seq,$qual);
	my $seq2_rc=reverse($seq2);
	$seq2_rc=~tr/ACGTacgt/TGCAtgca/;
	for (my $i=0;$i<=length($seq2_rc)-$Seed;$i++)
	{
		my $seedSeq=substr($seq2_rc,$i,$Seed);
		if ($seq1=~/.+?($seedSeq)/i)
		{
			my $pos=$-[1];
			if ($pos<length($seq1)-$Seed)
			{
				my $mo=($pos+1)."-".($pos+$Seed)." ".(length($seq2_rc)-$i-$Seed)."-".(length($seq2_rc)-$i).":rc";
				my $f5seq=substr($seq1,0,$pos+$Seed);
				my $f3seq=substr($seq1,$pos+$Seed,length($seq2_rc)-$i-$Seed);
				my $r3seq=substr($seq2_rc,$i+$Seed,length($seq2_rc)-$i-$Seed);
				$seq="$f5seq$r3seq";
				if (defined $qual1 && defined $qual2)
				{
					my $qual2_r=reverse($qual2);
					my $comapre_qual=MaxQual(substr($qual1,$pos+1,$Seed),substr($qual2_r,$i,$Seed));
					return 0 if (!defined $comapre_qual);
					$qual=(substr($qual1,0,$pos)).$comapre_qual.(substr($qual2_r,$i+$Seed,length($seq2_rc)-$i-$Seed));
				}
				if (length($seq)<=length($seq1))
				{
					if (uc($f3seq) ne uc($r3seq))
					{
						my $diff=0;
						my @f3=split '',uc($f3seq);
						my @r3=split '',uc($r3seq);
						for (my $j=0;$j<@f3;$j++)
						{
							$diff++ if ($f3[$j] ne $r3[$j] && $f3[$j] ne "N" && $r3[$j] ne "N");
							last if ($diff>$Diffbase);
						}
						return 0 if ($diff>$Diffbase);
						if (defined $qual1 && defined $qual2)
						{
							my $qual2_r=reverse($qual2);
							my $f3qual=substr($qual1,$pos+$Seed,length($seq2_rc)-$i-$Seed);
							my $r3qual=substr($qual2_r,$i+$Seed,length($seq2_rc)-$i-$Seed);
							my $compare_qual=MaxQual($f3qual,$r3qual);
							return 0 if (!defined $compare_qual);
							if ($compare_qual eq $f3qual)
							{
								$seq="$f5seq$f3seq";
								$qual=substr($qual1,0,$pos+length($seq2_rc)-$i);
							}
							else
							{
								$qual=(substr($qual1,0,$pos+$Seed)).$r3qual;
							}
						}
						else
						{
							$seq="$f5seq$f3seq";
						}
					}
				}
				else
				{
					my $f_ol=substr($seq1,$pos,length($seq1)-$pos);
					my $r_ol=substr($seq2_rc,$i,length($seq1)-$pos);
					if (uc($f_ol) ne uc($r_ol))
					{
						my $diff=0;
						my @fo=split '',uc($f_ol);
						my @ro=split '',uc($r_ol);
						for (my $j=0;$j<@fo;$j++)
						{
							$diff++ if ($fo[$j] ne $ro[$j] && $fo[$j] ne "N" && $ro[$j] ne "N");
							last if ($diff>$Diffbase);
						}
						return 0 if ($diff>$Diffbase);
						if (defined $qual1 && defined $qual2)
						{
							my $qual2_r=reverse($qual2);
							my $foqual=substr($qual1,$pos,length($seq1)-$pos);
							my $roqual=substr($qual2_r,$i,length($seq1)-$pos);
							my $compare_qual=MaxQual($foqual,$roqual);
							return 0 if (!defined $compare_qual);
							if ($compare_qual eq $foqual)
							{
								$seq=$seq1.(substr($seq2_rc,$i+length($seq1)-$pos,length($seq2_rc)-$i+length($seq1)+$pos));
								$qual=$qual1.(substr($qual2_r,$i+length($qual1)-$pos,length($qual2_r)-$i+length($qual1)+$pos));
							}
							else
							{
								$qual=(substr($qual1,0,$pos)).$qual2_r;
							}
						}
						else
						{
							$seq="$f5seq$f3seq";
						}
					}
				}
				return ($mo,$seq,$qual);
			}
		}
	}
}

sub MaxQual
{
	my ($q1,$q2)=@_;
	my $q="";
	my @qual1=split '',$q1;
	my @qual2=split '',$q2;
	my $qq1=0;
	my $qq2=0;
	my $lowq1=0;
	my $lowq2=0;
	for (my $i=0;$i<@qual1;$i++)
	{
		$qq1+=ord($qual1[$i])-$qualth{$Format};
		$lowq1++ if (ord($qual1[$i])-$qualth{$Format}<=10);
		$qq2+=ord($qual2[$i])-$qualth{$Format};
		$lowq2++ if (ord($qual2[$i])-$qualth{$Format}<=10);
	}
	if (($qq1/@qual1>=15 || $qq2/@qual2>=15)&&($lowq1<5 || $lowq2<5))
	{
		$q=($qq1>$qq2)?$q1:$q2;
		return $q;
	}
}

sub sam2fq
{
	Transformers(19);
	die("Usage: \nTransformers sam2fq <IN:aln.sam> <-p:output_prefix> [-chr delimited_by_reference]\n") if ((@ARGV==0)||(!defined $Prefix));
	my $cmd;
	if (defined $Chr)
	{
		$cmd=qq(awk 'BEGIN{c=""}{if ((\$10!="")){c=\$3;if (c=="*"){c="none"}print "@"\$1"\\n"\$10"\\n+\\n"\$11 >> "$Prefix\."c"\.fq"}}' $ARGV[0]);
	}
	else
	{
		$cmd=qq(awk '{print "@"\$1"\\n"\$10"\\n+\\n"\$11}' $ARGV[0] > $Prefix.fq);
	}
	system("$cmd") && die("** fail to run command '$cmd'");
}

sub sam2pe
{
	Transformers(20);
	die("Usage: \nTransformers sam2pe <IN:aln.sam> <-p:output_prefix> [-chr delimited_by_reference]\n") if ((@ARGV==0)||(!defined $Prefix));
	my $cmd;
	if (defined $Chr)
	{
		$cmd=qq(awk '
BEGIN{p="";a="";b="";c="";d=""}\\
{\\
if (and(\$2,0x0040)==64){a="@"\$1"\/1\\n"\$10"\\n+\\n"\$11;c=\$3}\\
if (and(\$2,0x0080)==128){b="@"\$1"\/2\\n"\$10"\\n+\\n"\$11;d=\$3}\\
if ((\$10!="")&&(p==\$1)&&(and(\$2,0x0001)==1)&&(a!="")&&(b!="")){\\
if (c=="*"){c="none"}if(d=="*"){d="none"}\\
if (c==d){print a >> "$Prefix."c"\_PE1.fq";print b >> "$Prefix."d"\_PE2.fq";a="";b=""}\\
else{print a >> "$Prefix."c"-"d"\_PE1.fq";print b >> "$Prefix."c"-"d"\_PE2.fq";a="";b=""}}\\
p=\$1}' $ARGV[0]);
	}
	else
	{
		$cmd=qq(awk 'BEGIN{p="";a="";b=""}{if (and(\$2,0x0040)==64){a="@"\$1"\/1\\n"\$10"\\n+\\n"\$11}if(and(\$2,0x0080)==128){b="@"\$1"\/2\\n"\$10"\\n+\\n"\$11}if ((p==\$1)&&(and(\$2,0x0001)==1)&&(a!="")&&(b!="")){print a >> "$Prefix\_PE1.fq";print b >> "$Prefix\_PE2.fq";a="";b=""}p=\$1}' $ARGV[0]);
	}
	system("$cmd") && die("** fail to run command '$cmd'");
}

sub sam2vcf
{
	Transfomers(27);
	die "Usage: \nTransformers sam2vcf <IN:aln.sam> <-p:output_prefix>\n" if ((@ARGV==0)||(!defined $Prefix));
	open (IN,$ARGV[0]) || die $!;
	while(<IN>)
	{
		
	}
	close IN;
}

sub check_fafq
{
	my $file=shift;
	my $answer=0;
	open (IN,$file) || die $!;
	$_=<IN>;
	if ($_=~/^\>/)
	{
		<IN>;
		$_=<IN>;
		if ($_=~/^\>/)
		{
			$answer=1;
		}
		else
		{
			die "\nunidentified format of input sequence!\n";
		}
	}
	elsif ($_=~/^\@/)
	{
		<IN>;<IN>;<IN>;
		$_=<IN>;
		
		if ($_=~/^\@/)
		{
			$answer=2;
		}
		else
		{
			die "\nunidentified format of input sequence!\n";
		}
	}
	else
	{
		die "\nunidentified format of input sequence!\n";
	}
	close IN;
	return $answer;
}

sub check_454PE
{
	my ($name,$seq,$qual)=@_;
	$seq=uc($seq);
	my $head=$1 if ($name=~/(\S+)/);
	my ($pe1,$pe2);
	my ($strand1,$strand2)=split '',$PE_strand;
	if ((defined $seq)&&(defined $qual))
	{
		my ($qual_left,$qual_right)=("","");
		if ($seq=~/(\w+)$titanium_r(\w+)/)
		{
			my ($left,$right)=($1,$2);
			if ((length($left)>=20)&&(length($right)>=20)&&($right!~/$titanium_f/)&&($right!~/$titanium_r/))
			{
				RC(\$right) if ($strand1 ne "F");
				RC(\$left) if ($strand2 ne "R");
				#@GKZCRQK02J02YM rank=0000912 x=3995.0 y=3244.0 length=56
				if (($Format!~/\d+/)&&($name=~/\S+\s+rank\=(\S+)\s+x\=(\S+)\s+y\=(\S+)\s+length\=(\d+)/))
				{
					my ($a,$b,$c,$d)=($1,$2,$3,$4);
					if ($Format=~/Solexa/i)
					{
						$pe1="\@$head:$a:$b:$c:$d\/1\n$right\n+\n";
						$pe2="\@$head:$a:$b:$c:$d\/2\n$left\n+\n";
					}
					else
					{
						$left=base2col($left);
						$right=base2col($right);
						$pe1="\@$head\_$a\_$b\_$c\_$d\_$strand1"."3\n$right\n+\n";
						$pe2="\@$head\_$a\_$b\_$c\_$d\_$strand2"."3\n$left\n+\n";
					}
				}
				else
				{
					$pe1="\@$head$strand1 template=$head dir=$strand1 library=pairlib\n".$right."\n+\n";
					$pe2="\@$head$strand2 template=$head dir=$strand2 library=pairlib\n".$left."\n+\n";
				}
				my @q=();
				if ($qual=~/\d+\s+\d+/)
				{
					@q=split /\s+/,$qual;
					for (my $i=0;$i<@q;$i++)
					{
						if ($i<length($left))
						{
							$qual_left.=$conv_table[$q[$i]+64];
						}
						if ($i>=length($left)+length($titanium_r))
						{
							$qual_right.=$conv_table[$q[$i]+64];
						}
					}
				}
				else
				{
					@q=split "",$qual;
					for (my $i=0;$i<@q;$i++)
					{
						if ($i<length($left))
						{
							$qual_left.=$q[$i];
						}
						if ($i>=length($left)+length($titanium_r))
						{
							$qual_right.=$q[$i];
						}
					}
				}
				$qual_right=reverse($qual_right) if ($strand1 ne "F");
				$qual_left=reverse($qual_left) if ($strand2 ne "R");
				$pe1.=$qual_right."\n";
				$pe2.=$qual_left."\n";
				return ($pe1,$pe2);
			}
			else
			{
				my @q=();
				my $q_out="";
				if ($qual=~/\d+\s+\d+/)
				{
					@q=split /\s+/,$qual;
					map{$q_out.=$conv_table[$_+64];}@q;
				}
				else
				{
					$q_out=$qual;
				}
				$pe1="\@$name\n$seq\n+\n$q_out\n";
				return $pe1;
			}
		}
		elsif ($seq=~/(\w+)$titanium_f(\w+)/)
		{
			my ($left,$right)=($1,$2);
			if ((length($left)>=20)&&(length($right)>=20)&&($right!~/$titanium_f/)&&($right!~/$titanium_r/))
			{
				RC(\$left) unless ($strand1 ne "F");
				RC(\$right) unless ($strand2 ne "R");
				if (($Format!~/\d+/)&&($name=~/\S+\s+rand\=(\S+)\s+x\=(\S+)\s+y\=(\S+)\s+length\=(\d+)/))
				{
					if ($Format=~/Solexa/i)
					{
						$pe1="\@$head:$1:$2:$3:$4\/1\n$left\n+\n";
						$pe2="\@$head:$1:$2:$3:$4\/2\n$right\n+\n";
					}
					else
					{
						$left=base2col($left);
						$right=base2col($right);
						$pe1="\@$head\_$1\_$2\_$3\_$4\_$strand1"."3\n$left\n+\n";
						$pe2="\@$head\_$1\_$2\_$3\_$4\_$strand2"."3\n$right\n+\n";
					}
				}
				else
				{
					$pe1="\@$head$strand1 template=$head dir=$strand1 library=pairlib\n".$left."\n+\n";
					$pe2="\@$head$strand2 template=$head dir=$strand2 library=pairlib\n".$right."\n+\n";
				}
				my @q=();
				if ($qual=~/\d+\s+\d+/)
				{
					@q=split /\s+/,$qual;
					for (my $i=0;$i<@q;$i++)
					{
						if ($i<length($left))
						{
							$qual_left.=$conv_table[$q[$i]+64];
						}
						if ($i>=length($left)+length($titanium_r))
						{
							$qual_right.=$conv_table[$q[$i]+64];
						}
					}
				}
				else
				{
					@q=split "",$qual;
					for (my $i=0;$i<@q;$i++)
					{
						if ($i<length($left))
						{
							$qual_left.=$q[$i];
						}
						if ($i>=length($left)+length($titanium_r))
						{
							$qual_right.=$q[$i];
						}
					}
				}
				$qual_left=reverse($qual_left) unless ($strand1 ne "F");
				$qual_right=reverse($qual_right) unless ($strand2 ne "R");
				$pe1.=$qual_left."\n";
				$pe2.=$qual_right."\n";
				return ($pe1,$pe2);
			}
			else
			{
				my @q=();
				my $q_out="";
				if ($qual=~/\d+\s+\d+/)
				{
					@q=split /\s+/,$qual;
					map{$q_out.=$conv_table[$_+64];}@q;
				}
				else
				{
					$q_out=$qual;
				}
				$pe1="\@$name\n$seq\n+\n$q_out\n";
				return $pe1;
			}
		}
		else
		{
			my @q=();
			my $q_out="";
			if ($qual=~/\d+\s+\d+/)
			{
				@q=split /\s+/,$qual;
				map{$q_out.=$conv_table[$_+64];}@q;
			}
			else
			{
				$q_out=$qual;
			}
			$pe1="\@$name\n$seq\n+\n$q_out\n";
			return $pe1;
		}
	}
	elsif (defined $seq)
	{
		if ($seq=~/(\w+)$titanium_r(\w+)/)
		{
			my ($left,$right)=($1,$2);
			if ((length($left)>=20)&&(length($right)>=20)&&($2!~/$titanium_f/)&&($right!~/$titanium_r/))
			{
				RC(\$right) if ($strand1 ne "F");
				RC(\$left) if ($strand2 ne "R");
				if (($Format!~/\d+/)&&($name=~/\S+\s+rand\=(\S+)\s+x\=(\S+)\s+y\=(\S+)\s+length\=(\d+)/))
				{
					if ($Format=~/Solexa/i)
					{
						$pe1="\>$head:$1:$2:$3:$4\/1\n$right\n";
						$pe2="\>$head:$1:$2:$3:$4\/2\n$left\n";
					}
					else
					{
						$left=base2col($left);
						$right=base2col($right);
						$pe1="\>$head\_$1\_$2\_$3\_$4\_$strand1"."3\n$right\n+\n";
						$pe2="\>$head\_$1\_$2\_$3\_$4\_$strand2"."3\n$left\n+\n";
					}
				}
				else
				{
					$pe1="\>$head$strand1 template=$head dir=$strand1 library=pairlib\n$right\n";
					$pe2="\>$head$strand2 template=$head dir=$strand2 library=pairlib\n$left\n";
				}
				return ($pe1,$pe2);
			}
			else
			{
				$pe1="\>$name\n$seq\n";
				return $pe1;
			}
		}
		elsif ($seq=~/(\w+)$titanium_f(\w+)/)
		{
			my ($left,$right)=($1,$2);
			if ((length($left)>=20)&&(length($right)>=20)&&($2!~/$titanium_f/)&&($right!~/$titanium_r/))
			{
				RC(\$left) unless ($strand1 ne "F");
				RC(\$right) unless ($strand2 ne "R");
				if (($Format=~/Solexa/i)&&($name=~/\S+\s+rand\=(\S+)\s+x\=(\S+)\s+y\=(\S+)\s+length\=(\d+)/))
				{
					if ($Format=~/Solexa/i)
					{
						$pe1="\>$head:$1:$2:$3:$4\/1\n$left\n";
						$pe2="\>$head:$1:$2:$3:$4\/2\n$right\n";
					}
					else
					{
						$left=base2col($left);
						$right=base2col($right);
						$pe1="\>$head\_$1\_$2\_$3\_$4\_$strand1"."3\n$left\n+\n";
						$pe2="\>$head\_$1\_$2\_$3\_$4\_$strand2"."3\n$right\n+\n";
					}
				}
				else
				{
					$pe1="\>$head$strand1 template=$head dir=$strand1 library=pairlib\n$left\n";
					$pe2="\>$head$strand2 template=$head dir=$strand2 library=pairlib\n$right\n";
				}
				return ($pe1,$pe2);
			}
			else
			{
				$pe1="\>$name\n$seq\n";
				return $pe1;
			}
		}
		else
		{
			$pe1="\>$name\n$seq\n";
			return $pe1;
		}
	}
}

sub RC
{
	my $str=shift;
	$$str=~s/^\s+//;
	$$str=~s/\s+$//;
	if ($$str=~/^\d+\s+\d+/)
	{
		my @t=split;
		for (my $i=@t-1;$i>=0;$i--)
		{
			$$str.="$t[$i] ";
		}
		chomp $$str;
	}
	else
	{
		if ($$str=~/^[ACGTN]+$/)
		{
			$$str =~ tr/ACGTacgt/TGCAtgca/ ;
		}
		$$str=reverse($$str);
	}
}

sub cds2aa
{
	Transformers(25);
	$Transet ||= 'nuclear';
	die("Usage: \nTransformers cds2aa <-S Seq.fa> <-p:output_prefix> [-translate nuclear]\n") if ((!defined $Seq)||(!defined $Prefix));
	warn("-- The default translate table is set to $Transet. Use '-translate' at the command line to change the default.\n");
	my $output=$Prefix.".aa";
	open (IN,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
	open (OUT,">$output") || die "Fail to create FASTQ file:$output for writing\n";
	my ($head,$seq,$phase)=("","",0);
	while (<IN>)
	{
		if (/^>(.*)/)
		{
			if (length($seq)>0)
			{
				RC(\$seq) if (defined $RC);
				my $prot = na2aa($seq,$phase,$CODE{$Transet});
				if ((defined $prot)&&(length($prot)>0))
				{
					Display_str(\$prot,100);
					print OUT ">$head [translate_table: $Transet]\n".$prot;
				}
			}
			$seq = "";
			$head = $1;
			$phase = ($head =~ /\s+phase[:\s]+([012])\s+/i) ? $1 : 0 ;
		}
		else
		{
			s/\s+$//;
			$seq.=$_;
		}
	}
	if (length($seq)>0)
	{
		RC(\$seq) if (defined $RC);
		my $prot = na2aa($seq,$phase,$CODE{$Transet});
		if ((defined $prot)&&(length($prot)>0))
		{
			Display_str(\$prot,100);
			print OUT ">$head [translate_table: $Transet]\n".$prot;
		}
	}
	close IN;
	close OUT;

}

sub na2aa
{
	my $seq = shift;
	my $phase = shift || 0;
	my $translate_p = shift;
	
	$seq =~ s/\s//g;
	$seq = uc($seq);
	
	my $len = length($seq);
	
	my $prot;
	for (my $i=$phase; $i<$len; $i+=3) {
		my $codon = substr($seq,$i,3);
		last if(length($codon) < 3);
		$prot .= (exists $translate_p->{$codon}) ? $translate_p->{$codon} : 'X';
	}
	$prot =~ s/U$// if (defined $prot);
	$prot .= "*";
	return $prot;

}

sub sort2pe
{
	Transformers(20);
	die "No inputfile!\nExample: Transformers sort2pe <in_reads1> <in_reads2> -p out [-t fq|fa] [-o 0|1] [-len 30]\n-t  output format, default: fa\n-o  output format, 0 will shuffle the reads into *.pair file, 1 will keep in single format, default: 0\n-len is used to defined the theshold length \(the shortest\) of reads\n" if ((@ARGV<2)||(!defined $Prefix));
	my $outformat = (defined $Outset) ? $Outset : 0;
	my $type = (defined $Typeset) ? $Typeset : 'fa';
	my $len_th= (defined $Lenset) ? $Lenset : 30;
	#print STDERR "-len set value of: $len_th\n";
	my ($in1,$in2)= @ARGV;
	open (IN2,$in2) || die "Can't open $in2 for reading\n";
	my %hash=();
	while(<IN2>)
	{
		if ((/^\@(\S+)\/[12]/)||(/^\@(\S+)\s+\S+/))
		{
			my $id=$1;
			$_=~s/^\@/\>/ if ($type eq "fa");
			my $name=$1 if ($_=~/^(\S+)/);
			$_=<IN2>;
			chomp;
			if (($_ ne "")&&(length($_)>=$len_th))
			{
				$hash{$id}="$name\n";
				$hash{$id}.="$_\n";
				if ($type eq "fq")
				{
					$hash{$id}.=<IN2>;
					$hash{$id}.=<IN2>;
				}
			}
		}
		elsif ((/^\>(\S+)\/[12]/)||(/^\>(\S+)\s+\S+/))
		{
			my $name=$1 if ($_=~/^(\S+)/);
			$_=<IN2>;
			chomp;
			if (($_ ne "")&&(length($_)>=$len_th))
			{
				$hash{$1}="$name\n";
				$hash{$1}.="$_\n";
			}
		}
	}
	close IN2;
	open (IN1,$in1) || die "Can't open $in1 for reading\n";
	if ($outformat==0)
	{
		open (PE,">$Prefix.pair");
	}
	else
	{
		open (OUT1,">$Prefix\_PE1.$type");
		open (OUT2,">$Prefix\_PE2.$type");
	}
	open (SE,">$Prefix.single");
	while(<IN1>)
	{
		if ((/^\@(\S+)\/[12]/)||(/^\@(\S+)\s+\S+/))
		{
			my $id=$1;
			my $name=$1 if ($_=~/^(\S+)/);
			$_=<IN1>;
			chomp;
			if (exists $hash{$id})
			{
				if ($outformat == 0)
				{
					if ($type eq "fq")
					{
						if (($_ ne "")&&(length($_)>=$len_th))
						{
							print PE "$name\n";
							print PE "$_\n";
							$_=<IN1>;
							print PE $_;
							$_=<IN1>;
							print PE $_;
							print PE $hash{$id};
							undef $hash{$id};
							delete $hash{$id};
						}
						else
						{
							print SE $hash{$id};
							undef $hash{$id};
							delete $hash{$id};
						}
					}
					else
					{
						if (($_ ne "")&&(length($_)>=$len_th))
						{
							$name=~s/^\@/\>/;
							print PE "$name\n";
							print PE "$_\n";
							print PE $hash{$id};
							undef $hash{$id};
							delete $hash{$id};
						}
						else
						{
							print SE $hash{$id};
							undef $hash{$id};
							delete $hash{$id};
						}
					}
				}
				else
				{
					if ($type eq "fq")
					{
						if (($_ ne "")&&(length($_)>=$len_th))
						{
							print OUT1 "$name\n";
							print OUT1 "$_\n";
							$_=<IN1>;
							print OUT1 $_;
							$_=<IN1>;
							print OUT1 $_;
						}
						print OUT2 $hash{$id};
						undef $hash{$id};
						delete $hash{$id};
					}
					else
					{
						if (($_ ne "")&&(length($_)>=$len_th))
						{
							$name=~s/^\@/\>/;
							print OUT1 "$name\n";
							print OUT1 "$_\n";
						}
						print OUT2 $hash{$id};
						undef $hash{$id};
						delete $hash{$id};
					}
				}
			}
			else
			{
				if ($type eq "fq")
				{
					if (($_ ne "")&&(length($_)>=$len_th))
					{
						print SE "$name\n";
						print SE "$_\n";
						$_=<IN1>;
						print SE $_;
						$_=<IN1>;
						print SE $_;
					}
				}
				else
				{
					if (($_ ne "")&&(length($_)>=$len_th))
					{
						$name=~s/\@/\>/;
						print SE "$name\n";
						print SE "$_\n";
					}
				}
			}
		}
		elsif ((/^\>(\S+)\/[12]/)||(/^\>(\S+)\s+\S+/))
		{
			my $id=$1;
			my $name=$1 if ($_=~/^(\S+)/);
			$_=<IN1>;
			chomp;
			if (exists $hash{$id})
			{
				if ($outformat == 0)
				{
					if (($_ ne "")&&(length($_)>=$len_th))
					{
						print PE "$name\n";
						print PE "$_\n";
						print PE $hash{$id};
						undef $hash{$id};
						delete $hash{$id};
					}
					else
					{
						print SE $hash{$id};
						undef $hash{$id};
						delete $hash{$id};
					}
				}
				else
				{
					if (($_ ne "")&&(length($_)>=$len_th))
					{
						print OUT1 "$name\n";
						print OUT1 "$_\n";
					}
					print OUT2 $hash{$id};
					undef $hash{$id};
					delete $hash{$id};
				}
			}
			else
			{
				if (($_ ne "")&&(length($_)>=$len_th))
				{
					print SE "$name\n";
					print SE "$_\n";
				}
			}
		}
	}
	close IN1;
	close PE if ($outformat==0);
	close OUT1 if ($outformat==1);
	close OUT2 if ($outformat==1);
	foreach my $k(keys %hash)
	{
		print SE $hash{$k} unless (($hash{$k} eq "")||($hash{$k} =~ /^\s+/));
	}
	close SE;
}

sub bed2gff
{
	die("Usage: \nTransformers bed2gff <infile.bed> [-p:output_prefix] [-f gene]\n") if (@ARGV==0);
	my $source = $ARGV[0];
	open (IN,$source) || die $!;
	open (OUT,">$Prefix.gff")  if (defined $Prefix);
	while(<IN>)
	{
		chomp;
		my @t=split /\s+/,$_;
		$t[10]=~s/\,$//;
		$t[11]=~s/\,$//;
		my @blockSizes=split /\,/,$t[10];
		my @blockStarts=split /\,/,$t[11];
		if (defined $Prefix)
		{
			print OUT ("$t[0]\t$source\tframe\t",($t[1]+1),"\t",($t[2]+1),"\t$t[4]\t$t[5]\t\.\tID=$t[3]\; Parent=$t[3]\n");
		}
		else
		{
			print ("$t[0]\t$source\tframe\t",($t[1]+1),"\t",($t[2]+1),"\t$t[4]\t$t[5]\t\.\tID=$t[3]\; Parent=$t[3]\n");
		}
		for (my $i=0;$i<$t[9];$i++)
		{
			if (defined $Prefix)
			{
				print OUT ("$t[0]\t$source\tblock\t",($t[1]+1+$blockStarts[$i]),"\t",($t[1]+$blockStarts[$i]+$blockSizes[$i]),"\t$t[4]\t$t[5]\t\.\tParent=$t[3]\n");
			}
			else
			{
				print ("$t[0]\t$source\tblock\t",($t[1]+1+$blockStarts[$i]),"\t",($t[1]+$blockStarts[$i]+$blockSizes[$i]),"\t$t[4]\t$t[5]\t\.\tParent=$t[3]\n");
			}
		}
	}
	close IN;
	close OUT if (defined $Prefix);
}

sub instruction
{
  print qq(
==> FASTQ <==
FASTQ format is first used in the Sanger Institute, and therefore
we take the Sanger specification as the standard FASTQ. Although
Solexa/Illumina reads file looks pretty much like the standard
FASTQ, they are different in that the qualities are scaled
differently. In the quality string, if you can see a character
with its ASCII code higher than 90, probably your file is in the
Solexa/Illumina format.

Sometimes we also use an integer, instead of a single character,
to explicitly show the qualities. In that case, negative
qualities indicates that Solexa/Illumina qualities are used.

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

SOLiD color int bases can be use a chain matrix code to convert reciprocally with fasta bases

decode:
"A0" => "A",
"C0" => "C",
"G0" => "G",
"T0" => "T",
"A1" => "C",
"C1" => "A",
"G2" => "A",
"A2" => "G",
"A3" => "T",
"T3" => "A",
"C2" => "T",
"T2" => "C",
"C3" => "G",
"G3" => "C",
"G1" => "T",
"T1" => "G",
"A4" => "N",
"A." => "N",
"A5" => "N",
"C4" => "N",
"C." => "N",
"C5" => "N",
"G4" => "N",
"G." => "N",
"G5" => "N",
"T4" => "N",
"T." => "N",
"T5" => "N",
"N5" => ( A C T or G ),
"N." => "N",
"N6" => "N",
"N0" => "N",
"N1" => "N",
"N2" => "N",
"N3" => "N",

code:
"AA" => 0,
"AC" => 1,
"AG" => 2,
"AT" => 3,
"CA" => 1,
"CC" => 0,
"CG" => 3,
"CT" => 2,
"GA" => 2,
"GC" => 3,
"GG" => 0,
"GT" => 1,
"TA" => 3,
"TC" => 2,
"TG" => 1,
"TT" => 0,

==> SAM <==
read_28833_29006_6945	99	chr20	28833	20	10M1D25M = 28993	195	AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG <<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<	NM:i:1 RG:Z:L1
read_28701_28881_323b	147	chr20	28834	30	35M	=	28701	-168	ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA	<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<	MF:i:18 RG:Z:L2

SAM is TAB-delimited. Apart from the header lines,  which  are  started
with the `\@' symbol, each alignment line consists of:

+----+-------+----------------------------------------------------------+
|Col | Field |                       Description                        |
+----+-------+----------------------------------------------------------+
| 1  | QNAME | Query (pair) NAME                                        |
| 2  | FLAG  | bitwise FLAG                                             |
| 3  | RNAME | Reference sequence NAME                                  |
| 4  | POS   | 1-based leftmost POSition/coordinate of clipped sequence |
| 5  | MAPQ  | MAPping Quality (Phred-scaled)                           |
| 6  | CIAGR | extended CIGAR string                                    |
| 7  | MRNM  | Mate Reference sequence NaMe (`=' if same as RNAME)      |
| 8  | MPOS  | 1-based Mate POSistion                                   |
| 9  | ISIZE | Inferred insert SIZE                                     |
|10  | SEQ   | query SEQuence on the same strand as the reference       |
|11  | QUAL  | query QUALity (ASCII-33 gives the Phred base quality)    |
|12  | OPT   | variable OPTional fields in the format TAG:VTYPE:VALUE   |
+----+-------+----------------------------------------------------------+

Each bit in the FLAG field is defined as:

+-------+--------------------------------------------------+
| Flag  |                   Description                    |
+-------+--------------------------------------------------+
|0x0001 | the read is paired in sequencing                 |
|0x0002 | the read is mapped in a proper pair              |
|0x0004 | the query sequence itself is unmapped            |
|0x0008 | the mate is unmapped                             |
|0x0010 | strand of the query (1 for reverse)              |
|0x0020 | strand of the mate                               |
|0x0040 | the read is the first read in a pair             |
|0x0080 | the read is the second read in a pair            |
|0x0100 | the alignment is not primary                     |
|0x0200 | the read fails platform/vendor quality checks    |
|0x0400 | the read is either a PCR or an optical duplicate |
+-------+--------------------------------------------------+

==> VCF <==
VCF is a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome.
There is an option whether to contain genotype information on samples for each position or not.
Example:

##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
20     1230237 .         T      .       47   PASS   NS=3;DP=13;AA=T                   GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
20     1234567 microsat1 GTCT   G,GTACT 50   PASS   NS=3;DP=9;AA=G                    GT:GQ:DP    0/1:35:4       0/2:17:2       1/1:40:3

The header line names the 8 fixed, mandatory columns. These columns are as follows:
#CHROM
POS
ID
REF
ALT
QUAL
FILTER
INFO

==> BED <==
chr17	42927836	42928156	n104	982	+	42927836	42928156	0	1	320,	0,
chr1	231948705	231956006	n1114	988	+	231948705	231956006	0	15	1137,837,765,13,23,1098,523,1164,62,9,39,916,196,175,324,	0,1141,1990,2756,2770,2793,3891,4414,5578,5640,5649,5688,6605,6801,6977,

The first three required BED fields are:
1.chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
2. chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
3. chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
The 9 additional optional BED fields are:
4.name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
5.score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray score in range  � 166, 167-277, 278-388, 389-499, 500-611, 612-722, 723-833, 834-944, � 945
6.strand - Defines the strand - either '+' or '-'.
7.thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
8.thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
9.itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
10.blockCount - The number of blocks (exons) in the BED line.
11.blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
12.blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.

==> GFF <==
Here is a brief description of the GFF fields:
1.seqname - The name of the sequence. Must be a chromosome or scaffold.
2.source - The program that generated this feature.
3.feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
4.start - The starting position of the feature in the sequence. The first base is numbered 1.
5.end - The ending position of the feature (inclusive).
6.score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
7.strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
8.frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
9.group - All lines with the same group are linked together into a single item.

);

}

sub example {
  my $exam_scarf = '
USI-EAS50_1:4:2:710:120:GTCAAAGTAATAATAGGAGATTTGAGCTATTT:23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 19 23 23 23 18 23 23 23
USI-EAS50_1:4:2:690:87:GTTTTTTTTTTTCTTTCCATTAATTTCCCTTT:23 23 23 23 23 23 23 23 23 23 23 23 12 23 23 23 23 23 16 23 23 9 18 23 23 23 12 23 18 23 23 23
USI-EAS50_1:4:2:709:32:GAGAAGTCAAACCTGTGTTAGAAATTTTATAC:23 23 23 23 23 23 23 23 20 23 23 23 23 23 23 23 23 23 23 23 23 12 23 18 23 23 23 23 23 23 23 23
USI-EAS50_1:4:2:886:890:GCTTATTTAAAAATTTACTTGGGGTTGTCTTT:23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23
USI-EAS50_1:4:2:682:91:GGGTTTCTAGACTAAAGGGATTTAACAAGTTT:23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 20 23 23 23 23 23 23 23 23 23 23 23 18 23 23 23 23
USI-EAS50_1:4:2:663:928:GAATTTGTTTGAAGAGTGTCATGGTCAGATCT:23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23
';

  my $exam_fqint = '
@4_1_912_360
AAGGGGCTAGAGAAACACGTAATGAAGGGAGGACTC
+4_1_912_360
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 21 40 40 40 40 40 40 40 40 40 26 40 40 14 39 40 40
@4_1_54_483
TAATAAATGTGCTTCCTTGATGCATGTGCTATGATT
+4_1_54_483
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 16 40 40 40 28 40 40 40 40 40 40 16 40 40 5 40 40
@4_1_537_334
ATTGATGATGCTGTGCACCTAGCAAGAAGTTGCATA
+4_1_537_334
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 21 29 40 40 33 40 40 33 40 40 33 31 40 40 40 40 18 26 40 -2
@4_1_920_361
AACGGCACAATCCAGGTTGATGCCTACGGCGGGTAC
+4_1_920_361
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 39 40 40 40 40 40 40 40 40 31 40 40 40 40 40 40 15 5 -1 3
@4_1_784_155
AATGCATGCTTCGAATGGCATTCTCTTCAATCACGA
+4_1_784_155
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 31 40 40 40 40 40
@4_1_595_150
AAAGACGTGGCCAGATGGGTGGCCAAGTGCCCGACT
+4_1_595_150
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 30 40 40 40 40 40 40 40 40 40 20 40 40 40 40 40 14 40 40
';

  my $exam_solexa = '
@SLXA-B3_649_FC8437_R1_1_1_610_79
GATGTGCAATACCTTTGTAGAGGAA
+SLXA-B3_649_FC8437_R1_1_1_610_79
YYYYYYYYYYYYYYYYYYWYWYYSU
@SLXA-B3_649_FC8437_R1_1_1_397_389
GGTTTGAGAAAGAGAAATGAGATAA
+SLXA-B3_649_FC8437_R1_1_1_397_389
YYYYYYYYYWYYYYWWYYYWYWYWW
@SLXA-B3_649_FC8437_R1_1_1_850_123
GAGGGTGTTGATCATGATGATGGCG
+SLXA-B3_649_FC8437_R1_1_1_850_123
YYYYYYYYYYYYYWYYWYYSYYYSY
@SLXA-B3_649_FC8437_R1_1_1_362_549
GGAAACAAAGTTTTTCTCAACATAG
+SLXA-B3_649_FC8437_R1_1_1_362_549
YYYYYYYYYYYYYYYYYYWWWWYWY
@SLXA-B3_649_FC8437_R1_1_1_183_714
GTATTATTTAATGGCATACACTCAA
+SLXA-B3_649_FC8437_R1_1_1_183_714
YYYYYYYYYYWYYYYWYWWUWWWQQ
';
  my $exam_csfasta = '
# Title: BARB_20071114_2_YorubanMP-BC3_
>3_16_150_F3
T0220100010131232212020122
>3_16_157_F3
T3321200210021113013033210
';

  my $titanium_linkers = '
>titanium_r
TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG
>titanium_f
CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA  
';

  print qq(
solexa
======
$exam_solexa
fqint
=====
$exam_fqint
scarf
=====
$exam_scarf
csfasta
======
$exam_csfasta
);
}
sub Transformers
{
	my $TF=shift;
	if ($TF eq "BENM")
	{
my $Transformers = <<TRANSFORMERS;
##########################################################################
#  Copyright (c) 2009-2013 BENM(Binxiao) Feng                            #
#  All Rights Reserved                                                   #
#  Send all comments to BENM - BinxiaoFeng\@gmail.com                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  (at your option) any later version.                                   #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#  GNU General Public License for more details.                          #
#  You should have received a copy of the GNU General Public License     #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. #
##########################################################################

                          \@
                        \@            \@\@\@\@\@\@\@\@\@\@\@
                      \@\@      \@\@NNNW&YyAB&V\@\@
                    B7\@   \@AQcYyVyXXAXWH\@N\@\@\@\@\@\@\@\@\@\@
                 \@2YN\@\@7vy&Y&VXX&ccvcYVYVYYVVY&&XXAXX&AXQN\@\@
                Q\%VQ\@<vYVVXXQVvcvy&VV&YV&XXQQAAXAAQH\@\@\@
               \@\(&&Av\%VAAXAX2&VYYYYYYYYVYYYYYYYYYVAW&\@
               N7yAYWXVv777y&XXXYV&&YYYYYY&&YY&VcvvYXXA\@
              \@A7AbvAVcyYV&\%v\(\%YXYYy\%Y&YVYyYcYV2cyVV&y&Q&\@
            \@WVAvNNX&vv\(\(\(\%YYcc7vyAy2\(YYY&VYV\%\%Y&V&YYYYY&&N
          \@H7V\@YvcYv\%c\(v\%\%y<cvyV&7yVY\%2\(Y&VY&&&27YVYYYYYYVyA\@
          VcYYWy\(y2v\(vV<&YYYV\(&VYyV\%YYYV\(VVYY&&Y&2\%V&X&YV\%XVYA\@
         N&X7YXAv\(Yv&\(Vy<7yYVV\(YVYYYy\%&&&\%VYVY&YVYYyyY&XAXXVyYYX\@
        \@\%X\(YVy&Y77YV&7YAv\(VV&&7\%YYYYY\%Y&&XY&yYYVVY&X&YV&&AX&X&VV\@
        VW/7V2VHX&\(YV2ycWY\%vY&&&Yy&V&YV&YAVX&YV&VY&VYYX&VVYA\@  \@\@QW\@
       \@X&7&\%\(YQHWV\(Y2V\%XWVycYX&Y\%&VVYVYVAAXcV&XVYX\%&A&XXV&V&\@      \@\@
       \@AAYy&v&W\@\@XyVVVVY&QYVc&WXVV&AYYVYAXX2VVXYY&yYXYVVWYYYY\@
      \@B\@v&y&vyAWc\@A&&&VVABQ&YYbXXV&yAYYYXAWX&VX&&VVV&XY&&AWY\%A\@
      \@\@ &XVV<yQH7vQAXvVVAX\@WXY&VQWX&YAYYAAAWXYXAVYVVYXV&VWW\@QY&\@
      \@  YXV&7vWYv7vWXQV&XQ\@&HNX&YW\@A&XQYXAQB&&XY&YYVYAXY\%BN\@ \@Q\%
      \@  QA&XA7A\%v</v7&WAAXWN77HBXXAWBXWQ&AAQWYA&&XQYYWA&&ANN   \@\@
         \@WQvQYAAXV\%</cWQ&XQQ\%\(7v\@QXXW\@\@Q&AAAWYXVVX&Y&XAA&XNN     \@
         \@YWyXXA\@X\@\@&X/ \%\@XXAH7 <\%AXHQWXWNAQQbXXVYAYVXXHHYYN\@
          W&WXQX\@y\@\@y77////\@WHc/<X\@\@\@\@V/XWW\@N&XYV&AYVAB\@\@QyHB
          \@XQHXB7\%/\@\( 7/   / /B/7/\@\@\@NV <A7v\@QXY&Q\@YVQB   XB\@
           \@B\@\@XW\@v  \(y        < 7      A/<v&WV&XWAyXQ\@    W\@
            \@ \@A\%7A7V<           /7XX&&&//;WAVYXWXWYA\@\@    \@\@
              \@WV//    //             /// \%7y&AH&N\@WBH      \@
             \@ \@Vv</ /                // v7vVA\@A/\@ \@\@\@
               \@\@  \@/               ////\(\(;&&AW/// v\@
               \@      \@7 //\(      //\(\(\(\(\(/&Q\(  //  X\(\@
               \@ < /\(/v /7<7;\(\(<\(\(\(\(\(\(\(\(/V</         /N
                \@/   /7/\(\(/ \(v///\(\(\(\(/   /
               //  /\(<\(/\(7                       /     //
               \@    \(7\(/\(\(///   //  /      /< / //     /<\@
              \@   /\(\(\(//7v//               </ </ /     \( <
              \@//<\(\(//7\(\(\(7//    /<<///<<<// /\(7;      <//\@
                <////y7\(\(\(;\(/   //////;< / /<v\(       /7/<\@
               \( \(7\(<<//;/\(\(<;< ///7 / //<//<        /7<///
              \@ y<</    //<<;\%\(/\(<;\(/\(/<7<\(< <       /\(// \(
             \@7/\(/       /<\(/c<///////<<7\(;\(v /     /c\(/  <
             \@</\(/        /<\(v/      ///;77\(;\(y/7</ vv</  <
              \@<</        /<\(\(/      / /<\(\(;</;V/<7;;//  /N
               \@c\(      //<<;\%<</   / /<;\(v\(\(;/<;Yv\(///   \@
                 \@y/     /<<\(7;/<    ///<\(\(7\(<</<;Y<</   /\@
                   7/ /    /<\(7</    / / <<2\(\(<<<<\(<//   /\@
                    \@v/     /<c\(</   / //<<v\(\(<<<<7<;/  //
                      /      /;2\(/   ////<<\(\(\(\(<</\(\(/</    \@
                      \( / /<<<722<   /  /<;\(\(\(\(<<<\(\(// <\(/<
                      \@   /<\(\(<;<</   /<;/<<7\(\(\(\(\(7v2cv   v
                       \@v;7\(\(<<\(c\%// /\(\(<<<<\%77\(;\(cy
                       \@/  ////2<//<    //\(\(\%Y7cc\%\@
                       c</   v\(c\(;/ \(///<\(7\(77\(
                    7<<  //\(// <\(\(\%v //<;\(\(V\(
                    y\(7/\(7v\@  < // /  /  /\(\@
                            \@; \(/ \(/ // /\@
                             \@\@ \@\@ \@\@

                                  BENM
                         BinxiaoFeng\@gmail.com
                            1981.10-2079.12
                             Made in China
                           All rights reserved
TRANSFORMERS
	print STDERR $Transformers;
	exit();
	}
	if ($TF==0)
	{
my $Transformers0 = <<TRANSFORMERS;
                                                                    A
                AAA0000AAA               00      0         C       0
       00000A0000CCC1CCCC0000000000      AC0      C0      CA     00C
       A0000 000C0      0CCC00C000A       CCCC    01C   0CC     CCC0
        0C 0CC CCC11C01CCC000C0000        CCCCCCC10111111110CCCCCCC
        CC0C0 CCC01CC1CC CCC C000A        CCC101C1 111 C11A1110 CCC
        ACCC0C101CC0C000110C00CCC          C11111111 10111 C11CCCC0
         1011CC 11C1111 CCA101C00          C11C0 C 0111111CC  0C1CA
         A1C011C11 111C0C1110 1CA           1111C111     C11C1C11A
         A11C     111110      11           A CC0  A111CC11C  1CCA1
          1111C1111111C111111111            10CC111 11 111  1CCA1A
          C111111111111111C111CC            111111111111111111 11
          CCCC111C1111CC11A1CCCC            ACC1A111111111111 1CC
          ACCCC1C1111C11110CCCCA             CCC10C11C1CC111CCCC0
            ACCCC1     ACC0CCA               CC0C1ACCCC1CC1CCCA0A
              ACC0CCCCC00C0                  0C0CC1 CCCCC0CCCCCC
                 000AAA0A                        C0C CCCACCC0
                                                    AACAA


   AC111 11 111111  A 11   0 1 111A1111
      11 A 1   C    11A        111A
       1 0     1 A 1 1          110AAAA
      C1 A    A   11 C  1    1  11  1 11
      C1 0 AAA  AAAAAAAAAAAAAA AAAAAAAAAA
      C 11       C1   C1 1 1   1   1111 111C11111  C11111111111CA
             CAAA0  AA10A1 AA    11A01 A   1  AAA  11A0 1C11AAAA
                C   A   01   1      01  11 1C   C  11    01    1C
             1          A  A1       0   C     1    111 1 CCC1C 1
                     1  0  A 1      0   1      11 1 1AA  A1 1 11
             1
            AA                 C1111C1111111 1
                                ACCCC0C0CCCCCC1111111111CC1C
                                 A0ACA0A00CC AA01111AAAAAAC
                                        01CC   A1111
                                        0CCC   A1111 111C
                                        CCC1   AC1C1AAA0
                                        CC11   A 111
                                        CCC1   A1  1
                                        CCCC   ACCC1
                                        ACCC   A C1C
                                         ACC   A 10
                                          A1   A1C
                                           AC  AC
                                            A  C
TRANSFORMERS
		print STDERR $Transformers0;
	}
	elsif ($TF==1)
	{
my $Transformers1 = <<TRANSFORMERS;
                                          11C 1
                                    11 1 CC1A1C1
                           1        CAA01 CA 111
                          C  11 11 1  CCAAAA0011
                         1 AACCCC1CCA101A 00A1
                00 A1A   AA10A011A00AAA A1A1AC
                CCA11   A1101 1111C0AC    1
                1AAAC0A1  C  C 0CCCA11
                1A1CA0AC1 AC01AAAAA1C
                    1    ACAAA101000A
                          111  ACAAA1C
                        1AAA C111CC1AC
                      0A011ACAC0AC A111 1
                        A1CC01C    AA111C
                 110 AAAAA0A          1 1 1 A
                   10C1101 10         C  A01
                    1111AA 1            1 1A AA  C
                    1A11 00C          11AC AAC1AA11
                     CC 0C1           A10AAAAA1 A  A
                      C0A0A111            A   011  1
                     11 ACCAC0A1         0A1A11CCCC  1
                     1 C11C01  A1        AAC 110 CC0AC11111
                       1CAA  1 1 11  11  AA11   11 CA
                   1  1CCA11 1A11111   CAC1  1 C 1111 1
                                         1111  1A    111
                                                 1 C1
TRANSFORMERS
		print STDERR $Transformers1;
	}
	elsif($TF==2)
	{
my $Transformers2 = <<TRANSFORMERS;
                                      C1     1
                                      ACAC   1A
                                 1    A10A C1ACAC10C
                                ACACACC11C1CCA0000AA1
                             0C0A00CCCC0A1CC0 0CAACA
                            CCC A00C0AAAA0ACAA10 AA10CC
                              CAA 0CC00C 1AAA0    CCA0CC
                            CAA1C  CCAC111AA        C1C
                            AC0A1A0 1A0A 0A       AA000C
                           C0A1AC 0AC   AC10CC0   CAAAA
                       1A AA1C10 AA1C0AAAAC1 0A 0 ACC A
                      1AA  1A A CC1 11    CA1 C0AAAAC
                     A0A A0A0  A0CCC       A0100A
                    AA1 AAC10  AC00AC       AACAA
                   CA   CC A   AA0C          00AA
                  AA AC 1     1AA001        C0AAAA1
                    0A 01     AAAACC        CCCAAAC
                   0 C0       0A0AC1C       C0AACAC
                  C C         ACAAACC       A011AAAC
                 1 A         AAAA11AC       CC1A 1A0
                11A          1C ACAA0       CA0AAAA11
                 C          1C1110A11       CA0C11  10
                           AAAC1AA1 1        CAAAA A0C
                           0A1AAAAA          CAA0AAA1
                           C   A 0ACC11      AAAACA1A0
                         ACA1 1C010          00AA1 11A
                              1                01   AAA
                                                0C1
TRANSFORMERS
	print STDERR $Transformers2;
	}
	elsif($TF==3)
	{
my $Transformers3 = <<TRANSFORMERS;
                          1A   1
                       C  CC1 AA
                         1 A
                          C 0A1
                          A111        0       1
                            0A      CC  A1  1
                             A 1  ACA1A A1
                        A   1 A1CAA 1A0A AC   1
                AA       10A 0AA1C1   A   A110
                     AAAACA0C01A1CAA          A
                      A  C11A 11C  0
                         AA      AA01AA
                      A1 A       A1 A A
                              11 A     1  A
                              AAC C      AAA
                              1 A      AAAAA
                               AA       1A1
                             A1A         AA
                              A           AA 0
                             A             AC
                                            A1C
                             A                 AA
                            11                  AA
                           C1CA               1AACAA
                          A11 1C                 A 11A
                           11 1C1            1A C1111A0
                                                   AC
TRANSFORMERS
		print STDERR $Transformers3;
	}
	elsif($TF==4)
	{
my $Transformers4 = <<TRANSFORMERS;
                        1         10C
                      111       1 1AA
                       1      1 C A C1 1 ACC C
                      11 0     11AACCA  C CA A1
                        01 0ACAA AA0A0 C001ACAA11C
                       1A0  A     A 10AA01A1 1 1AAA
                       CAA1C      11C0AA A0CCA   01A
                     111C        1C AC00C11ACA 1CA01   1 11
                       1  1A   10 1C1 0AA1 0 A1  A0C000
                        C011     AA ACC0A1 A0     0AA0C1AAC
                                   A11 CC1 101C      AA1C0 C
                               1 11CAAA0A1  ACAA      0A00C1
                              10 11 1   0C1CAAC11          1
                           A C A01         AC11 10       0
                          C  1AC             0A A AC1  A 1
                         01 A10  C          CACAAAA11
                        AA0ACC 1A            1AAAAA11
                    1  0111C1 C 11         011AC11C
                    1C111CC1  111           11  11AC
                    0CCAA 01A1                10A10A
                 C1 C11 111 CA               11CCAA1101
                  0A11111 1 0A               1AC1CA 0C0
                  AC   1A1AA10A A1         1C  AAAA111
                   1   A  111AA0 11       11C 01  1AA0A1
             11 1CA1    A1  111C1           11 1   11111
TRANSFORMERS
		print STDERR $Transformers4;
	}
	elsif($TF==5)
	{
my $Transformers5 = <<TRANSFORMERS;
                                        1
                                       1
                                     11C1
                                  0   011
                                  A1 A11 11
                                  1  0A
                          1      11A          1A111
                     0C11C AC00ACA1ACC      1CCC
                    1C  0C1C  A0  AA00A0 A1   A
                  10C0AACA       0A 1 0A1C10CAAA   1
                  0C 1CAA         A A0CAA0A11AA01C      1   1
                   C  A           1 00A11A00A0AAC   1C0A011
                  A1   1        1AA000A0C C1CAACAA0AA0C A1
                                   AAAAAA1C11 AC CACA0001
                            011AAAAAAACACA1C01CA    C1CAA
                            AA0CA00AAAAA0CA0AAC    11A
                            AC0CAA   AAA AC1CAC     A
                      1  A0 0AAA            CC1    1A
                     100AAAA0          001CCCC1CA1C1 1
                     1AA0A0 C       1A0 A AA AAA1AC1C
                     AAAA100     A  1AAAAAA1AACAAAAA
                      ACCA0     CA 00C1CCAA  ACAACAA
                       AAA        AAA AAACA AAA
                 A   11CA                   0010 0
                  0 C0A                    A0A11A
                1001AACCC                   CA1 A C
              1AACCC1CCCC         11111111111111C11
TRANSFORMERS
		print STDERR $Transformers5;
	}
	elsif($TF==6)
	{
my $Transformers6 = <<TRANSFORMERS;
                                             1
                                            1 0
                                           111C
                            1             ACC01101
                         C1C 0    C      1 AA  C1
                         C10101  C1 1C1 1C   00111
                          A11 0  C1 CAA  0C   1AA
                        A 0A0A   00 1111101A010A0
                       ACCA0A0AA0A0AAAA10A1ACAAAAAC
                    AAA 0C1CACAAA1AAA1C0A0AA  AAAA A
                   AACC0AA0C   0AA0CA0AACCC      11CC
                     0A0C1AA    AAAC1AA1AA0      AC1C1
                    0A0CC1 A0   CA100C11ACAA     0 A0C
                      C1 A1A1A00CCAAAAA00AAA     A 1AA0
                     11C  1  CACA1C     0 11     AA1AA1
                      1 0   01 CCA      1 C000    1AACCC
                       1     A111     1 00000   01 C00
                            0C11C     AC01A000     A 0
                            AA0CC 0    C00AACA   1 AA
                           A0AACC A    0A0AA1
                            AAAAAA0C   A00AA
                              A1AAAC    0AC11
                            C   C110    111CC
                             C1000C1    CC1C0
                             0000A1CA  CAC11A
                           10AAAA00   C010C0C1C
                                      A 00001CC
                                      1111CA0
TRANSFORMERS
		print STDERR $Transformers6;
	}
	elsif($TF==7)
	{
my $Transformers7 = <<TRANSFORMERS;
                       CC1                          CAA1AA0
                  1A C1AA01            1CC0    0A111 AACAA0
                11AAACA0AA1C         AC 0C   1CA01 CCA0CAAA
                 1A1CAA01  111   C   AAC0  111A0AA0100 AC1
                AA10C0AAA 11  A0 CC0 C10   A CA0CAAA01
                 1 C1CACAA11ACAACC00A00C0AAA00C 1C1AA
                     0C000AC100AA0AAA1C0AA0 CAA111C111
                       1CCC AA000CAC0 0A0AA1AAAAA CA1
                       0A1AAAAC1111CA1AA1 0 101AA
                         0C01C101CACCAC1
                        1CC01 1CA0CA1C
                       000CCA1C0ACC00A
                      0  C0A0C11CAC00A0CC
                        CAC0AA AAA0AC1AC0
                        AACAAAA 11CAC01CCC
                     CACA01AC  1 0CAA0CAC
                   10AAAAACC0    A10100C
                  CCC0C0AAAAC1 CAAACC1A0
                   CC0A111C0    AACC11AA0
                     C11C0 10   0AA0A1CAAC
                     0 A   A    0 01C0A
                       A   A    1AACA110C
                          111    C1A11CCAC1
                              111000CC0A1 C 0
                                01011A1A0A0C10
                               01111CCCAC1CC
                              0C001    C   0A1

TRANSFORMERS
		print STDERR $Transformers7;
	}
	elsif($TF==8)
	{
my $Transformers8 = <<TRANSFORMERS;
                          A             0
                         1C            AA
                        1CAA           C1
                        1AAAC     C   CC0
                       CCA AC     A AA0AA1
                       ACAACC001CCAA0ACAA1
                        CAACCA0001A0AAAA1
                         AC0A      CAAAA110A
                          0C0 C     0AAA01CAACC
                            A00     AA01AAA0AAA0AC0A0A
                                 0  ACA0           0A
                                 A0AAAA1AA1
                                AAAA01A0A0A1
                               0A00      0CAC
                               AA01       AAAC
                               CAC         A010
                              C1AA         AC00
                           0AACA         ACAAA
                          A0CA          10CC
                           AA            C1
                           1 A1           A0
                         0AAA A         AA0AA
                        1AAA0CCC1C1011C11AA AA
                        AA  AA000AA001011CA 1 0
                        CCC1111CC11111111111111
                        0CC1AC11CCC1111111C1111
                    1C1C11  11100CC01111111111CCC1
                    010CCC11CC0CC11A 11C11C111C101

TRANSFORMERS
		print STDERR $Transformers8;
	}
	elsif($TF==9)
	{
my $Transformers9 = <<TRANSFORMERS;
                                                A
                                             C  0 0
                                               CA00A0 C
                                             AAAAAA 00CC
                                          1C1AA0AA0    00CC
                                              0A   A0C   A C
                                            1CAAA   ACC0
                                                CA   AA0AC
                                                 1A   0AA0 C
                                                        AC
                  1C  10 1                              A1
            A    C 0000C10010                          CCC
        C    C C0  1ACC1CC AACC C1    A0CC           CCC
        1 C C0ACACA0A01CC01CAC0AA00A01A01CAAAA10AAAAA1C
        0AAC0110A0 AAAAAACA11 AA0CC01C 0CCC0A00AA
        C0CAAACAAA0AAAAAAA0AC1CAA0AC00 1CC0A 0A0C  C  1
       1AA1AAA1A10AAA0A01110001AC0C11 1CC01ACAA00C     C
        0AA001  AA0CA0C1 C CCAAA  1111CCA11  AAA0A0000  0C
       1 C0     AAAA0C1AC0AAA  C1CAC1CCCCCC C1A0A0A1A10011CA1
               AA00  0AA0AAA  11CCCC1AC 0 C1C0 AAAA0000AC0000
             0AACA1 110AC 0    A1CA0    011C0 A1  A      CC0A0
            0AAA0 1C00AAA         1    ACCA               0AA0C
        000AA00   CC0A          100C1CC                    0A0A0A
       A0CC00CCA0  0 0A        0A 1CA0CAAA0    1        0C0C0AAC0CA1
      0 0111AA111             0  AACA0AC  A0           A1A1C1AAA11AA
TRANSFORMERS
	print STDERR $Transformers9;
	}
	elsif($TF==10)
	{
my $Transformers10 = <<TRANSFORMERS;
                                 01CCC
               11           C  C 11C0CC   010
               10C11       C1C00C111CCCAACCC 00    1 C1
                CC  1   CA0C1000C0CC1CCA00C0ACCC1C 1C C0C
              AC0C1CC0C011ACA00CCCC00A0A1111C0101 0AC1CCCCAC
              1000AC1AAA0 CCCCCCC1C1AA00CAA1C11ACCC100A0A1 0
            1C00CA110A00CC0CCC10C C0A00A0CCC0C00AAAA00CA1 CC1
           C0CCC AACAAC0ACC1   C100AA0AA0C0C0CC0A0AAA11CC0000
           1A0CAAA0C 100C0C 11 1C0AAAAAAAAAACCA000AAA10A10AAA0C
            A00AAAC0 C0CA110AC ACAAAA0AAAAA0A    AAAAAAAA00C100C1
             1AA1A1  A0C0001000000AAA0A0A0A 0     0CAAA0AA1CC1 0AC
          111C1A0A0  0CCC11C00A000AAAA00CA        A0C0A 0AC 11CCA0C
        10 AAAA0AA   0C11CC0A      00AAAAC         AAAAA00C111CCC0C
       001AAA00C00   01C000         CAAA000       AAAAAAA1CC1C0AAA0
       0CC0CC111C                 AC   0C0        1000ACAC ACCA  C10
      0C1  0C1C                                       C    AACC    0C0
           A                                               1A01    C0C
                                                            A0A1
                                                             0A
TRANSFORMERS
		print STDERR $Transformers10;
	}
	elsif($TF==11)
	{
my $Transformers11 = <<TRANSFORMERS;
               C1CCC C11C111             CA0AA      C01
           11A1  A A10C A  C 1          CCCCC0 1 0C00C111 1C
            11CCAA1 CC1 01C00 C        0A0CCC111111CC0A1A0 1C
           111A00111 CACACC1011       C0C 1A0AC11CC0AAA0AC 0C
          CA0011AACCC10AAA 01CA1  1     C0C 0AACC111110A111100
        1 0AA111AA11C0A0AAAC010   1    C00 A00 0A1CAAACAA0 C CA
       110A1 1A111CCCCAA0A  111111   C1CC0    0CACCCC0A0AC11000C0 C
     11A A   CACCCCC1 A 0  1A01 C1 CC0A1       AAACCCAACCC    0ACC01
   11100C     0CA0C1CC0C   01C01   0C1C1      0CA0AA000C0       0ACC1
   CA11AA0   CA0CC1A0101C  1A1C    CCCC1      CCCCC000AAA 1      0ACC1
   C0CCC0   C AAA   00ACC          0C0ACC   C1C 0      A0AC1C    C0CAA1
    CCCA   1  A      110             1A0     11A        AC101C  001AA10
           1 CA      111                 C111CCA C     C0ACC01  1001A0
          01A0CAA    000CA0              CC1C0CAAC    0AAAAA1A
         A000AAA1CA AAA11CAC             CACAA0CC1  C1C000AACCC
         AAAAAAAAAA AA0CCCAC            0A0C0A0C10   00AA00CC0C
       C11AAA01CA0A1111CA0A0          CACAAACCAA A   CA0C00AACCCAA
    11C111 A0AAC1 000C1A0A0        111A0A0AA0A1CC111111AA00ACCC000AACCC
TRANSFORMERS
		print STDERR $Transformers11;
	}
	elsif($TF==12)
	{
my $Transformers12 = <<TRANSFORMERS;
            1                      1 C       1
          11  1            11     01  1      11
          11  111 11        11    C 0   1     00
              111111111111111111111 0 1C 111   AA1
                   11111111CC1111C101AA0 CCA1C 10A1
                    111 00111C0AAC000C00000CC1 111C
                 11 C1C0A0C0AAAAAAAAC00AC0CCA0C0  C11C0C1
            1111C01001110CC0A10ACAACAA0AAAAAAAC1 1C1C1111 1
          111101A 11C0CA0C10AA1AAAA A00A0C0CAAAA01  01C
        111111C11C 111 C000CAAAAAA0000A00AA00AAAAAA11C01CCC0001
        11CC01  10111CAA00A0AAA0AAA0CCCCC000C0AAA00AAACCCA00
       1111CA011111AAAAAA0AA111CCC1A1A0C10C1CAAA00A1   11   A
         11CC0A110ACA110C11111C0ACA0AC 11 C00AAA C100AAA 10 111
          111CA10ACAA1C    11CCC0ACA 11C1CC ACCAAAA11C00C CAA   1
             1ACCCC         CCC01CCA1A CA0 10A01 11111A0 A1C0A0C1C1
                             CC11110000AA10CCAC11111  AA1  CACA1  1
                                 11C001A00A0CC11111C1   AAA CCCCC
                                   11C1ACC00CCC0AAAAAAC01 A10AC
                                 1CCCCC1C0C1CC CAA1CCAAA0C   C
                                CCCC1 CCA 1CCC0C011CAACCAAC01
                              1CCC1011CCC01101C1CCACACCCA0C1
                             CCC1 1AA 101 1C1C1CCACCC1CCC
                            CCC1111C  1 1C111ACA010  A
                            CCCA0A1CA101111AC 1CC
                              00CC1CCC1111CAC11C
                                0C0C111C1C
                                  0CCC11
                                     C
TRANSFORMERS
		print STDERR $Transformers12;
	}
	elsif($TF==13)
	{
my $Transformers13 = <<TRANSFORMERS;
                           1C0
                          1CC0C
                         1C11AA  0C0      C11C
                        A11 C 1AC1A01A   AAC 1C0 C00
                    111AA1C CCC01AACCAAA1A011C111CC1   C
                 C1A11A1AA111CC 1AAC0AA  010A 10C1CC0C0C0
                 1CA11ACAC0C11C0 A100A1 1AAA CACA1CCAA0CA
              C11CA1A0CC0AA0AAACCACAA1A CCAAAA101C C  CC 00  0
             1C1CACACC11A  C011A C1A11 C1100CC         0C00A1 0
              ACA1110C0AC     C0  010 0100C1C0C           AC
              C0A 1AA1111         1C01 AA1110CC           AAC C
            1C00 CC0 11C           01011CCC1A0            C  1
           C1AC     11C            11A0000AA0A              C1
          CACC    CC             CA ACCC0C0000C101
          111                     1ACA1CA0ACC0CCAAC
          11C                   11 1AAA110 CA11C1
          11                    0A1AA0C C   01 11
          C                     1C10ACA     C1  1C 1
                               C1  1A0C    0CAA 111C0
                                11110CC    CCA001C1C
                                CCCC0AC1   A0C0AAA00 1
                                CC10A A   01AAAACC1CA0
                     111111111C0CAA1ACA1110C1CAAAA0C00
             11111111111  111 11AAAA1CC1111CCC0000C0111
           CC11   111           CCAAC0CCCCCCC01C0C0     1
          C                    C01AA00CCC1AA00AACC1
                              1 1A0A1A1  CCAAC10C0111
                               0C000CC    11111111
                                     1
TRANSFORMERS
		print STDERR $Transformers13;
	}
	elsif($TF==14)
	{
my $Transformers14 = <<TRANSFORMERS;
                                         CC
                                        1AC00C0A0    A00
                                      CAAC C1C0A111AC00A010
                                     1C A0AAC1A         A00 C1
                                    AAACAA1C0             A0C
                              11CAC0AAAAA0A1             C
                             ACACA01C A0CCA
                             1 1CA1A11CC000
                         0CAA01AA00C 1AAAAA
                       1CCA00AAA10 C0C10AAA
                     AC0C00AC00C0CA0CC1AACA    0
                     0C0A01CC C AA0AAACCA    C0CCC
                    CCCAC C1 0AA1AA00AA0A  CA0CA1C
                       0CA01111 AAAAA0AC00CA1A1CA0
                       A1C0CAA  C0 1C11AA111A0
                       AA10A10CC00C00111CAA0C
                       AA10C111AC001101 ACA00C
                      AA1CCAA0000100AA00AA0CC0CA
                    1CA11A0A              CC0100
                  A01C0A CA              0CCAA00111111 1      11111
              0A1AA0A1 CA          111111C0AA0111
              A C1CACC0    1111111111111CC1010110C0CCCC
            11A C1C        1CCCCCCCCCCCCC1 A000A0CCC    C
           1A0AAAA       1CCCCC1 CC111       0ACC010C 10C
         AA1A0A0  CC0000000                CAAAAACC1CC0CC
         CC0A0A00AAC                                    A0
        11 CAAA
       0    A
TRANSFORMERS
		print STDERR $Transformers14;
	}
	elsif($TF==15)
	{
my $Transformers15 = <<TRANSFORMERS;
                                                     A
                                                     T
                                                     GAC      C
                                                     GAG TAC
                                                      AGA
                                      CC       ATG TGCCTT
                                     CAAAC  A T TTCC  G C
                                     TACGCGGG CAAAT CAGCC
                          CGG  T C   T AGTTACC GGCCAATCC
                           C ACCCAGGCAAATACGCAGGCC
                            TCG CGAA AGCCTGCAAG C
            C             AACGACCA GCCATTGTTCGGC
              CGAAGC    AC ACC CTTACAAACAACGCA C
                  GGGGGCCACT       C   CCA TCCCG
                   GAACTCT C        C CTCAC AATC AC
                      CGATAAC     C   CGCGAA ACCTAC
                     TCC     CC   AAACA CAAAC ACC TC A
                  CT   CC      GGGAC GT CC     CC A AGT
                C     T G   GTAACACGCC          GCAG CCA
                             CCACAAC C          CCTCCCCC
                           T   TGCCAC            ACCGC
                           T  TA CG C           CGTGTCC
                              CCTCCC           CGCAGCC
                            CA GC ACCC         CCCCACAC
                           TCTAAG AGC C         CCCCGC CC
                      CCT C ACCC AGT            GCACCC ACC
                      CGG  AA AAAAAAGAC         G AG  CGT C
                     ATAAGATGTCTAA   ATG       ACGAAAAACA C
                    CA  ATGCGAAAAAAGGGTCCCCCCCCCAAAAAAAAAGC
TRANSFORMERS
		print STDERR $Transformers15;
	}
	elsif($TF==16)
	{
my $Transformers16 = <<TRANSFORMERS;
                                                C
                             CC                CC
                           C   CC               C
                          ATC AA CC AGCGCAA    AAC C
                            AAGAC A AGGACCA         G A
                           A CAAAGTGATCGAACCC      ACC   GAA CA
                            CA TC CCAGCGT GA     GCA A      C
                             ACC   ATCAGTCC       CC CAC   A
                       TCC     A    GG GCCACC C   A   C
                         AGC                 ACC      CC
                             TAGCAA       CAAA       CGG
                        GA    CAAA      CC AT  G
                              GA     CCC GAATCCC      C
                   CG     C  A    C      CAC A  AA A   A
                     CCAA            C       AAAGC AAAT
                      A            C C            AAAA  A
                   CGGC             C               AA
                                     CC            CGA
                                   ACCG  C         A A
                                   C C             C G
                                  C   C      C   CCACA
                                 T     C   ACC GGACA  C
                               C                 CCC
                           C  A
                          C C
                   CGCCCACA CCCC
                   CAACCGCCCC TACC
                         TA
                           C
TRANSFORMERS
		print STDERR $Transformers16;
	}
	elsif($TF==17)
	{
my $Transformers17 = <<TRANSFORMERS;
                                   AT         TG
                                   A0         0G
                                   1A0GT    GAT0       C
                                   GT0T0TT0CGCGCC     11T
                                C1TG GG00GA  1A GTGC  GTT AA0CA0
                                   GT0AG0G  CCAAAAGAAA0GA0GAAC
                                011AGGGATAAGGTA00GGCAACGAGGCT
                          10    T GCAAT111TTGACAGGGGAAAAAAG G
                          G1GG  CAAAAAGAAAAAGTGACGAAGA0AG0AG
                    T1CA0AAGT 0AAAAA00GGAGAAAAAT00  GAAAG0G T0
                      TAGGACTT0AAAAGG00AAAAGGAA100   0AT A   1
                        AA00AATAAA GG11GG011AGGGC0G   A
                         GAAGAAAAA A0CAGGGG1GGATGGCGCTA1
                         0GGA0A GCAGTAGAGC1AGAGAAGC    GAG
                         C0 GGAC     AA0AAGGGG0 AA     CGA
                         1    AG      AAAAGA  A GGC     GAAG
                                      GAAAA    AAG
                                        AGAG  AAA
                                         AGC  G0TC1
                                        GAGA  A0T    CCC
                                              0AA  CC1
                                            CCCCTC
                                      TTT CTTCC1
                                  CTTTTCCCCCCCC
                               CCCCCCCCCCCCCCCC
                         C1CCCCCCCCCCCCCCCCCC
                        11CCCCCCCCCCCCCCCCC11CC
                       11111111111111111111111111
                  1  11111   111111111111    1111 1
                 1111111111     1111        11111111
             1111111111                      11111111
                 1                           1
TRANSFORMERS
		print STDERR $Transformers17;
	}
	elsif($TF==18)
	{
my $Transformers18 = <<TRANSFORMERS;
                                             CA    GC1
                              TT1           CCGAA   GC1      011CAAG1
                        1   C110G11A       C11C1AA   111  CC1A GGA1CC1A
                 AG0111ATAA1A    ACGC      CCC11CA      A0AGG     C1AA
             1GA C11AAAAAG 11                CCCC C 11AGC01A1      AG
           G1111AGGGA0GCTCG1                 0  GGCA11ACCGC
           C1111011AAAAACC A                 0AAG0G1G0ACC1T1
   11C     1GGGCGT1AG01CG1                  1AT1G1011AGG111
     AG0 1TAAAGGG1C1CCCG00 AG             AC11   0C 1GCA1C1
   TT0 C11TAC    G1ACGAGA0A0               11T    ACC000A
      0T1GGA     AACT1AGAGAAC0              0CT  GAG1G A1
     1ACGA01    1GTTT   AAAAG AA                T1CA0  T1
      1CC        0C1AA1 GGAACC                  GGAAGC  A
                AGTAAA    AAA                   0GGGA GAAAGCA
                 GACAA                         GG 00  GCTAAGA
                   GA                          G ACC GAAT0AA
                   GT                                   G
                                               GCGGT 0GGGG0G1
    111  11111111111CCCCCC                   C1 11111111111
TRANSFORMERS
		print STDERR $Transformers18;
	}
	elsif($TF==19)
	{
my $Transformers19 = <<TRANSFORMERS;
                                     C
                              C    1AC
                             1A   ACAA1C
                            C11C00TGC TT
                          1CTCTGAA0C11CG0AGCAC
                    AA AAA AA1CC11CC 1T1GC AGG1
                   C11GAAAAA1CTACTAA111CCTA1C1
                   C    1AA     TGGGGGGG 1GAA0G
                     AA01111  C        1AACC1AAC1TC
                     TA1AAGATCTC    110AA01C 11ACT11
                    CA  AAC1CAAGAAAAA1TAG11G    AAAGTT
                   CGA   CCA     AA C1C0G0G1        G
                       G  G0GAA1C1T0 CA1AA0G
                       1GAACGG1C 1GCAGC1AAAGC1
                        GAAAGGCA00C0AAAA0 AC0G
                         AAGGACG10CCATGAGT1C1CTC1
                         AAAGTT         G0GTC 1AG0
                        GACAAAA11         A0A1 CA
                        GAATAAGC1           1TCCCA1
                        GAG0G0CTTC           GA0AA C0A
                        TAA0A0AGCG           CCGTAAAA
                        CAATA1C11C           GGAAAAC
                        CAA11A1CA0          1GGAAAC1
                        CCG0AA1TAC           AAAA11CG
                       1TG0CGTC0G            AA00GA1
                       AGTAC1ACCAT           1AAG11C
                        AC11111GA           1AA 1CT10
                        ACA1  1111
                     1TAA0AGG 1C1
TRANSFORMERS
		print STDERR $Transformers19;
	}
	elsif($TF==20)
	{
my $Transformers20 = <<TRANSFORMERS;
                       1AA
                      1TTAAA1                  ATTTT1
                      AAA1AT11A             A111 TAAA
                         AAAAT111    AA 1 AAAAAAAA
                             AAA   TT11A1AAAAAA
                    1111A1   A1A1 AAA AAAA1AAA
                 1 AAAA1A AA 1111 1AAAAAAAAAA1A1  1 1AAA
                1A1 1111A11AAAAAA1AAA1AA1AAAAAAA1  1A1 1A
                 A AA1AAAAAT A AA A11A1 1 AAA A1AA1AAAAA11
                  111AAAA1  A  1A 11A11AA1 AA 1 TAAAAAAAA11
                    1AAA1  1 1AT1A11AAAT1AA11  A1AAA T AAA
                       A1 AT   A1AAAAAAT1AA  A  A 1A
                          1T  1A T11A111AA AA
                                  1A11AAA AA11
                                AAAAA1AAAAAA A1
                              11AA11A  AAA1A11AA
                            A11AAA  T    AA1AT1T1
                           1A AAA1A      A1A ATAA
                           AAA A 1         1AAAA
                      1  11A1AAA1A1          1A1
                      AA11A1AAAAA1       1AA11AAA1
                      AAAAA1TA1          AAA111A1AA
                   1A  111 11A1             1A1AAA
                  T AA AA1 A1              A1TAA1A
                    A1AATA                AAA11A 1
                   AA111 1A             AAA1 AT111111
                      111AA             AA  A11A1AAA1
                                           1A111AAAA1
                                        AA1A1111A11AAA
                                       1ATA1111  1
TRANSFORMERS
		print STDERR $Transformers20;
	}
	elsif($TF==21)
	{
my $Transformers21 = <<TRANSFORMERS;
                                1 1
                                A11  1      C01
                                  1 A  C 1  11
                              A1A1AC 1AAA1  G1111
                             1A1C  A 1AA AA A 11
                            1AT T1A0AAAC 1AC111
                           TAAA1A AATA1GAA    11C1
                         111A01C0G CCTAACAC 1 A  01
                    1   0A1A1  A1A1ACA1A A11A10C11
                    1 TA 1A0   1 GAGAA11  11ACG1
                       GA1ACG A C T1A C1111A 1
                        C 1G1 AC1CGAA   1 C1
                       1A10A0 1T0GA1 AGC11
                              ATCAA 1AA111
                             110CGT   11
                             110C0A1  C  1
                             11G1G011G111AA0
                               11CG11     C
                              1AA11A   1   11
                               CACC1 1 AAG 0
                          1 111A1CC111A1AA 11
                        111GAT1GAGGA111C11 1A11111
                     1111111111CTAC1111111AAA
                          1A11C1111CA11 0  1 1
                         1 AA1CC   C111AGA  1 1
                                     1 AGAGC1A1
                                       TA  11G1
                                       CC  CCC
                                      C11  C1CC
                                            1
TRANSFORMERS
		print STDERR $Transformers21;
	}
	elsif($TF==22)
	{
my $Transformers22 = <<TRANSFORMERS;
                         C
                        1
                      11C11        1
                      1 01 1  CC   1 C
                      CT0C111  TA 011C
                      GCCTCCGGAAAAACTGATA        0
                        AA1C 1GAG11 AAAC0A  1   11
                       11A A00CCC  A 1AC1A  A   C1
                     1 TAA CAA0CGAA ACG0T0   GAGCC
                       CC0    ACA1TTAGCATA AATA
                    AGACA      C0AATCGTCAAA0AC
                       0       C1TCCA1AC
                   C0 TC       1GTCCC1AA
                    10AG       A1GACCGATC
                    1  T     C11110AA0CAC C
                    0 A11   111 11TAAG1GCGCC
                    T CGC   11  1TTGG0C1CGC1
                     1CG C  11  1CT 0CCAGT 1A1
                       TC   A0  T0A  TGA1A10CAC
                           AAA  TAC    GCTGCTCA1
                          1A TTGG A    G A A111A1
                         1G11A0ATG      0T0GCC1
                        1C11CA1TA0     GGAC011A
                        AGAAAAT1CG    100ACA1111
                        AAAAGAGA00   CGAAAA1 AAC
                       ACGAAATTA       AG001G0CAC
                       C1CCC011A        0GGC110A
                       1   GGCC1       TGT0CC1GA
                       110CGG1C1     GGGATGA1G1G11G1
                      1CT  GATG00   0A00AGAGACA11CCC
                      1ATC00GTCGA11CGGGTTTTG0ATA1CG111
                         1T1111
TRANSFORMERS
		print STDERR $Transformers22;
	}
	elsif($TF==23)
	{
my $Transformers23 = <<TRANSFORMERS;
                           1AAAAAAGGGGGGAAAAAAC
                     1AAGG00000TTTCCCTCCTT0TT0GGGAAA
        AG00GG0GAAAAGGGGG0T0CTCCCCCC1CCCCCTTTTT0GGGGGAAGAAGGGAGAGA
        A00000GG00GAGG0GGG0TTT00CCC11111CCTTTT00T00G0GA0GGGGGGGGG1
        1GT0TT0T0TCAA00T00CCGTTC111111111CCCGGTTT00G0AGG0000G00GA
         A000GT0T0T0G0GT0TCTCCCG         0GCCCCCTTTGGGTTTT00GT0GA
         A0T0CGGTT00GAG0G0TCC1111GG    GGCCCCCTCGGGAGG0T00GGC0GG1
          GTTC0GAG00GGGAA0CC11111111GG1111C11CT1TAGT0TTTGAAG0G0A
          ATTCGC1GA0000TCGA01C111111111111C111GG0CCCC0GAGTGGG0GA
          ATTC1GGTCGGA0TCCC00G11111111111C1CAGTCCCT0AG000GAC0GG1
           GTCCTGG0TTCGAT1CC11AGC1CC1CCC1GGGC1CCCTG0TCTGAAG00GA
           ACCCCCTGATC11CG01110C1TACT11AGC00111GGTCCT0AGG0G0C0A
           A0CCTTCC1TGTCC1T111CC1110AA0T1CGTCC1CCCTGAGTT000TTG1
            1GGTCC111CT0AC 111CA111111111CACC1C1TA0CCT1CCCC0G
             TGG0C1111C110A1111G111111111CGC110A0111C111CGGGG
             ACC0GG111CCC11T111TC1111111C0TC11GCC1C1C11GAGC1A
             AC11CCAGGTC11C1111CT1 111111GC11111110G0GGC111CA
             A1111TC        10CCA    11T AGGG1         11111A
             AC11111          GCA     111GGG         G  1111A
             ACC1111CCG   1GC111A 1    11GC11GA1   0111111111
             0111111111CA1111111A    11 CG111111TA1111  111T
             11111111111AC111111A  111 1GGC111110A1 11  111A
              GC11111111A1111111A1111111CG111111CG1111111CCA
              AC11111111AC111111A11111 11G1 1111CG111111CCCA
              AC1C111111A1111111AC111111CGC111110A 111111CTA
              ACC11CC111A1111111AC1111110GC1111C0G1111CCCCCA
              ACC11CCCCCAC111CCCAAAAAAAAAGC111C1TG1C1C1CTTTA
              AACCCCCCTCAC11CCCCCCCC111CCC1111110ACTTTCCCCA1
                AACC000CACC11CCC1CCCC111C11C11CCGATGG00TA1
                  AAGGG0A0CC111GAAAAAAAAAAGCCCCCGATGGGA
                    AGG0AGCCC1CAT11    1 G0TC11CGATTA1
                      A0A0CCC1GGC11111CCCCGGCCCCGAG1
                        AGTTCGAT1CCCCCTCCCGA0CTTAC
                          A0CGGTTTTTT0G00TCGGTGA
                            AGAAAAAAAAAAAAAAAG
TRANSFORMERS
		print STDERR $Transformers23;
	}
	elsif($TF==24)
	{
my $Transformers24 = <<TRANSFORMERS;
         T
          A1              C                    A              CA
          AAA             AA                  AG             AAA
          AGGA             GA               1GG             AAAA
          C000G            G0G1            AG0A           1GGGA
           G0G0GG          GCC0A          GTTTA          AGGTGG
           GT00TGG         G0111G       1GCCT0          G0TGCCA
           ATTTCTTGAAC      C11110GAAAGAGCCCCG      TAA0TTT0CCA
           1TTCCCT00GGTGAAG1T1111111111CCC1CCA 0AGTT0000TCC00G
            0CCCCTT00G000C1CA111111111111111CA0CCCCTTT00TT0GGG
            TCCCCCTC0T0CCC11A1111 C111C111110AC1CCCCTTC1C0GGGA
            ATCCCC1AAGT111110T1  11G  011111A011CCTTGAG1TGGGGA
            10TCCCCC1CCAAAATCG11 11CTA11111CAC0AAAA0CTTTGGGGG
             CTCCCCC111C111C1A1  1 1GG111C10A1C1C1CCCTT00GG0A
             G000CC 011111111G111  1C1111C1A0C111CC1CG100GGGA
             AGGTC110GAAGTC 1CG   1111 1111AC1 0GGAAGG1CCTTGG
             1CCC1111111CGAAA0AA        11AAAAGGTC1C11C11CTG
              A11C1CCCCCC111111CA1     1TAGC1111CC11CC111CCA
              CAC1CCC10 GC11  1 10C    0G11111111C C1C11CAAA
              ACACCCCCCG1TG11    1CG 1AT111111CAC10CCCCCGA1A
               C1A1CTTC00    G1 1  1AA1 1111A    0C1CCTGA1G
               AC1AGCTCCC1CT  T11      1111G  GTCCCCCTAA 1A
               A111AGCCC111  C10        11A1C111111TCAA111A
               G1CC AGCCC11111  T      1 01111111CCCAA 1110
                1C11 AGC1111111   1    11 1111C1111GG 11CG
                ACC11 AA1111111   11    1111111111GG 1CCCA
                ACCCCC AA111111111111111111111111GG 11CCCA
                TCCCCCC AG1C1CCC11C11111CC11C111G0 CCCTTC1
                 TTCCCCC1AGCT0TTCCCCCCCCCCC1111AG CCTT0TG
                 ATT0T00C1A0C0TTTT0CCC1C1CCCC1GG1TT0GGTCA
                 ATT0GGGTC1AGC0T0GG0TCCCCCCT1AG1TGGAAAA0A
                 GCTTGGG0CC1AGCCT0GGTCCCCTTCAA1CTTGGGGGG1
                  G000GGTCCC1A0CCTTCCCTTT01AA1CTTT00T00A
                   CGGGG00CTTCAACTT0CTCT01AAC0TTT0000G1
                      1A00TT0C1AACT00G00CAAC00T000G
                          AG0TTCAATGGGGTAA0GGGGA
                             TAAGAATGGGAAGAA1
                                 AAAT0AAA
                                    AA
TRANSFORMERS
		print STDERR $Transformers24;
	}
	elsif($TF==25)
	{
my $Transformers25 = <<TRANSFORMERS;
11111111111111111111111111111111111111111111111111111111111111111111111111
11111111111111111111111111111111111111111111111111111111111111111111111111
11111111111111111111110CC111111111111111111T 1C111111111111111111111111111
111111111111111111111  1C11111CT1 C1111111111 1111111111111111111111111111
11111111111111111111  1111111111C 11 111111     C1111111111111111111111111
11111111111111111T      0C11   1  C 1   1111 1   1111111111111111111111111
11111111111111111      CC11      11        C11     C1111111111111111111111
1111111111111111 C1  AT0G        11        CATC     1111111111111111111111
1111111111111C 1   GCGTA              C 1   GGAT01   111111111111111111111
111111111111C1  T0GAC0G                C     000CGG 1 C1111111111111111111
11111111111  0GGGGGTGATC1               1 111GGGG0GG0G1  C1111111111111111
111111111 GGGGG0GT00CGG11                1 1C00G0GG0G000G 1111111111111111
1111111C 1000000CCT11                        1C1GGGGT0000C  C1111111111111
11111111  GTTTTTT110111 1                  G 1CC01CCTTTTTC  C1111111111111
11111110   TTCCCG1101111G1          1    11111TCT 1GCCCTG  1C1111111111111
1111111C   TTCC11  TT1111CT              1111CT1C1 A1CCC   1C1111111111111
1111111101        1CT  A111     GCA01   0G11C  CC1        10C1111111111111
11111111C0        1 T111 GG 1 A0ACCAAG  C11 CG1010 1      CCC1111111111111
1111111CCCG   AC CC11GGC01G AG1T1  1AAA 1GT0GA0  C C1    CCCC1111111111111
11111111CCC01TC11G1110C1C1 0GAAG1111CGAAGA1AGA 1C1A1CC0CCCCC11111111111111
111111C1CC1C0CC1ACAA 1C0TTGTCCGC11C1GAACGCT0G 1 AAT1CCTCCCC111111111111111
111111111CCCCCTCACG CAACC1TG01 CG 1TC11AGG1A1 AC GGGCCC0CC1111111111111111
11111111111TCTTAG1 111G0TGA0T1 C11 111ATAAAGCC  1  GGCCCC11111111111111111
11111111111CCAAG  11  11GA01ACA11C1 CACA1GGGG   01  GATCC11111111111111111
1111111CTT111AG1  11A1TATCA01TAAGTCTA C1TA GACG011A  GA11CTT11111111111111
11111GCGTC0111G110CA 1AG 11 11  0GA0 111 C CAC  TCCGAA11 CT00CCC111111111C
11111CA00CCG 10GCT1  1ATT0T0C    11    0GAAAA1  11C0A11 TTTCGG1111111111CC
1111111AGAC0T 1ACG11111CCC A1 G  1  011CATCAC A 0TCG11 0TCGAG11111111C1111
1111111C00GAT0 1GT G 10CG0GCAGA0CGTCTAGA0C0G11 AA 0CCCTCAGGG11111111111CC
11111111GC1GGGG0CGCG0   CGAG 1G1CCC0GGC GA1T  10C0C A0AGT1GA111111111C11
11111111G111CGGGA000GG   1T1TACAAAAAGGA01C1  C0GGGCGAGG111CC11111111CCC
111111111111111A0AGGTGGGGT AA0G0GAAAGAGAGCAGGG0GAAGC0111111111111111CT1
1111111111111CC10GGGACGGC1CT01 TCC11G1CG1CT0GCC00GG111111111111T TGAC11  G
C1AA     1001GAAAAAAGT00T1C1C0C1CC0G1TC01C1CGCC0GAAG1TC     CA1C1 CA1   10
GA       CC10ACAAGAAG1GA1 11C0CGG  AG0C111110CGAGAAGC111      CG01C110 C0G
         1TC0A0GAGAAGAGT11 11GGCG0GA0CC111CTGTAAATGGC111         CTTAC1 CG
          T0CGAGGAAAAAAA0GC  1CC11 0 C1 1CGGAAAAAAG1TC1             0C0A
         1CCCGA1GAGGA0AAAAGC1 0CC 1AC1 CTGGAAAGATAGT1G11           1   GTA
 1      110CA0GGGCGAAGC0GA0AT  CA11A0 C0AGGACCGGGTCAT1C       11         C
      1CC011AAGCTGT0GG1AATTGG0 AGG0T CGTT0CA01GGT0GGG1AC
      11GC1AAAAGC0GGGGCGGGGAGAAGGGGCTAGGGGGGGAATCAAAA11C1
     1 TACGGGCAAG10AAGGGAAG0AGGGCCT0AGAT0GGGGGCCAA0GGA11TC
    11CT11AA1GGAATCTAAAGGC0CTA0G0TCG0TACTGAG0C1AAAG1AGC10T1
 C1  C0011A01AAAAAC10AGG11CG1CGGAGGAGA1TGGA011AAAAGG1G1C1GT1 G1
  CCCCC0C11GTGAAAAA CCAAGGGGGGGAAAGGAGAGGG0C AAAAAA0GCCCTTCCCC
C011CCCCGGC10AAA CA C1CGGAAGTCGGAAAGAACCAA00 AA AGAA1CCGCC1CC1TC1G1
TRANSFORMERS
		print STDERR $Transformers25;
	}
	elsif($TF==26)
	{
my $Transformers26 = <<TRANSFORMERS;
                               1 C   A
                            TC01  A11111C  11CA
                            CAGC1  CA11 AGAAGAC
                              AAA A  AG AAT 11
                           1CGTA CC GC11GAGC10
                      1     G     11AA1  GA    11A1
                            A0A01 TCGA1 0T1111  1
                      1C0CAGTA   A CT1A 1A0AAATT 1A
                    G    A0GA0 G110CGC11CTCCA 0CGA    C 111
                 1   1 1A1CCA 1CC1AAGAT1GA0A TCT1     0CGT
                   A  1CCT  TCG0A1CCTACAAAA11AG10C  G 1TAC
                   AC0CG AG ATACAAATAAA0 11 1A  CG TG0 G1
                   1AC1TG     TGAAA01AGGGAT1 G    CAA1C
                  11G10CA    1TG T1GCAAG1A11       GT11
                   1AAAAA     0C CGAA A 111    C GG
                   1CA0ACG    1G  AG1AG  GC1   T 1GA
                   10GAGGGC  TTTC1GC  G AAG0CA AC0
                      CGA G     TCTAG A1G0AC
                      AA1A  1 CAGTC    AG C
                      CC0T 1C CAAA1G  GA1GAT  1
                           110TCTCC   GATTAA1 1
                           1TA0CGG     AA1CT0G
                           CT11AA1      1GCT0C
                          11CC G1C1     C GCC C
                            ATC1A1C   C A1TCGA 1
                            1 0 0T0   AAGG1A A
                             CGG TC    AACCAC1 T 1
                        AA1C1G1   A   CAA G1 1CAC0
                       AGA10ACGTCC111111C0 T11ACT
                          CC1 1A111     11AA11C 11
                                          C111
TRANSFORMERS
		print STDERR $Transformers26;
	}
}