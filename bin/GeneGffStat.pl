#!/usr/bin/perl -w
=head1 Name

GeneGffStat.pl -- Gene Finding Statistical Program based on GFF3 file format

=head1 Version

Author: BENM <binxiaofeng@gmail.com>

Version: v0.4 Aurora, Dec-12th-2011

=head1 Usage

  --gff <infile>		input annotation GFF3 file
  --cmpref <infile>		input compared reference GFF3 file, need disposeOverlap.pl
  --cmpid <str>			set compared region ID, default: CDS
  --genome <infile>		input genome fasta
  --size <int>			input genome size
  --ignore <str>		ingore 'Exon' or 'CDS', if both exist, default: Exon:gene

=head1 Example

perl GeneGffStat.pl -gff W14_V7.2.gth.gff3 -genome W14_V7.2.fasta -ignore BENM > V7.2_GeneStat.txt 

=cut
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my %opts;
my ($GFF,$CmpGFF,$CmpID,$Genome,$genomeSize,$Ignore,$Verbose,$Help);
GetOptions(
	\%opts,
	"gff:s"=>\$GFF,
	"cmpref:s"=>\$CmpGFF,
	"cmpid:s"=>\$CmpID,
	"genome:s"=>\$Genome,
	"size:i"=>\$genomeSize,
	"ignore:s"=>\$Ignore,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if ((!defined $GFF)||($Help));
$CmpID ||= "CDS";
$Ignore ||= "Exon:gene";
my @IgnoreID = split /\:/,$Ignore;

#my %Abbrev = (
#	#'A' => [ 'A' ],
#	#'C' => [ 'C' ],
#	#'G' => [ 'G' ],
#	#'T' => [ 'T' ],
#	'M' => [ 'A', 'C' ],
#	'R' => [ 'A', 'G' ],
#	'W' => [ 'A', 'T' ],
#	'S' => [ 'C', 'G' ],
#	'Y' => [ 'C', 'T' ],
#	'K' => [ 'G', 'T' ],
#	'V' => [ 'A', 'C', 'G' ],
#	'H' => [ 'A', 'C', 'T' ],
#	'D' => [ 'A', 'G', 'T' ],
#	'B' => [ 'C', 'G', 'T' ],
#	'X' => [ 'A', 'C', 'G', 'T' ],
#	'N' => [ 'A', 'C', 'G', 'T' ]
#);
#
#my %Heterozygous;
#while ( my ($key, $value) = each(%Abbrev) )
#{
#	my $het1=join ",",@$value;
#	my $het2=join "/",@$value;
#	$Heterozygous{$het1}=$key;
#	$Heterozygous{$het2}=$key;
#}

my %fasta = ();
my %gene = ();
my %ref = ();
my %intergenic = ();
my %gene_strand = ();
my %sstype = ();

######### need information ############
$genomeSize ||= 0;
my $countGC = 0;
my $genomeGC = 0;
my $geneNumber = 0;
my $geneDensity = 0;
my $codingRegionLength = 0;

my $geneLength = 0;
my $exonNumber = 0;
my $exonLength = 0;
my $exonMeanLen = 0;
my $exonMinLen = 0;
my $exonMaxLen = 0;
my $meanExonPerGene = 0;
my $exonGC = 0;

my $intronNumber = 0;
my $intronLength = 0;
my $intronMeanLen = 0;
my $intronMinLen = 0;
my $intronMaxLen = 0;
my $meanIntronPerGene = 0;
my $intronGC = 0;
my $averLenOfIntergenic = 0;
my $minLenOfIntergenic = 0;
my $maxLenOfIntergenic = 0;

########################################

if ((defined $Genome)&&(-f $Genome)){
	read_geneGff($GFF,\%gene) if ((defined $GFF)&&(-f $GFF));
	read_fasta($Genome,\%fasta);
	statisticGene();
print "Annotated genome size:\t\t$genomeSize\n";
print "GC(\%):\t\t\t\t$genomeGC\n";
print "Gene Number:\t\t\t$geneNumber\n";
print "Gene Length:\t\t\t$geneLength\n";
print "Coding Region Length:\t\t$codingRegionLength\n";
print "Gene Density(\%):\t\t$geneDensity\n";
print "Mean Length of Intergenic:\t$averLenOfIntergenic\n";
print "Minimum Length of Intergenic:\t$minLenOfIntergenic\n";
print "Maxmium Length of Intergenic:\t$maxLenOfIntergenic\n";

print "Exon Number:\t\t\t$exonNumber\n";
print "Exon Number/Gene:\t\t$meanExonPerGene\n";
print "Exon Length:\t\t\t$exonLength\n";
print "Mean Length of Exon:\t\t$exonMeanLen\n";
print "Minimum Length of Exon:\t\t$exonMinLen\n";
print "Maxmium Length of Exon:\t\t$exonMaxLen\n";
print "GC(\%) of Exon:\t\t\t$exonGC\n";

print "Intron Number:\t\t\t$intronNumber\n";
print "Intron number/Gene:\t\t$meanIntronPerGene\n";
print "Intron Length:\t\t\t$intronLength\n";
print "Mean Length of Intron:\t\t$intronMeanLen\n";
print "Minimum Length of Intron:\t$intronMinLen\n";
print "Maxmium Length of Intron:\t$intronMaxLen\n";
print "GC(\%) of Intron:\t\t$intronGC\n";

print "\nIntron Splice site subtypes Statistics:\n";
print "type\tfrequency\n";
foreach my $k(sort{$sstype{$b}<=>$sstype{$a}}keys %sstype){
	print "$k\t$sstype{$k}\n";
}
##Splice site subtypes
## GT--AG
## GC--AG
## AT--AC
}

compare_gff($GFF,$CmpGFF) if ((defined $CmpGFF)&&(-f $CmpGFF)&&($genomeSize>0));

my ($TP,$TN,$FP,$FN)=(0,0,0,0);
#my $sn=cal_Sn($TP,$TN,$FP,$FN);
#my $sp=cal_Sp($TP,$TN,$FP,$FN);
#my $cc=cal_CC($TP,$TN,$FP,$FN);
#my $acp=cal_ACP($TP,$TN,$FP,$FN);
#my $ac=cal_AC($acp);

sub read_geneGff{
	my ($file,$r_hash) = @_;
	#warn "reading $file...\n";
	open IN,"$file" or die "$!\n";
	while(my $line = <IN>){
		chomp $line;
		my @info = split /\t+/,$line;
		my $i=0;
		next if (($line eq "")||(@info<8));
		foreach my $ig(@IgnoreID){if($info[2] =~ /$ig/i){$i++}}
		next if ($i>=1);
		if(($info[2] =~ /mRNA/i)||($info[2] =~ /gene/i)){
			push @{$intergenic{$info[0]}},[$info[3],$info[4]];
			$geneLength += $info[4] - $info[3] + 1;
			next;
		}
		my ($geneID) = $1 if ($info[-1] =~ /Parent\=([^\;]+)/);
		push @{$$r_hash{$info[0]}{$geneID}},[$info[3],$info[4]];
		$gene_strand{$geneID} = $info[6];
	}
	close IN;
	#warn "\t\tdone!\n";
}

sub read_fasta{
	my ($file,$r_hash) = @_;
	#warn "reading $file ...\n";
	open (IN,$file) or die "Error in reading file [$file]:$!\n";
	#my ($Chr,$Seq)=("","");
	my $Chr="";
	while(<IN>)
	{
		if (/^\>(\S+)/)
		{
			$Chr = $1;
		}
		else
		{
			s/\s*//g;
			$$r_hash{$Chr} .= uc($_);
			$genomeSize += ($_=~s/([ACTG])/$1/ig);
			$countGC += ($_=~s/[GC]//ig);
		}
	}
	close IN;
	#warn "\t\tdone!\n";
}

sub compare_gff{
	my ($gff,$ref)=@_;
#ID      FirstTableID    Start   End     Length  OverlapNumber   OverlapSize     OverlapRate     SecondTableID:Start,End,Length,OverlapSize...
	my $overlap1=qq(disposeOverlap.pl --E O --i1 $ref --f1 0,6-2-3-4 --i2 $gff --f2 0,6-2-3-4 --M $CmpID | awk 'BEGIN{tp=0;fn=0}{if(\$1!~/ID/){if (\$5>\$7){tp+=\$7;fn+=\$5-\$7}else{tp+=\$5}}}END{print tp,fn}'); #TP,FN
	open (OL1,"$overlap1|");
	my $stat1=<OL1>;
	chomp $stat1;
	close OL1;
	my $overlap2=qq(disposeOverlap.pl --E O --i1 $gff --f1 0,6-2-3-4 --i2 $ref --f2 0,6-2-3-4 --M $CmpID | awk 'BEGIN{fp=0}{if(\$1!~/ID/){if (\$5>\$7){fp+=\$5-\$7}}}END{print fp}'); #FP
	open (OL2,"$overlap2|");
	my $stat2=<OL2>;
	chomp $stat2;
	close OL2;
	my ($TP,$FN)=split /\s+/,$stat1;
	my  $FP=$stat2;
	my $TN=$genomeSize-$TP-$FN-$FP;
	my $SN=cal_Sn($TP,$TN,$FP,$FN);
	my $SP=cal_Sp($TP,$TN,$FP,$FN);
	my $CC=cal_CC($TP,$TN,$FP,$FN);
	my $ACP=cal_ACP($TP,$TN,$FP,$FN);
	my $AC=cal_AC($ACP);
	print "\n$gff compare to $ref\n";
	print "True Positive(TP):\t\t$TP\n";
	print "True Negative(TN):\t\t$TN\n";
	print "False Positive(FP):\t\t$FP\n";
	print "False Negative(FN):\t\t$FN\n";
	print "Sensitivity(Sn):\t\t$SN\n";
	print "Specificity(Sp):\t\t$SP\n";
	print "Correlation Coefficient(CC):\t\t$CC\n";
	print "Approximate Correlation Probability(ACP):\t$ACP\n";
	print "Approximate Correlation(AC):\t\t$ACP\n";
}

#Sn=TP/(TP+FN)
sub cal_Sn{
	my ($TP,$TN,$FP,$FN)=@_;
	return($TP/($TP+$FN));
}

#Sp=TP/(TP+FP)
sub cal_Sp{
	my ($TP,$TN,$FP,$FN)=@_;
	return($TP/($TP+$FP));
}

#CC=((TP*TN)-(FN*FP))/sqrt((TP+FN)*(TN+FP)*(TP+FP)*(TN+FN))
sub cal_CC{
	my ($TP,$TN,$FP,$FN)=@_;
	return (($TP*$TN)-($FN*$FP))/sqrt(($TP+$FN)*($TN+$FP)*($TP+$FP)*($TN+$FN));
}

#ACP=1/4*(TP/(TP+FN)+TP/(TP+FP)+TN/(TN+FP)+TN/(TN+FN))
sub cal_ACP{
	my ($TP,$TN,$FP,$FN)=@_;
	return (1/4*($TP/($TP+$FN)+$TP/($TP+$FP)+$TN/($TN+$FP)+$TN/($TN+$FN)));
}

#AC=(ACP-0.5)*2
sub cal_AC{
	my $acp=shift;
	return (($acp-0.5)*2);
}

sub statisticGene{
	#warn "start statistic ... \n";
	my @geneLen = ();
	my @exonLen = ();
	my @intronLen = ();
	my $exonSeq = '';
	my $intronSeq = '';
	#my @intergenic = ();
	my $lengthOfIntergenic = 0;
	my $intergenicCount = 0;
	foreach my $contig(keys %gene){
		my @codingList = sort{$a->[0] <=> $b->[0]}@{$intergenic{$contig}};
		if(@codingList>1){
			for(my $i=0;$i<@codingList-1;$i++){
				my ($s1,$e1) = @{$codingList[$i]};
				my ($s2,$e2) = @{$codingList[$i+1]};
				my $intergenic_len = ($s2>$e1) ? $s2 - $e1 -1 : 0;
				$lengthOfIntergenic += $intergenic_len;
				$minLenOfIntergenic = $intergenic_len if (($intergenic_len<$minLenOfIntergenic)||($minLenOfIntergenic==0));
				$maxLenOfIntergenic = $intergenic_len if ($intergenic_len>$maxLenOfIntergenic);
				$intergenicCount++;
				#$geneLength += $e1 - $s1 + 1;
			}
			#$geneLength += $codingList[-1][1] - $codingList[-1][0] + 1;
		}
		if (!exists $fasta{$contig})
		{
			warn "$contig doesn't exist!\n";
		}
		foreach my $geneID(keys %{$gene{$contig}}){
			$geneNumber++;
			my @cdsStartEnd = ($gene_strand{$geneID} eq "+") ? sort{$a->[0]<=>$b->[0]}@{$gene{$contig}{$geneID}} : sort{$b->[0]<=>$a->[0]}@{$gene{$contig}{$geneID}};
			for (my $i=0;$i<@cdsStartEnd;$i++){
				my ($s,$e) = @{$cdsStartEnd[$i]};
				$exonNumber++;
				my $exonlen = $e-$s+1;
				$codingRegionLength += $exonlen;
				$exonMinLen = $exonlen if (($exonlen<$exonMinLen && $exonlen>=3)||($exonMinLen==0));
				$exonMaxLen = $exonlen if ($exonlen>$exonMaxLen);
				my $sublen=(length($fasta{$contig})-$s>$exonlen)?$exonlen:length($fasta{$contig})-$s;
				my $subseq = substr($fasta{$contig},$s-1,$sublen);
				if ($gene_strand{$geneID} eq "-"){
					$subseq = reverse($subseq);
					$subseq =~ tr/ACGTacgt/TGCAtgca/;
				}
				$exonSeq .= $subseq;
				if ($i>0)
				{
					my ($s1,$e1) = @{$cdsStartEnd[$i-1]};
					my ($s2,$e2) = @{$cdsStartEnd[$i]};
					if ($gene_strand{$geneID} eq "-")
					{
						my ($s,$e)=($s1,$e1);
						($s1,$e1)=($s2,$e2);
						($s2,$e2)=($s,$e);
					}
					$intronNumber++;
					my $intronlen = ($s2>$e1) ? ($s2-$e1-1) : 0;
					$intronLength += $intronlen;
					$intronMinLen =  $intronlen if (($intronlen<$intronMinLen)||($intronMinLen==0));
					$intronMaxLen =  $intronlen if ($intronlen>$intronMaxLen);
					$sublen=(length($fasta{$contig})-$e1>$intronlen)?$exonlen:length($fasta{$contig})-$e1;
					my $cutSeq = substr($fasta{$contig},$e1,$intronlen);
					if ($gene_strand{$geneID} eq "-"){
						$cutSeq = reverse($cutSeq);
						$cutSeq =~ tr/ACGTacgt/TGCAtgca/;
					}
					$intronSeq .= $cutSeq;
					$sstype{"$1--$2"}++ if ($intronlen > 4 && $cutSeq=~/^(\w{2}).*(\w{2})$/);
				}
			}
		}
	}
	$averLenOfIntergenic = sprintf("%.2f",$lengthOfIntergenic/$intergenicCount);

	$exonLength = $codingRegionLength;
	$exonMeanLen = sprintf("%.2f",$exonLength/$exonNumber);
	$meanExonPerGene = sprintf("%.2f",$exonNumber / $geneNumber);

	$intronMeanLen = sprintf("%.2f",$intronLength/$intronNumber);
	$meanIntronPerGene = sprintf("%.2f",$intronNumber / $geneNumber);
	
	$geneDensity = sprintf("%.2f",$exonLength/$genomeSize*100);
	$exonGC = calGC(\$exonSeq);
	$intronGC = calGC(\$intronSeq);
	$genomeGC = sprintf("%.2f",$countGC/$genomeSize * 100);
	#warn "\t\tdone!\n";
}

sub calGC{
	my ($seq) = @_;
	my $tmp = length $$seq;
	my $aCount = ($$seq =~ s/[Aa]/A/g);
	my $tCount = ($$seq =~ s/[Tt]/T/g);
	my $gCount = ($$seq =~ s/[Gg]/G/g);
	my $cCount = ($$seq =~ s/[Cc]/C/g);
	my $gcContent = sprintf("%.2f",($gCount+$cCount)/($gCount+$cCount+$aCount+$tCount) * 100);
	return $gcContent;
}

#W14_V6.2.scf00008_GC=31.79_Length=62092 EVM     gene    6682    7633    .       -       .       ID=evm.TU.W14_V6.2.scf00008_GC=31.79_Length=62092.1; Name=EVM%20prediction%20W14_V6.2.scf00008_GC%3D31.79_Length%3D62092.1; 
#W14_V6.2.scf00008_GC=31.79_Length=62092 EVM     mRNA    6682    7633    .       -       .       ID=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1; Parent=evm.TU.W14_V6.2.scf00008_GC=31.79_Length=62092.1; 
#W14_V6.2.scf00008_GC=31.79_Length=62092 EVM     exon    7541    7633    .       -       .       ID=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1.exon1; Parent=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1
#W14_V6.2.scf00008_GC=31.79_Length=62092 EVM     CDS     7541    7633    .       -       0       ID=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1.cds1; Parent=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1
#W14_V6.2.scf00008_GC=31.79_Length=62092 EVM     exon    7383    7424    .       -       .       ID=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1.exon2; Parent=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1
#W14_V6.2.scf00008_GC=31.79_Length=62092 EVM     CDS     7383    7424    .       -       0       ID=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1.cds2; Parent=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1
#W14_V6.2.scf00008_GC=31.79_Length=62092 EVM     exon    7212    7295    .       -       .       ID=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1.exon3; Parent=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1
#W14_V6.2.scf00008_GC=31.79_Length=62092 EVM     CDS     7212    7295    .       -       0       ID=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1.cds3; Parent=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1
#W14_V6.2.scf00008_GC=31.79_Length=62092 EVM     exon    6777    6931    .       -       .       ID=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1.exon4; Parent=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1
#W14_V6.2.scf00008_GC=31.79_Length=62092 EVM     CDS     6777    6931    .       -       0       ID=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1.cds4; Parent=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1
#W14_V6.2.scf00008_GC=31.79_Length=62092 EVM     exon    6682    6697    .       -       .       ID=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1.exon5; Parent=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1
#W14_V6.2.scf00008_GC=31.79_Length=62092 EVM     CDS     6682    6697    .       -       2       ID=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1.cds5; Parent=evm.model.W14_V6.2.scf00008_GC=31.79_Length=62092.1