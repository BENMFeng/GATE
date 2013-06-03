#!/usr/bin/perl -w
#for Cassava Project

=head1 Name

Copyright (c) 2008-2012 BENM(Binxiao) Feng                            
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

=head1 Description

=head1 Version
Author: BENM <binxiaofeng@gmail.com>
Date: Dec 20, 2011
Update: Aug 29, 2012
Version: 0.5a, Luna

=head1 Usage

	--dir <directory>				input output directory, default: ./addinfo
	--gff <infile>					input gene prediction gff3 file
	--add <infile[:int,int-int-int-int]>		input addictional gene structure file:chr,strand-id-sstart-end
	--repeat <infile[:int-int-int-int]>		input repeat info file:chr-id-sstart-end
	--support <infile[:int,int-int-int-int]>	input cDNA/EST or protein info file:chr,strand-id-sstart-end
	--annotation <infile>				input anntattion file, such as GO info
	--abi						abinitio add structure (start codon, stop codon, 5'-UTR, 3'-UTR),
							need reference sequence input
	--ref <infile>					input reference seuqnece
	--verbose					ouput running process and error info to error screen
	--help						help

=head1 Exmple

$ perl GeneGffMasker.pl -dir ./overlap -gff gene_prediction.gff -add \
  structure.txt:0,1-2-3-4 -repeat repeat.txt:0-1-2-3 -support \
  support.txt:0-1-2-3

=cut

use strict;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Getopt::Long;

#Chr:strand	ID	Start	End
#GFF
#Here is a brief description of the GFF fields:
#1.seqname - The name of the sequence. Must be a chromosome or scaffold.
#2.source - The program that generated this feature.
#3.feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
#4.start - The starting position of the feature in the sequence. The first base is numbered 1.
#5.end - The ending position of the feature (inclusive).
#6.score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
#7.strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
#8.frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
#9.group - All lines with the same group are linked together into a single item.

my ($Dir,$Gff,$addStructure,$Repeat,$Support,$Annotation,$AbInitio,$RefFa,$Verbose,$Help);
my %opts;
GetOptions (\%opts,
	"dir:s"=>\$Dir,
	"gff:s"=>\$Gff,
	"add:s"=>\$addStructure,
	"repeat:s"=>\$Repeat,
	"support:s"=>\$Support,
	"annotation:s"=>\$Annotation,
	"abi"=>\$AbInitio,
	"ref:s"=>\$RefFa,
	"verbose!"=>\$Verbose,
	"help!"=>\$Help
);
die `pod2text $0` if ( (defined $Help) || (!defined $Gff) );

if (!defined $Dir)
{
	`mkdir -p ./addinfo`;
	$Dir="./addinfo";
}

my $start_codon="ATG";
my @stop_codon=("TAG","TAA","TGA");
my $Gfffile=(split /\//,$Gff)[-1];
my $region_column="0-1-2-3";
if (defined $AbInitio)
{
	if (!defined $RefFa)
	{
		print STDERR "no reference sequence input\n";
		die `pod2text $0`;
	}
	my %ORFsite=predictORF($Gff,$RefFa);
	open (IN,$Gff) || die $!;
	open (OUT,">$Dir/$Gff.abinitio.addutr.gff") || die $!;
	my %infocount=();
	my %Info=();
	my $id=0;
	my $geneid="";
	while(<IN>)
	{
		s/\s*$//;
		my @t=split /\t+/,$_;
		$t[1]="BENM";
		if (@t>8) ## && $t[8]=~/ID\=\S+($t[0]\.\d+)/
		{
			if ($t[2] =~ /gene/i || $t[2] =~ /mRNA/i)
			{
				if ($t[2]=~/gene/i)
				{
					delete @infocount {keys %Info};
					for my $k(keys %Info)
					{
						$Info{$k}=0;
						delete $Info{$k};
					}
					%Info=();
					$id++;
					if ($id>1)
					{
						print OUT "\n";
					}
				}
				my $b=($t[2]=~/gene/i)?"g":"m";
				if ($t[8]=~/ID\=\S+\.([^\s\.\;]+\.\d+)/)
				{
					$geneid=$1;
				}
				else
				{
					$geneid="$t[0]\.$id";
				}
				$t[8]="ID\=fbx.$geneid.$b; Name\=fbx.$geneid\;";
			}
			else
			{
				$Info{$t[2]}++;
				$t[8]="ID\=fbx.$geneid.$t[2]$Info{$t[2]}; Parent\=fbx.$geneid\;";
			}
			if (exists $ORFsite{$id})
			{
				if ($t[2] =~ /gene/i || $t[2] =~ /mRNA/i)
				{
					delete @infocount {keys %infocount};
					for my $k(keys %infocount)
					{
						$infocount{$k}=0;
						delete $infocount{$k};
					}
					%infocount=();
					print OUT ((join "\t",@t),"\n");
				}
				elsif ($t[2] =~ /exon/i || $t[2]=~ /UTR/i)
				{
					if ($ORFsite{$id}{start} < $t[3])
					{
						if ($t[6] eq "+")
						{
							if ($ORFsite{$id}{end}>=$t[3]+2)
							{
								my @tmp=@t;
								my $mark=lc($t[2]);
								$tmp[2]="CDS";
								$tmp[4]=($ORFsite{$id}{end}>$t[4])?$t[4]:$ORFsite{$id}{end};
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
							if ($ORFsite{$id}{end}>=$t[3]+2 && $ORFsite{$id}{end}<=$t[4])
							{
								my @tmp=@t;
								my $mark=lc($tmp[2]);
								$tmp[2]="stop_codon";
								$tmp[3]=$ORFsite{$id}{end}-2;
								$tmp[4]=$ORFsite{$id}{end};
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
							if ($ORFsite{$id}{end}>=$t[3] && $ORFsite{$id}{end}<$t[4])
							{
								my @tmp=@t;
								my $mark=lc($tmp[2]);
								$tmp[2]="3-UTR";
								$tmp[3]=$ORFsite{$id}{end}+1;
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
						}
						else
						{
							my @tmp=@t;
							my $mark=lc($tmp[2]);
							$tmp[2]="5-UTR";
							$infocount{$tmp[2]}++;
							my $subinfo="$tmp[2]$infocount{$tmp[2]}";
							$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
							print OUT ((join "\t",@tmp),"\n");
						}
					}
					elsif ($ORFsite{$id}{start}>=$t[3]&&$ORFsite{$id}{start}<=$t[4])
					{
						if ($t[6] eq "+")
						{
							if ($ORFsite{$id}{start}>$t[3])
							{
								my @tmp=@t;
								my $mark=lc($tmp[2]);
								$tmp[2]="5-UTR";
								$tmp[4]=$ORFsite{$id}{start}-1;
								$infocount{$tmp[2]}++;
								my $subinfo="5-UTR$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
							if ($ORFsite{$id}{start}<=$t[4]-2)
							{
								my @tmp=@t;
								my $mark=lc($tmp[2]);
								$tmp[2]="start_codon";
								$tmp[3]=$ORFsite{$id}{start};
								$tmp[4]=$ORFsite{$id}{start}+2;
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
							if ($ORFsite{$id}{end}>$t[3]+2)
							{
								my @tmp=@t;
								$tmp[2]="CDS";
								$tmp[3]=$ORFsite{$id}{start};
								$tmp[4]=($ORFsite{$id}{end}>$t[4])?$t[4]:$ORFsite{$id}{end};
								my $mark=lc($t[2]);
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
							if ($ORFsite{$id}{end}<=$t[4])
							{
								my @tmp=@t;
								my $mark=lc($tmp[2]);
								$tmp[2]="stop_codon";
								$tmp[3]=$ORFsite{$id}{end}-2;
								$tmp[4]=$ORFsite{$id}{end};
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
							if ($ORFsite{$id}{end}<$t[4])
							{
								my @tmp=@t;
								my $mark=lc($tmp[2]);
								$tmp[2]="3-UTR";
								$tmp[4]=$ORFsite{$id}{end}+1;
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
						}
						else
						{
							if ($ORFsite{$id}{start}<$t[4])
							{
								my @tmp=@t;
								my $mark=lc($tmp[2]);
								$tmp[2]="5-UTR";
								$tmp[3]=$ORFsite{$id}{start}+1;
								$infocount{$tmp[2]}++;
								my $subinfo="5-UTR$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
							if ($ORFsite{$id}{start}>=$t[3]+2)
							{
								my @tmp=@t;
								my $mark=lc($tmp[2]);
								$tmp[2]="start_codon";
								$tmp[3]=$ORFsite{$id}{start}-2;
								$tmp[4]=$ORFsite{$id}{start};
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
							if ($ORFsite{$id}{end}<$t[4]-2 && $ORFsite{$id}{start}>=$t[3]+2)
							{
								my @tmp=@t;
								$tmp[2]="CDS";
								$tmp[3]=($ORFsite{$id}{end}<$t[3])?$t[3]:$ORFsite{$id}{end};
								$tmp[4]=$ORFsite{$id}{start};
								my $mark=lc($t[2]);
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
							if ($ORFsite{$id}{end}>=$t[3])
							{
								my @tmp=@t;
								my $mark=lc($tmp[2]);
								$tmp[2]="stop_codon";
								$tmp[3]=$ORFsite{$id}{end};
								$tmp[4]=$ORFsite{$id}{end}+2;
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
							if ($ORFsite{$id}{end}>$t[3])
							{
								my @tmp=@t;
								my $mark=lc($tmp[2]);
								$tmp[2]="3-UTR";
								$tmp[4]=$ORFsite{$id}{end}-1;
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
						}
					}
					elsif($ORFsite{$id}{start}>$t[4])
					{
						if ($t[6] eq "+")
						{
							my @tmp=@t;
							my $mark=lc($tmp[2]);
							$tmp[2]="5-UTR";
							$infocount{$tmp[2]}++;
							my $subinfo="$tmp[2]$infocount{$tmp[2]}";
							$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
							print OUT ((join "\t",@tmp),"\n");
						}
						else
						{
							if ($ORFsite{$id}{end}<$t[4]-2)
							{
								my @tmp=@t;
								$tmp[2]="CDS";
								my $mark=lc($t[2]);
								$tmp[3]=($ORFsite{$id}{end}<$t[3])?$t[3]:$ORFsite{$id}{end};
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
							if ($ORFsite{$id}{end}>=$t[3] && ($ORFsite{$id}{end}<=$t[4]-2))
							{
								my @tmp=@t;
								my $mark=lc($tmp[2]);
								$tmp[2]="stop_codon";
								$tmp[3]=$ORFsite{$id}{end};
								$tmp[4]=$ORFsite{$id}{end}+2;
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
							if ($ORFsite{$id}{end}>$t[3])
							{
								my @tmp=@t;
								my $mark=lc($tmp[2]);
								$tmp[2]="3-UTR";
								$tmp[3]=$t[3];
								$tmp[4]=($ORFsite{$id}{end}-1>$t[4])?$t[4]:$ORFsite{$id}{end}-1;
								$infocount{$tmp[2]}++;
								my $subinfo="$tmp[2]$infocount{$tmp[2]}";
								$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
								print OUT ((join "\t",@tmp),"\n");
							}
						}
					}
				}
			}
			else
			{
				print OUT ((join "\t",@t),"\n") if ( $t[2] ne "CDS");
			}
		}
		#else
		#{
		#	print OUT $_;
		#}
	}
	close IN;
	close OUT;
	&check_gff("$Dir/$Gff.abinitio.addutr.gff","$Dir/$Gff.abinitio.addutr.gff.check");
}
elsif (defined $addStructure)
{
	$region_column=get_region($Gff,"$Dir/$Gfffile.region","1,7,3,4,5");
	print STDERR "get_region\n" if (defined $Verbose);
	print STDERR "addStructure\n" if (defined $Verbose);
	my ($structure_table,$structure_column)=("","");
	&parse_table($addStructure,\$structure_table,\$structure_column);
	&overlap("$Dir/$Gfffile.region",$region_column,$structure_table,$structure_column,0,"$Dir/$Gfffile.region_overlap_addStructure.txt");
	my %sturcture_hash=();
	&get_overlap_info("$Dir/$Gfffile.region_overlap_addStructure.txt",\%sturcture_hash);
	&add_structure($Gff,"$Dir/$Gfffile\.addStructure",\%sturcture_hash);
	for (keys %sturcture_hash)
	{
		delete $sturcture_hash{$_};
	}
	&check_structure("$Dir/$Gfffile\.addStructure","$Dir/$Gfffile\.addStructure.check","1,3,4,5");
	&identify_utr("$Dir/$Gfffile\.addStructure.check","$Dir/$Gfffile\.addStructure.check.addutr");
	&check_gff("$Dir/$Gfffile\.addStructure.check.addutr","$Dir/$Gfffile\.addStructure.check.addutr.check");
}

my $infile=(-f "$Dir/$Gfffile\.addStructure.check.addutr.check")?"$Dir/$Gfffile\.addStructure.check.addutr.check":"$Dir/$Gfffile";
if (defined $Repeat)
{
	print STDERR "add repeat info\n" if (defined $Verbose);
	$region_column=get_region($infile,"$infile.region","1,3,4,5") unless (-f "$infile.region");
	my ($repeat_table,$repeat_column)=("","");
	&parse_table($Repeat,\$repeat_table,\$repeat_column);
	&overlap("$infile.region","$region_column",$repeat_table,$repeat_column,"0.1-small","$infile.overlaprepeat.txt");
	my %repeat_hash;
	&get_overlap_info("$infile.overlaprepeat.txt",\%repeat_hash,"gene");
	&add_support($infile,"$infile.addrepeat",\%repeat_hash,"Repeat");
	&filter_repeat("$infile.addrepeat","$infile.filterrepeat");
}

$infile=(-f "$infile.filterrepeat")?"$infile.filterrepeat":$infile;
if (defined $Support)
{
	print STDERR "add support" if (defined $Verbose);
	$region_column=get_region($infile,"$infile.region","1,7,3,4,5") unless (-f "$infile.region");
	my ($support_table,$support_column)=("","");
	&parse_table($Support,\$support_table,\$support_column);
	&overlap("$infile.region",$region_column,$support_table,$support_column,0,"$infile.region.overlapsupport.txt");
	my %support_hash;
	&get_overlap_info("$infile.region.overlapsupport.txt",\%support_hash,"gene",1);
	&add_support($infile,"$infile.addsupport",\%support_hash,"Support");
}

sub predictORF
{
	my ($gff,$seq,$eff)=@_;
	my %hash;
	$eff ||= 0.8;
	my $Length=length($seq)*$eff;
	open (IN,$gff) || die $!;
	my $id=0;
	while(<IN>)
	{
		s/\s*$//;
		my @t=split /\t+/,$_;
		next if ($_ eq "" || @t<8);
		if ($t[2]=~/gene/i)
		{
			$id++;
		}
		elsif ($t[2]=~/exon/i || $t[2]=~/UTR/i)
		{
			push @{$hash{$t[0]}{$id}},[$t[6],$t[3],$t[4]];
		}
	}
	close IN;
	my ($Chr,$Seq) = ("","");
	my %ORF=();
	open (IN,$seq) || die $!;
	while(<IN>)
	{
		s/\s*$//;
		if (/^>/ || eof)
		{
			my $name=$1 if ($_=~/^\>(\S+)/);
			if ($_ ne "" && $_!~/^\>/)
			{
				$Seq.=uc($_);
			}
			if ($Chr ne "" && $Seq ne "" && exists $hash{$Chr})
			{
				my $hash_p=$hash{$Chr};
				foreach my $gene(keys %$hash_p)
				{
					my $exon_seq="";
					my @array=@{$hash{$Chr}{$gene}};
					if (${$hash_p->{$gene}}[0][0] eq "+")
					{
						@array=sort {$a->[1]<=>$b->[1]}@array;
						for (my $i=0;$i<@array;$i++)
						{
							my ($S,$E)=($array[$i][1],$array[$i][2]);
							my $tmp=substr($Seq,$S-1,$E-$S+1);
							$exon_seq.=$tmp;
						}
					}
					else
					{
						@array=sort {$b->[1]<=>$a->[1]}@array;
						for (my $i=0;$i<@array;$i++)
						{
							my ($S,$E)=($array[$i][1],$array[$i][2]);
							my $tmp=reverse(substr($Seq,$S-1,$E-$S+1));
							$tmp=~ tr/ACGTacgt/TGCAtgca/;
							$exon_seq.=$tmp;
						}
					}
					my ($ORF_start,$ORF_end)=scanORF($exon_seq,${$hash_p->{$gene}}[0][0]) if ($exon_seq ne "");
					($ORF_start,$ORF_end)=($ORF_start >= $ORF_end) ? ($ORF_end,$ORF_start): ($ORF_start,$ORF_end);
					my $pre_len=0;
					my ($start_codon_site,$stop_codon_site)=(-1,-1);
					for (my $i=0;$i<@array;$i++)
					{
						if ($pre_len+$array[$i][2]-$array[$i][1]+1>=$ORF_start && $pre_len<=$ORF_start)
						{
							$start_codon_site=$array[$i][1]+$ORF_start-$pre_len;
						}
						if ($pre_len+$array[$i][2]-$array[$i][1]+1>=$ORF_end && $pre_len<=$ORF_end)
						{
							$stop_codon_site=$array[$i][1]+$ORF_end-$pre_len;
						}
						$pre_len+=$array[$i][2]-$array[$i][1]+1;
					}
					if ($start_codon_site != -1 && $stop_codon_site != -1 && abs($stop_codon_site-$start_codon_site)>$Length)
					{
						$ORF{$gene}{start}=$start_codon_site; ##(${$hash_p->{$gene}}[0][0] eq "-")?$stop_codon_site:$start_codon_site;
						$ORF{$gene}{end}=$stop_codon_site; ##${$hash_p->{$gene}}[0][0] eq "-")?$start_codon_site:$stop_codon_site;
						if ( (${$hash_p->{$gene}}[0][0] eq "+" && $ORF{$gene}{end}<$ORF{$gene}{start}) ||
						   (${$hash_p->{$gene}}[0][0] eq "-" && $ORF{$gene}{start}<$ORF{$gene}{end}) )
						{
							$ORF{$gene}{start}=$stop_codon_site;
							$ORF{$gene}{end}=$start_codon_site;
						}
					}
				}
			}
			$Seq="";
			$Chr=$name;
		}
		else
		{
			$Seq.=uc($_);
		}
	}
	close IN;
	return %ORF;
}

sub scanORF
{
	my ($seq,$strand)=@_;
	my $pos=0;
	my $Length = length($seq)/2;
	my $pre_end=0;
	my ($start_site,$end_site)=(0,0);
	while($seq=~/$start_codon/g)
	{
		my $max_len=0;
		my $min_len;
		my $S=index($seq,$start_codon,$pos);
		$pos=$S+1;
		my $orf_seq="";
		foreach my $end(@stop_codon)
		{
			my $tmp=substr($seq,$S,length($seq)-$S);
			my $orf;
			while($tmp=~/^($start_codon\w*$end)/g)
			{
				$orf=$1;
			}
			next if (!defined $orf);
			my $orf_len=length($orf);
			if (($orf_len>=$Length)&&(length($orf)%3==0))
			{
				next if ($S+$orf_len-1<$pre_end);
				($max_len,$orf_seq)=($orf_len>$max_len)?($orf_len,$orf):($max_len,$orf_seq);
			}
		}
		if ($orf_seq ne "")
		{
			$pre_end=$S+$max_len;
			my ($s_site,$e_site)=check_orf($orf_seq);
			if (abs($e_site-$s_site)>abs($end_site-$start_site))
			{
				if ($strand eq "+")
				{
					($start_site,$end_site)=($S+$s_site,$S+$e_site);
				}
				else
				{
					($start_site,$end_site)=(length($seq)-$S-($e_site-$s_site+1),length($seq)-$S-1+$s_site);
				}
			}
		}
	}
	return ($start_site,$end_site);
}

sub check_orf
{
	my $Seq=shift;
	my @array=();
	for my $i(0..2)
	{
		my $shift_Seq=substr($Seq,$i,length($Seq)-$i);
		my ($S,$E)=(-1,-1);
		my $an=0;
		while ($shift_Seq=~/(\w{3})/g)
		{
			my $aa=$1;
			if (uc($aa) eq $start_codon)
			{
				$S=$an if ($S==-1);
			}
			else
			{
				foreach my $end (@stop_codon)
				{
					if ((uc($aa) eq $end)&&($an>$S))
					{
						$E=($E!=-1)?$E:$an;
						push @array,[$S+$i,$E+$i+2] unless (($S==-1)||($E==-1)||($S>$E));
					}
				}
			}
			$an+=3;
		}
	}
	my @combine=combine_overlap(\@array);
	my ($codon_start,$codon_end,$max_len)=(0,0,0);
	for (my $i=0;$i<@combine;$i++)
	{
		if ($combine[$i][1]-$combine[$i][0]>$max_len)
		{
			$max_len=$combine[$i][1]-$combine[$i][0];
			($codon_start,$codon_end)=($combine[$i][0],$combine[$i][1]);
		}
	}
	return ($codon_start,$codon_end);
}

sub combine_overlap
{
	my $array_p=shift;
	my @array=sort{$a->[0]<=>$b->[0]}@$array_p;
	my @combine=();
	my ($S,$E)=($array[0][0],$array[0][1]);
	for (my $i=1;$i<@array;$i++)
	{
		if ($array[$i][1]>$E)
		{
			push @combine,[$S,$E];
			($S,$E)=($array[$i][0],$array[$i][1]);
		}
		if ($i==@array-1)
		{
			if (($S==$array[$i][0])&&($E==$array[$i][1]))
			{
				push @combine,[$S,$E];
			}
		}
	}
	@combine=@array if (@combine==0);
	return @combine;
}

sub parse_table
{
	my ($file,$table,$column)=@_;
	($$table,$$column)=split /\:/,$file;
	$$column ||= "0-1-2-3";
}

sub get_region
{
	my ($input,$output,$region)=@_;
	my @column=split /[\:\,]/,$region;
	my $cmd="";
	my $region_column="";
	if (@column==5)
	{
		$cmd=qq(awk 'BEGIN{FS="[\\t;]";OFS="\\t"}(\$$column[0] != ""){print \$$column[0],\$$column[1],\$$column[2],\$$column[3],\$$column[4]}' $input > $output);
		$region_column="0,1-2-3-4";
	}
	elsif (@column==4)
	{
		$cmd=qq(awk 'BEGIN{FS="[\\t;]";OFS="\\t"}(\$$column[0] != ""){print \$$column[0],\$$column[1],\$$column[2],\$$column[3]}' $input > $output);
		$region_column="0-1-2-3";
	}
	system $cmd;
	return $region_column;
}

sub overlap
{
	my ($i1,$f1,$i2,$f2,$ol,$out)=@_;
	my $cmd=qq(disposeOverlap.pl --E O --i1 $i1 --f1 $f1 --i2 $i2 --f2 $f2 --OL $ol -M 1,gene > $out);
	system $cmd;
	#my $sed=qq(sed 's/\:[\+\-]//' $out > tmp && mv tmp $out);
	#system $sed;
}

sub get_overlap_info
{
	my ($file,$hash,$pattern,$trim)=@_;
#W14_V7.2.ctg000004:+	ID=evm.TU.W14_V7.2.ctg000004.1;	6359	8787	2429	2	6	0.00247015232606011	start_codon:6359,6361,3,3	stop_codon:8785,8787,3,3
	open (IN,$file)||die $!;
	while(<IN>)
	{
		chomp;
		my @t=split /\t/,$_;
		next if ((defined $pattern)&&($t[1] !~ /$pattern/i));
		$t[0]=~s/[\:\|][\+\-]$// if (defined $trim);
		$t[0]=~s/\|/\:/;
		my @arr=();
		for (my $i=8;$i<@t;$i++)
		{
			push @arr,[$1,$2,$3] if ($t[$i]=~/([^\:\|]+)[\:\|](\d+)\,(\d+)/);
		}
		@{$$hash{"$t[0]:$t[2]-$t[3]"}}=($t[0]=~/[\:\|]\+$/)?(sort{$a->[1]<=>$b->[1]}@arr):(sort{$b->[1]<=>$a->[1]}@arr);
	}
	close IN;
}

sub add_support
{
	my ($infile,$outfile,$hash,$mark)=@_;
	open (IN,$infile) || die $!;
	open (OUT,">$outfile") || die $!;
	while(<IN>){
		chomp;
		my @t=split /\t/,$_;
		if ((@t>2)&&($t[2] =~ /gene/i)){
			s/\s*$//;
			$_.=";" if ($_ !~ /\;$/);
			if (exists $$hash{"$t[0]:$t[3]-$t[4]"}){
				my @arr=@{$$hash{"$t[0]:$t[3]-$t[4]"}};
				my $addinfo="";
				foreach my $i(@arr)
				{
					$addinfo.="$$i[0]:$$i[1]-$$i[2],";
				}
				$addinfo=~s/\,$//;
				print OUT ("$_ $mark=$addinfo\n");
			}
			else{
				print OUT "$_\n";
			}
		}
		else{
			print OUT "$_\n";;
		}
	}
	close IN;
}

sub add_structure
{
	my ($infile,$outfile,$hash)=@_;
	my @pre=();
	my @addsite=();
	open (IN,$infile)||die $!;
	open (OUT,">$outfile") || die $!;
	while(<IN>)
	{
		chomp;
		my @t=split /\t/,$_;
		if ((@t<8)||(eof))
		{
			if (@addsite > 0)
			{
				for (my $i=0;$i<@addsite;$i++)
				{
					if ($pre[6] eq "+")
					{
						my @out=@pre;
						$out[1]="BENM";
						$out[2]=$1 if ($addsite[$i][0]=~/^([^\:\|]+)/);
						$out[3]=$addsite[$i][1];
						$out[4]=$addsite[$i][2];
						$out[8]=~s/(ID\=[^\;]+)\.([^\.\;]+);/$1\.$addsite[$i][0]\;/;
						print OUT ((join "\t",@out),"\n");
					}
					else
					{
						my @out=@pre;
						$out[1]="BENM";
						$out[2]=$1 if ($addsite[$i][0]=~/^([^\:\|]+)/);
						$out[3]=$addsite[$i][1];
						$out[4]=$addsite[$i][2];
						$out[8]=~s/(ID\=[^\;]+)\.([^\.\;]+);/$1\.$addsite[$i][0]\;/;
						print OUT ((join "\t",@out),"\n");
					}
				}
			}
			print OUT "$_\n";
			next;
		}
		elsif (($t[2] eq "gene")||($t[2] eq "mRNA"))
		{
			@addsite=();
			@addsite=@{$$hash{"$t[0]:$t[6]:$t[3]-$t[4]"}} if (exists $$hash{"$t[0]:$t[6]:$t[3]-$t[4]"});
		}
		else
		{
			if (@addsite>0)
			{
				my $j=0;
				for (my $i=0;$i<@addsite;$i++)
				{
					if ($t[6] eq "+")
					{
						if ($addsite[$i][2]<$t[4])
						{
							my @out=@t;
							$out[1]="BENM";
							$out[2]=$1 if ($addsite[$i][0]=~/^([^\:\|]+)/);
							$out[3]=$addsite[$i][1];
							$out[4]=$addsite[$i][2];
							$out[8]=~s/(ID\=[^\;]+)\.([^\.\;]+);/$1\.$addsite[$i][0]\;/;
							print OUT ((join "\t",@out),"\n");
							$j++;
						}
					}
					else
					{
						if ($addsite[$i][1]>$t[3])
						{
							my @out=@t;
							$out[1]="BENM";
							$out[2]=$1 if ($addsite[$i][0]=~/^([^\:\|]+)/);
							$out[3]=$addsite[$i][1];
							$out[4]=$addsite[$i][2];
							$out[8]=~s/(ID\=[^\;]+)\.([^\.\;]+);/$1\.$addsite[$i][0]\;/;
							print OUT ((join "\t",@out),"\n");
							$j++;
						} 
					}
				}
				if ($j>0)
				{
					for(0..($j-1))
					{
						shift @addsite;
					}
				}
			}
		}
		my @tmp=@t;
		$tmp[1]="BENM";
		print OUT ((join "\t",@tmp),"\n");
		@pre=@t;
	}
	close IN;
	close OUT;
}

sub check_structure
{
	my ($infile,$outfile,$region)=@_;
	my @column=split /[\:\,\-]/,$region;
	my $cmd=qq(awk '{if(\$$column[1]!=""&&\$$column[1]!="gene"&&\$$column[1]!="CDS"&&\$$column[1]!="exon"&&\$$column[1]!="mRNA"){print \$$column[0]","\$$column[1]","\$$column[2]","\$$column[3]}}' $infile |sort |uniq -c |awk '\$1>1{print \$2,\$1}' > $Dir/structure_site_split_dup.txt);
	system $cmd;
	my %hash;
	open (IN,"$Dir/structure_site_split_dup.txt") || die $!;
	while(<IN>)
	{
		chomp;
		s/^\s+//;
		my @t=split /\s+/,$_;
		$hash{$t[0]}=$t[1];
	}
	close IN;
	my @gffinfo=();
	#my @gffcombine=();
	my $add=0;
	my $k=1;
	my ($pre_chr,$pre_strand)=("","");
	my %dup;
	open (IN,$infile) || die $!;
	open (OUT,">$outfile") || die $!;
	while(<IN>)
	{
		s/\s*$//;
		my @t=split /\t/,$_;
		if ((@t<8)&&(!eof))
		{
			next;
		}
		if ((eof)||((@t>=8)&&($t[2] eq "gene")))
		{
			push @gffinfo,@t if ((eof)&&(@t>=8));
			if (((eof)||($add==0)||(($pre_chr ne "")&&(($pre_chr ne $t[0]))||(($pre_strand ne "")&&($pre_strand ne $t[6]))))&&(@gffinfo>0))
			{
				my ($chr,$S,$E,$strand,$score)=("","","","",".");
				my @info=();
				my %dup=();
				my $start_loc=0;
				my $end_loc=1;
				my %frame_count=();
				foreach my $i(@gffinfo)
				{
					if (($$i[2]=~/gene/i)||($$i[2]=~/mRNA/i))
					{
						$S=(($S eq "")||($$i[3]<$S))?$$i[3]:$S;
						$E=(($E eq "")||($$i[4]>$E))?$$i[4]:$E;
						die ("strand error: $strand\t",(join "\t",@$i),"\n") if ((defined $strand)&&($strand ne "")&&($strand ne $$i[6]));
						die ("chr-id error: $chr\t",(join "\t",@$i),"\n") if ((defined $chr)&&($chr ne "")&&($chr ne $$i[0]));
						$chr=$$i[0];
						$strand=$$i[6];
						$score+=$$i[7] if (($$i[2]=~/gene/i)&&($$i[7] ne "."));
					}
					else
					{
						$start_loc = $$i[3] if ($$i[2]=~/start/i);
						$end_loc   = $$i[3] if ($$i[2]=~/stop/i);
						my $key=join ",",@$i[0,2,3,4];
						push @info,[@$i[0..7]] unless (exists $dup{$key});
						$dup{$key}++;
						$frame_count{$$i[2]}++;
					}
				}
				if ((($strand eq "+")&&($start_loc<$end_loc))||(($strand eq "-")&&($start_loc>$end_loc)))
				{
					@info=($strand eq "+")?(sort{$a->[3]<=>$b->[3]}@info):(sort{$b->[3]<=>$a->[3]}@info);
					print OUT "$chr\tBENM\tgene\t$S\t$E\t\.\t$strand\t$score\tID\=fbx.$chr\.$k.gene; Name\=fbx.$chr\.$k\;\n";
					print OUT "$chr\tBENM\tmRNA\t$S\t$E\t\.\t$strand\t$score\tID\=fbx.$chr\.$k.mrna; Parent\=fbx.$chr\.$k\;\n";
					my %infocount;
					foreach my $g(@info)
					{
						next if (($$g[2]=~/gene/i)||($$g[2]=~/mRNA/i));
						if ($$g[6] eq "+")
						{
							$infocount{$$g[2]}++;
						}
						else
						{
							if ((!exists $infocount{$$g[2]})||($infocount{$$g[2]}==0))
							{
								$infocount{$$g[2]}=$frame_count{$$g[2]};
							}
							$infocount{$$g[2]}--;
							$$g[7]=$infocount{$$g[2]} if ($$g[2]=~/CDS/i);
						}
						print OUT ((join "\t",@$g),"\tID\=fbx.$chr\.$k.".lc($$g[2]).".".$infocount{$$g[2]}."; Parent\=fbx.$chr\.$k\;\n");
					}
					$k++;
					print OUT "\n";
				}
				else
				{
					if (defined $Verbose)
					{
						print STDERR "annotation error!: start code, end code in the wrong strand\n";
						print STDERR "$chr\tBENM\tgene\t$S\t$E\t\.\t$strand\t$score\tID\=fbx.$chr\.$k.gene; Name\=fbx.$chr\.$k\;\n";
						print STDERR "$chr\tBENM\tmRNA\t$S\t$E\t\.\t$strand\t$score\tID\=fbx.$chr\.$k.mrna; Parent\=fbx.$chr\.$k\;\n";
					}
				}
				@gffinfo=();
			}
			$add=0;
			$k=1 if ((!eof)&&($pre_chr ne "")&&($pre_chr ne $t[0]));
			push @gffinfo,[@t];
		}
		else
		{
			my $key=join ",",@t[0,2,3,4];
			if (exists $hash{$key})
			{
				$hash{$key}--;
				if ($hash{$key}>0)
				{
					$add=1;
				}
			}
			push @gffinfo,[@t];
		}
		$pre_chr=$t[0];
		$pre_strand=$t[6];
	}
	close IN;
	close OUT;
}

sub identify_utr
{
	my ($infile,$outfile)=@_;
	open (IN,$infile) || die $!;
	open (OUT,">$outfile") || die $!;
	my @Gene=();
	my ($fs,$fe,$ts,$te)=("","","",""); #5-UTR 3-UTR
	while(<IN>)
	{
		chomp;
		my @t=split /\t/,$_;
		if ((@t<8)||(eof))
		{
			if (($fs eq "")&&($ts eq ""))
			{
				foreach my $info(@Gene)
				{
					print OUT ((join "\t",@$info),"\n");
				}
				@Gene=();
			}
			else
			{
				my %infocount=();
				foreach my $info(@Gene)
				{
					if (($$info[2]=~/gene/i)||($$info[2]=~/mRNA/i))
					{
						print OUT ((join "\t",@$info),"\n");
					}
					else
					{
						if (($$info[2]=~/exon/i)||($$info[2]=~/CDS/i))
						{
							if ($fs ne "")
							{
								if ($$info[6] eq "+")
								{
									if ($$info[4]<$fs)
									{
										next if ($$info[2]=~/CDS/i);
										my @tmp=@$info;
										my $mark=lc($tmp[2]);
										$tmp[2]="5-UTR";
										$infocount{$tmp[2]}++;
										#"fbx.".$tmp[0].".$k.".lc($tmp[2]);
										my $subinfo="5-UTR.".$infocount{$tmp[2]};
										$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
										print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
									}
									elsif (($$info[4]>$fs)&&($$info[3]<=$fe))
									{
										next if ($$info[2]=~/CDS/i);
										my @tmp=@$info;
										$tmp[4]=$fs-1;
										my $mark=lc($tmp[2]);
										$tmp[2]="5-UTR";
										$infocount{$tmp[2]}++;
										my $subinfo="5-UTR.".$infocount{$tmp[2]};
										$tmp[8]=~s/ID\=([^\;]+)$mark\d+\;/ID\=$1$subinfo\;/i;
										print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
										
										@tmp=@$info;
										$tmp[3]=$fs;
										$infocount{$tmp[2]}++;
										$subinfo=$infocount{$tmp[2]};
										$tmp[8]=~s/ID\=([^\;]+\d+\.$mark)\d+\;/ID\=$1$subinfo\;/i;
										print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
									}
									else
									{
										if ($ts eq "")
										{
											my @tmp=@$info;
											my $mark=lc($tmp[2]);
											$infocount{$tmp[2]}++;
											my $subinfo=$infocount{$tmp[2]};
											$tmp[8]=~s/ID\=([^\;]+\d+\.$mark)\d+\;/ID\=$1$subinfo\;/i;
											print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
										}
									}
								}
								else
								{
									if ($$info[3]>$fe)
									{
										next if ($$info[2]=~/CDS/i);
										my @tmp=@$info;
										my $mark=lc($tmp[2]);
										$tmp[2]="5-UTR";
										$infocount{$tmp[2]}++;
										#"fbx.".$tmp[0].".$k.".lc($tmp[2]);
										my $subinfo="5-UTR.".$infocount{$tmp[2]};
										$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
										print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
									}
									elsif (($$info[3]<$fe)&&($$info[4]>=$fs))
									{
										next if ($$info[2]=~/CDS/i);
										my @tmp=@$info;
										$tmp[3]=$fe+1;
										my $mark=lc($tmp[2]);
										$tmp[2]="5-UTR";
										$infocount{$tmp[2]}++;
										my $subinfo="5-UTR.".$infocount{$tmp[2]};
										$tmp[8]=~s/ID\=([^\;]+)$mark\d+\;/ID\=$1$subinfo\;/i;
										print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
										
										@tmp=@$info;
										$tmp[4]=$fe;
										$infocount{$tmp[2]}++;
										$subinfo=$infocount{$tmp[2]};
										$tmp[8]=~s/ID\=([^\;]+\d+\.$mark)\d+\;/ID\=$1$subinfo\;/i;
										print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
									}
									else
									{
										if ($ts eq "")
										{
											my @tmp=@$info;
											my $mark=lc($tmp[2]);
											$infocount{$tmp[2]}++;
											my $subinfo=$infocount{$tmp[2]};
											$tmp[8]=~s/ID\=([^\;]+\d+\.$mark)\d+\;/ID\=$1$subinfo\;/i;
											print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
										}
									}
								}
							}
							if ($ts ne "")
							{
								if ($$info[6] eq "+")
								{
									if ($$info[3]>$te)
									{
										next if ($$info[2]=~/CDS/i);
										my @tmp=@$info;
										my $mark=lc($tmp[2]);
										$tmp[2]="3-UTR";
										$infocount{$tmp[2]}++;
										#"fbx.".$tmp[0].".$k.".lc($tmp[2]);
										my $subinfo="3-UTR.".$infocount{$tmp[2]};
										$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
										print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
									}
									elsif (($$info[3]<$te)&&($$info[4]>=$ts))
									{
										next if ($$info[2]=~/CDS/i);
										my @tmp=@$info;
										my $mark=lc($tmp[2]);
										$tmp[4]=$te;
										$infocount{$tmp[2]}++;
										my $subinfo=$infocount{$tmp[2]};
										$tmp[8]=~s/ID\=([^\;]+\d+\.$mark)\d+\;/ID\=$1$subinfo\;/i;
										print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
										
										@tmp=@$info;
										$tmp[3]=$te+1;
										if ($tmp[4]>$tmp[3])
										{
											$tmp[2]="3-UTR";
											$infocount{$tmp[2]}++;
											$subinfo="3-UTR.".$infocount{$tmp[2]};
											$tmp[8]=~s/ID\=[^\;]+($mark\d+)\;/ID\=$subinfo\;/i;
											print OUT ((join "\t",@tmp),"\n");
										}
									}
									else
									{
										if ($fs eq "")
										{
											my @tmp=@$info;
											my $mark=lc($tmp[2]);
											$infocount{$tmp[2]}++;
											my $subinfo=$infocount{$tmp[2]};
											$tmp[8]=~s/ID\=([^\;]+\d+\.$mark)\d+\;/ID\=$1$subinfo\;/i;
											print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
										}
									}
								}
								else
								{
									if ($$info[4]<$ts)
									{
										next if ($$info[2]=~/CDS/i);
										my @tmp=@$info;
										my $mark=lc($tmp[2]);
										$tmp[2]="3-UTR";
										$infocount{$tmp[2]}++;
										#"fbx.".$tmp[0].".$k.".lc($tmp[2]);
										my $subinfo="3-UTR.".$infocount{$tmp[2]};
										$tmp[8]=~s/ID\=([^\;]+\d+\.)$mark\d+\;/ID\=$1$subinfo\;/i;
										print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
									}
									elsif (($$info[4]>$ts)&&($$info[3]<=$te))
									{
										next if ($$info[2]=~/CDS/i);
										my @tmp=@$info;
										my $mark=lc($tmp[2]);
										$tmp[3]=$ts;
										$infocount{$tmp[2]}++;
										my $subinfo=$infocount{$tmp[2]};
										$tmp[8]=~s/ID\=([^\;]+\d+\.$mark)\d+\;/ID\=$1$subinfo\;/i;
										print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
										
										@tmp=@$info;
										$tmp[4]=$ts-1;
										if ($tmp[4]>$tmp[3])
										{
											$tmp[2]="3-UTR";
											$infocount{$tmp[2]}++;
											$subinfo="3-UTR.".$infocount{$tmp[2]};
											$tmp[8]=~s/ID\=[^\;]+($mark\d+)\;/ID\=$subinfo\;/i;
											print OUT ((join "\t",@tmp),"\n");
										}
									}
									else
									{
										if ($fs eq "")
										{
											my @tmp=@$info;
											my $mark=lc($tmp[2]);
											$infocount{$tmp[2]}++;
											my $subinfo=$infocount{$tmp[2]};
											$tmp[8]=~s/ID\=([^\;]+\d+\.$mark)\d+\;/ID\=$1$subinfo\;/i;
											print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
										}
									}
								}
							}
							if (($fs ne "")&&($te ne ""))
							{
								my ($s,$e)=(sort{$a<=>$b}($fs,$fe,$ts,$te)[1,2]);
								if (($$info[4]<=$e)&&($$info[3]>=$s))
								{
									my @tmp=@$info;
									my $mark=lc($tmp[2]);
									$infocount{$tmp[2]}++;
									my $subinfo=$infocount{$tmp[2]};
									$tmp[8]=~s/ID\=([^\;]+\d+\.$mark)\d+\;/ID\=$1$subinfo\;/;
									print OUT ((join "\t",@tmp),"\n") if ($tmp[4]>$tmp[3]);
								}
							}
						}
						else
						{
							print OUT ((join "\t",@$info),"\n");
						}
					}
				}
				@Gene=();
			}
			print OUT "\n";
			($fs,$fe,$ts,$te)=("","","","");
		}
		else
		{
			push @Gene,[@t];
			($fs,$fe)=($t[3],$t[4]) if ($t[2] =~ /start/);
			($ts,$te)=($t[3],$t[4]) if ($t[2] =~ /stop/);
		}
	}
	close IN;
	close OUT;
}

sub check_gff
{
	my ($infile,$outfile)=@_;
	my @info=();
	my %hash;
	my %site;
	open (IN,$infile) || die $!;
	while(<IN>) {
		if (/CDS/||/exon/){
			$hash{$1}=1 if (/Parent\=([^\;]+)/);
		}
		if (/start/||/stop/){
			$site{$1}++ if (/Parent\=([^\;]+)/);
		}
	}
	close IN;
	
	my %count;
	open (IN,$infile) || die $!;
	open (OUT,">$outfile") || die $!;
	while(<IN>) {
		chomp;
		if (/Parent\=([^\;]+)/){
			my $id=$1;
			if (!exists $hash{$id} || !exists $site{$id} || $site{$1}<2){
				my @t=split /\t+/,$_;
				my $s="exon";
				$count{$id}{$s}++;
				$t[8]=~s/$t[2]\d+/$s$count{$id}{$s}/g;
				if ($t[2]=~/five/ || $t[2]=~ /5/){
					if ($t[6] eq "+"){
						$t[4]+=1;
					}else{
						$t[3]-=1;
					}
				}elsif(/three/ || $t[2]=~ /3/){
					if ($t[6] eq "+"){
						$t[3]-=1;
					}else{
						$t[4]+=1;
					}
				}
				$t[2]=$s;
				print OUT ((join "\t",@t),"\n");
			}
			else{
				print OUT "$_\n";
			}
		}
		else{
			print OUT "$_\n";
		}
	}
	close IN;
	close OUT;
}

sub filter_repeat
{
	my ($infile,$outfile)=@_;
	my @arr=();
	my $print=1;
	my $print_out="";
	open (IN,$infile) || die $!;
	open (OUT,">$outfile") || die $!;
	while(<IN>)
	{
		my @t=split /\t/,$_;
		if (@t<8)
		{
			$print_out.=$_;
			next;
		}
		if ($t[2] eq "gene")
		{
			@arr=();
			print OUT $print_out if (($print==1)&&($print_out ne ""));
			$print_out="";
			my $repeat=$1 if ($t[8]=~/Repeat\=([^\;]+)/);
			if (defined $repeat)
			{
				while($repeat=~/\:(\d+)\-(\d+)/g)
				{
					push @arr,[$1,$2];
				}
			}
			$print=1;
		}
		$print_out.=$_;
		if ($t[2] eq "exon")
		{
			for (my $i=0;$i<@arr;$i++)
			{
				my $overlap=overlap_size([$t[3],$t[4]],[$arr[$i][0],$arr[$i][1]]);
				$print=0 if (($overlap/($t[4]-$t[3]+1)>0)||($overlap>=30));
			}
		}	
	}
	print OUT $print_out if (($print==1)&&($print_out ne ""));
	close IN;
	close OUT;
}

sub overlap_size
{
	my ($block1_p,$block2_p) = @_;
	my $combine_start = ($block1_p->[0] < $block2_p->[0]) ?  $block1_p->[0] : $block2_p->[0];
	my $combine_end = ($block1_p->[1] > $block2_p->[1]) ?  $block1_p->[1] : $block2_p->[1];
	my $overlap_size = ($block1_p->[1]-$block1_p->[0]+1) + ($block2_p->[1]-$block2_p->[0]+1) - ($combine_end-$combine_start+1);
	return $overlap_size;
}

__END__
