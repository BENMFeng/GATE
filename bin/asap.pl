#!/usr/bin/perl -w
##########################################################################
#  Copyright (c) 2012-2013 - BENM(Binxiao) Feng                          #
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
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;


=head1 Name

asap.pl -- Alternative Splicing Application integrative Program

=head1 Version

Author: BENM <binxiaofeng@gmail.com>

Version: v1.7 alpha, Mar 21th, 2013

=head1 Option

  --ref <file>		input refGene (unique gene) GFF/GTF annotation file
  --tran <file>		input transcripts ab inito annotation structure file
  --bam <file>		input RNA-Seq alignment to genome BAM format file
  --fpkm <float>	set FPKM minimuim value, default: 0.1
  --cutoff <int>	set junction reads supported cut-off, default: 6
  --version		show update info
  --verbose		print out run status
  --help		show this info

=head1 Example

$ perl bin/asap.pl -ref refGene.gff -tran transcripts.gtf \
  -bam accepted_hits.bam -cutoff 10

=cut

my %opts;
my ($RefGff,$TranGff,$Bam,$Cutoff,$FPKM,$Verbose,$Version,$Help);
GetOptions(
	\%opts,
	"ref:s"=>\$RefGff,
	"tran:s"=>\$TranGff,
	"bam:s"=>\$Bam,
	"fpkm:s"=>\$FPKM,
	"cutoff:i"=>\$Cutoff,
	"version"=>\$Version,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
my $updateLog=qq(
asap.pl -- Alternative Splicing Application integrative Program
Latest Version: 1.4 alpha, Oct 10th, 2012
1.0 alpha, Sep 25th, 2012
1.1 alpha, Sep 27th, 2012
1.2 alpha, Sep 29th, 2012 on the train to LZ
1.3 alpha, Oct 8th, 2012
1.4 alpha, Oct 10th, 2012
1.5 alpha, Oct 17th, 2012
);
my $instruction=qq(
);

die $updateLog if (defined $Version);
die `pod2text $0` if ((!defined $RefGff)||(!defined $TranGff)||($Help));

$Cutoff ||= 6;
$FPKM ||= 0.1;

my %ASevent=();
my %ASstat=();
my $withchr="";
my %TranscriptStruct=parseTranscript($TranGff);
my %GeneStruct=getStructure($RefGff,\%TranscriptStruct);
detectASevent($Bam,\%GeneStruct,\%TranscriptStruct);


sub parseTranscript {
	my $file = shift;
#chr1	Cufflinks	transcript	34611	36081	1	-	.	gene_id "FAM138F"; transcript_id "NR_026820"; FPKM "0.0000000000"; frac "0.000000"; conf_lo "0.000000"; conf_hi "0.000000"; cov "0.000000";
#chr1	Cufflinks	exon	34611	35174	1	-	.	gene_id "FAM138F"; transcript_id "NR_026820"; exon_number "1"; FPKM "0.0000000000"; frac "0.000000"; conf_lo "0.000000"; conf_hi "0.000000"; cov "0.000000";
#chr1	Cufflinks	exon	35277	35481	1	-	.	gene_id "FAM138F"; transcript_id "NR_026820"; exon_number "2"; FPKM "0.0000000000"; frac "0.000000"; conf_lo "0.000000"; conf_hi "0.000000"; cov "0.000000";
#chr1	Cufflinks	exon	35721	36081	1	-	.	gene_id "FAM138F"; transcript_id "NR_026820"; exon_number "3"; FPKM "0.0000000000"; frac "0.000000"; conf_lo "0.000000"; conf_hi "0.000000"; cov "0.000000";
	my %transcripts=();
	open (IN,$file) || die $!;
	$_=<IN>;
	$withchr="chr" if (/^chr/);
	$_=<IN>;
	$withchr="chr" if (/^chr/);
	$_=<IN>;
	$withchr="chr" if (/^chr/);
	close IN;
	open (IN,$file) || die $!;
	my $i=0;
	my $skip=0;
	while(<IN>) {
		chomp;
		my @col = split /\t+/,$_;
		my ($chr,$id,$start,$end,$strand,$info)=@col[0,2..4,6,8];
		if ($id=~/transcript/) {
			if ($info=~/FPKM\s+\"([^\"]+)\"/) {
				if ($1<$FPKM) {
					$skip=1;
					next;
				} else {
					$skip=0;
				}
			}
			$i++;
			if ($info=~/transcript_id\s+\"([^\"\s]+)\"/ || $info=~/gene_id\s+\"([^\"\s]+)\"/ || $info=~/gene_name\s+\"([^\"]+)\"/) {
				push @{$transcripts{$chr}{$1}},[$i,$start,$end,$strand,$info];
				$transcripts{"$chr-$i"}=@{$transcripts{$chr}{$1}}-1;
			} else {
				push @{$transcripts{$chr}},[$i,$start,$end,$strand,$info];
				$transcripts{"$chr-$i"}=@{$transcripts{$chr}}-1;
			}
		}
		else {
			next if ($skip==1);
			if ($info=~/FPKM\s+\"([^\"\s]+)\"/) {
				if ($1>=$FPKM) {
					push @{$transcripts{$i}},[$start,$end];
				}
			}
		}
	}
	close IN;
	return %transcripts;
}

sub getStructure {
	my ($refGene,$transcripts) = @_;
#refGene.gff
#chr1    refGene mRNA    11874   14408   .       +       .       ID=NR_046018; name=DDX11L1;
#chr1    refGene exon    11874   12227   .       +       .       Parent=NR_046018;
#chr1    refGene exon    12613   12721   .       +       .       Parent=NR_046018;
#chr1    refGene exon    13221   14408   .       +       .       Parent=NR_046018;

#Alternative Splicing Type
#A3SS: alternative 3' splicing sites
#A5SS: alternative 5' splicing sites
#AFE: alternative first exons
#ALE: alternative last exons
#MXE: mutually exclusive exons
#RI: retained intron
#SE: skipped exons

	my %genes;
	my $n=0;
	my $count=0;
	my ($line,$mark1,$mark2)=(0,-1,-1);
	open (IN,$refGene) || die $!;
	while(<IN>) {
		chomp;
		next if ($_=~/^\#/ || $_ eq "");
		last if ($line>=100);
		$line++;
		my @col = split /\t+/,$_;
		if ($col[2] =~ /mRNA/i) {
			$mark1=0;
		} elsif ($col[2] =~ /gene/i && $mark1==0) {
			$mark1=1;
		}
		if ($col[2] =~ /exon/i) {
			$mark2=0;
		} elsif ($col[2] =~ /CDS/i && $mark2==0) {
			$mark2=1;
		}
	}
	close IN;
	my $pre_gene_id="";
	my ($pre_chr,$pre_name,$coor1,$coor2,$pre_strand,$pre_info)=("","","","","","");
	open (IN,$refGene) || die $!;
	while(<IN>) {
		chomp;
		next if (/^\#/);
		my @col = split /\t+/,$_;
		my ($chr,$id,$start,$end,$strand,$info)=@col[0,2..4,6,8];
		$chr="$withchr$chr";
		if ($col[2]=~/mRNA/i || $col[2]=~/gene/i || $col[2]=~/transcript/i) {
			next if ($mark1==1 && $col[2]=~/gene/i);
			my $k=-1;
			my %K=();
			if ($info=~/gene_name\s+\"([^\"]+)\"/ || $info=~/gene_id\s+\"([^\"\s]+)\"/ || $info=~/transcript_id\s+\"([^\"\s]+)\"/
			 || $info=~/name\=([^\;\s]+)/i || $info=~/ID\=([^\;\s]+)/i) {
				my $name=$1;
				if (!exists $transcripts->{$chr}{$name})
				{
					$name = $1 if ($info=~/gene_id\s+\"([^\"\s]+)\"/);
					if (!exists $transcripts->{$chr}{$name})
					{
						$name = $1 if ($info=~/transcript_id\s+\"([^\"\s]+)\"/);
					}
				}
				if (exists $transcripts->{$chr}{$name}) {
					if ($info=~/gene_name\s+\"([^\"]+)\"/ || $info=~/name\=([^\;\s]+)/i || $info=~/ID\=([^\;\s]+)/i) {
						$info="name=$1;";
					} else {
						$info="name=$name;";
					}
					print STDERR "find gene: $name\n" if ($Verbose);
					my @tr=@{$transcripts->{$chr}{$name}};
					for (my $i=0;$i<@tr;$i++) {
						if ($tr[$i][3] ne $strand || $tr[$i][1]>$end || $tr[$i][2]<$start) {
							$k=-1;
							next;
						} else {
							$k=$tr[$i][0];
							$K{$k}=1 if (exists $transcripts->{$k} && defined $transcripts->{$k});
						}
					}
				} else {
					$name = $1 if ($info=~/ID\=([^\;\s]+)/i);
					if (exists $transcripts->{$chr}{$name}) {
						print STDERR "find gene: $name\n" if ($Verbose);
						my @tr=@{$transcripts->{$chr}{$name}};
						for (my $i=0;$i<@tr;$i++) {
							if ($tr[$i][3] ne $strand || $tr[$i][1]>$end || $tr[$i][2]<$start) {
								$k=-1;
								next;
							} else {
								$k=$tr[$i][0];
								$K{$k}=1 if (exists $transcripts->{$k} && defined $transcripts->{$k});
							}
						}
					}
					elsif (exists $transcripts->{$chr} && ref($transcripts->{$chr}) eq "ARRAY") {
						my @tr=@{$transcripts->{$chr}};
						for (my $i=0;$i<@tr;$i++) {
							if ($tr[$i][3] ne $strand || $tr[$i][1]>$end || $tr[$i][2]<$start) {
								$k=-1;
								next;
							} else {
								$k=$tr[$i][0];
								$K{$k}=1 if (exists $transcripts->{$k} && defined $transcripts->{$k});
							}
						}
					}
				}
			} else {
				my @tr=@{$transcripts->{$chr}};
				for (my $i=0;$i<@tr;$i++) {
					if ($tr[$i][3] ne $strand || $tr[$i][1]>$end || $tr[$i][2]<$start) {
						$k=-1;
						next;
					} else {
						$k=$tr[$i][0];
						$K{$k}=1 if (exists $transcripts->{$k} && defined $transcripts->{$k});
					}
				}
			}
			if (keys %K>0) {
				$n++;
				push @{$genes{$chr}},[$n,$start,$end,$strand,$info,\%K];
				$genes{'chr'}{$chr}++;
				$count=1;
			} else {
				$count=0;
			}
		}
		if ($mark1==-1)
		{
			if ($info=~/gene_name\s+\"([^\"]+)\"/ || $info=~/gene_id\s+\"([^\"\s]+)\"/ || $info=~/transcript_id\s+\"([^\"\s]+)\"/
			 || $info=~/name\=([^\;\s]+)/i || $info=~/ID\=([^\;\s]+)/i) {
				my $name=$1;
				if (!exists $transcripts->{$chr}{$name})
				{
					$name = $1 if ($info=~/gene_id\s+\"([^\"\s]+)\"/);
					if (!exists $transcripts->{$chr}{$name})
					{
						$name = $1 if ($info=~/transcript_id\s+\"([^\"\s]+)\"/);
					}
				}
				if ( ($pre_gene_id ne $name && exists $transcripts->{$chr}{$name}) || eof ) {
					if ($count==1)
					{
						if (eof)
						{
							($coor1,$coor2)=(sort{$a<=>$b}($coor1,$coor2,$start,$end))[0,3];
						}
						my $k=-1;
						my %K=();
						if (exists $transcripts->{$pre_chr}{$pre_name}) {
							print STDERR "find gene: $pre_name\n" if ($Verbose);
							my @tr=@{$transcripts->{$pre_chr}{$pre_name}};
							for (my $i=0;$i<@tr;$i++) {
								if ($tr[$i][3] ne $pre_strand || $tr[$i][1]>$coor2 || $tr[$i][2]<$coor1) {
									$k=-1;
									next;
								} else {
									$k=$tr[$i][0];
									$K{$k}=1 if (exists $transcripts->{$k} && defined $transcripts->{$k});
								}
							}
						}
						elsif (exists $transcripts->{$pre_chr} && ref($transcripts->{$pre_chr}) eq "ARRAY") {
							my @tr=@{$transcripts->{$pre_chr}};
							for (my $i=0;$i<@tr;$i++) {
								if ($tr[$i][3] ne $strand || $tr[$i][1]>$coor2 || $tr[$i][2]<$coor1) {
									$k=-1;
									next;
								} else {
									$k=$tr[$i][0];
									$K{$k}=1 if (exists $transcripts->{$k} && defined $transcripts->{$k});
								}
							}
						}
						push @{$genes{$pre_chr}},[$n,$coor1,$coor2,$pre_strand,$pre_info,\%K] if (keys %K>0);
					}
					$genes{'chr'}{$chr}++;
					$n++;
					$pre_chr=$chr;
					$pre_name=$name;
					($coor1,$coor2)=($start,$end);
					$pre_strand=$strand;
					$pre_info=$info;
					$pre_info=~s/exon\_number\s+\"[^\"\s]+\"\;\s//;
					$count=1;
				}
				elsif ($pre_gene_id eq $name)
				{
					if ($count==1)
					{
						($coor1,$coor2)=(sort{$a<=>$b}($coor1,$coor2,$start,$end))[0,3];
					}
				}
				if (!exists $transcripts->{$chr}{$name})
				{
					$count=0;
				}
				$pre_gene_id=$name;
			}
		}
		if ($col[2]=~/exon/i || $col[2]=~/CDS/i || $col[2]=~/UTR/i) {
			next if ($mark2==1 && $col[2]=~/CDS/i);
			if ($count==1) {
				push @{$genes{$n}},[$start,$end]; #unless (exists $genes{$n} && @{$genes{$n}}>0 && $start<${$genes{$n}}[-1][1] && $end>${$genes{$n}}[-1][0]);
			}
		}
	}
	close IN;
	return %genes;
}

sub detectASevent {
	my ($bam,$genes,$transcripts)=@_;
	my @chrAry=sort keys %{$genes->{'chr'}};
	foreach my $chr(@chrAry) {
		print STDERR "$chr\: $genes->{'chr'}{$chr} genes related with transcripts\n" if (defined $Verbose);
		for (my $i=0;$i<@{$genes->{$chr}};$i++) {
			my %geneStr=();
			my ($j,$start,$end,$strand,$info,$K_p)=@{${$genes->{$chr}}[$i]};
			my $name="";
			if ($info=~/gene_name\s+\"([^\"]+)\"/ || $info=~/gene_id\s+\"([^\"\s]+)\"/ || $info=~/transcript_id\s+\"([^\"\s]+)\"/
			 || $info=~/name\=([^\;\s]+)/i || $info=~/ID\=([^\;\s]+)/i) {
				$name=$1;
				if (!exists $transcripts->{$chr}{$name})
				{
					$name = $1 if ($info=~/gene_id\s+\"([^\"\s]+)\"/);
					if (!exists $transcripts->{$chr}{$name})
					{
						$name = $1 if ($info=~/transcript_id\s+\"([^\"\s]+)\"/);
					}
				}
			}
			my @gexon=@{$genes->{$j}};
			@{$geneStr{mRNA}}=($chr,$start,$end,$strand,$info);
			@{$geneStr{exon}}=($strand eq "+")?@gexon:sort{$b->[0]<=>$a->[0]}@gexon;
			foreach my $k(sort{$a<=>$b}keys %$K_p) {
				my %transcriptStr=();
				next if (!exists $transcripts->{"$chr-$k"} || !exists $transcripts->{$k});
				my $tr_no=$transcripts->{"$chr-$k"};
				my ($trstart,$trend,$trinfo);
				if (defined $name && exists $transcripts->{$chr}{$name}) {
					($trstart,$trend,$trinfo)=@{${$transcripts->{$chr}{$name}}[$tr_no]}[1,2,4];
				} elsif (ref($transcripts->{$chr}) eq "ARRAY") {
					($trstart,$trend,$trinfo)=@{${$transcripts->{$chr}}[$tr_no]}[1,2,4];
				}
				my @trexon=@{$transcripts->{$k}};
				@{$transcriptStr{mRNA}}=($start,$end,$strand,$info);
				@{$transcriptStr{exon}}=($strand eq "+")?@trexon:sort{$b->[0]<=>$a->[0]}@trexon;
				my $out=classifyAS($bam,\%geneStr,\%transcriptStr);
				print $out;
			}
		}
	}
}

sub classifyAS {
	my ($bam,$geneStr,$transcriptStr)=@_;
	checkIndex($bam);
	my $out="";
	my $chr=${$geneStr->{mRNA}}[0];
	my $strand=${$geneStr->{mRNA}}[3];
	my $info=${$geneStr->{mRNA}}[4];
	my @gexon=@{$geneStr->{exon}};
	my @texon=@{$transcriptStr->{exon}};
	$out.=detectA5SS($bam,$chr,$strand,$info,\@gexon,\@texon);
	$out.=detectAFE($bam,$chr,$strand,$info,\@gexon,\@texon);
	$out.=detectMXSE($bam,$chr,$strand,$info,\@gexon,\@texon);
	$out.=detectALE($bam,$chr,$strand,$info,\@gexon,\@texon);
	$out.=detectA3SS($bam,$chr,$strand,$info,\@gexon,\@texon);
	$out.=detectRI($bam,$chr,$strand,$info,\@gexon,\@texon);
	return $out;
}

sub detectA5SS {
	my ($bam,$chr,$strand,$info,$gexon,$texon)=@_;
	$info=~s/[\;\s]+$//;
	my $out="";
	if (@$gexon>=3 && @$texon>=3) {
		my ($fst_gs,$fst_ge)=@{$gexon->[0]};
		my ($sec_gs,$sec_ge)=@{$gexon->[1]};
		my ($fst_ts,$fst_te)=@{$texon->[0]};
		my ($sec_ts,$sec_te)=@{$texon->[1]};
		if ($fst_gs==$fst_ts && $fst_ge==$fst_te) {
			return "";
		} elsif ($fst_ts<=$fst_gs && $fst_te>=$fst_gs) {
			my ($jun1,$jun2)=(0,0);
			my ($S,$E)=(sort{$a<=>$b}($fst_gs,$fst_gs,$sec_gs,$sec_ge,$fst_ts,$fst_te,$sec_ts,$sec_te))[0,-1];
			my ($alts,$alte)=(0,0);
			open (IN,"samtools view $bam $chr:$S-$E|") || die $!;
			while(<IN>) {
				my @col=split /\t+/,$_;
				my ($chr,$pos,$CIAGR)=@col[2,3,5];
				next if ($CIAGR !~ /N/);
				#next if ($_!~/NH\:i\:1/);
				my $loc=$pos;
				while($CIAGR=~/(\d+)(\D+)/g) {
					my ($A,$B)=($1,$2);
					if ($B eq "N") {
						$loc=$pos+$A;
						if ( ($strand eq "+" && $pos<=$fst_ge && $pos>=$fst_gs && $A<=abs($sec_gs-$fst_ge)+3 && $A>=abs($sec_gs-$fst_ge)-3) ||
						($strand eq "-" && $loc<=$fst_ge && $loc>=$fst_gs && $A<=abs($sec_ge-$fst_gs)+3 && $A>=abs($sec_ge-$fst_gs)-3) ){
							$jun1++;
						} elsif ($A<abs($sec_ge-$fst_gs) &&
						( ($strand eq "+" && $pos<$sec_ts && $pos>=$fst_gs)||
						  ($strand eq "-" && $loc>$sec_te && $loc<=$fst_ge)) ) {
							$jun2++;
							if ($strand eq "+") {
								if ($pos>$fst_ge) {
									$alts=($alts==0 || $fst_ge+1<$alts)?$fst_ge:$alts;
									$alte=($alte==0 || $pos>$alte)?$pos:$alte;
								} else {
									$alts=($alts==0 || $pos+1<$alts)?$pos+1:$alts;
									$alte=($alte==0 || $fst_ge>$alte)?$fst_ge:$alte;
								}
							} else {
								if ($loc<=$fst_gs) {
									$alts=($alts==0 || $loc<$alts)?$loc:$alts;
									$alte=($alte==0 || $fst_gs-1>$alte)?$fst_gs-1:$alte;
								} else {
									$alts=($alts==0 || $fst_gs<$alts)?$fst_gs:$alts;
									$alte=($alte==0 || $loc-1>$alte)?$loc-1:$alte;
								}
							}
						}
					} elsif ($B=~ /[MDS]/) {
						$pos+=$A;
					}
				}
			}
			close IN;
			if ($jun1>$Cutoff/2 && $jun2>=$Cutoff/2) {
				my ($S,$E)=(sort{$a<=>$b}($fst_gs,$fst_ge,$sec_gs,$sec_ge))[0,-1];
				if (!exists $ASevent{A5SS}{"$chr:$alts-$alte"}) {
					$ASstat{A5SS}++;
					$ASevent{A5SS}{"$chr:$alts-$alte"}=1;
					$out.="$chr\tAS\tA5SS\t$S\t$E\t\.\t$strand\t\.\tID=A5SS$ASstat{A5SS}; $info; Junctions: 1^2=$jun1 1-^2=$jun2;\n";
					$out.="$chr\tAS\tconstitutive\t$fst_gs\t$fst_gs\t\.\t$strand\t\.\tID=A5SS$ASstat{A5SS};\n";
					$out.="$chr\tAS\talternative\t$alts\t$alte\t\.\t$strand\t\.\tID=A5SS$ASstat{A5SS};\n";
					$out.="$chr\tAS\tconstitutive\t$sec_gs\t$sec_ge\t\.\t$strand\t\.\tID=A5SS$ASstat{A5SS};\n";
				}
			}
		}
	}
	return $out;
}

sub detectA3SS {
	my ($bam,$chr,$strand,$info,$gexon,$texon)=@_;
	$info=~s/[\;\s]+$//;
	my $out="";
	if (@$gexon>=3 && @$texon>=3) {
		my ($lst_gs,$lst_ge)=@{$gexon->[-1]};
		my ($sec_gs,$sec_ge)=@{$gexon->[-2]};
		my ($lst_ts,$lst_te)=@{$texon->[-1]};
		my ($sec_ts,$sec_te)=@{$texon->[-2]};
		if ($lst_gs==$lst_ts && $lst_ge==$lst_te) {
			return "";
		} elsif ($lst_ts<=$lst_gs && $lst_te>=$lst_gs) {
			my ($jun1,$jun2)=(0,0);
			my ($S,$E)=(sort{$a<=>$b}($lst_gs,$lst_gs,$sec_gs,$sec_ge,$lst_ts,$lst_te,$sec_ts,$sec_te))[0,-1];
			my ($alts,$alte)=(0,0);
			open (IN,"samtools view $bam $chr:$S-$E|") || die $!;
			while(<IN>) {
				my @col=split /\t+/,$_;
				my ($chr,$pos,$CIAGR)=@col[2,3,5];
				next if ($CIAGR !~ /N/);
				my $loc=$pos;
				while($CIAGR=~/(\d+)(\D+)/g) {
					my ($A,$B)=($1,$2);
					if ($B eq "N") {
						$loc=$pos+$A;
						if ( ($strand eq "+" && $loc<=$lst_ge && $loc>=$lst_gs && $A<=abs($lst_gs-$sec_ge)+3 && $A>=abs($lst_gs-$sec_ge)-3) ||
						($strand eq "-" && $pos<$lst_ge && $pos>=$lst_gs && $A<=abs($sec_gs-$lst_ge)+3 && $A>=abs($sec_gs-$lst_ge)-3) ){
							$jun1++;
						} elsif ( ($A < abs($sec_gs-$lst_ge)-3) &&
						(($strand eq "+" && $loc<=$lst_ge && $loc>$sec_ge) ||
						 ($strand eq "-" && $pos>=$lst_gs && $pos<$sec_gs)) ) {
							$jun2++;
							if ($strand eq "+") {
								if ($loc<$lst_gs) {
									$alts=($alts==0 || $loc<$alts)?$loc:$alts;
									$alte=($alte==0 || $lst_gs-1>$alte)?$lst_gs-1:$alte;
								} else {
									$alts=($alts==0 || $lst_gs<$alts)?$lst_gs:$lst_gs;
									$alte=($alte==0 || $loc-1>$alte)?$loc-1:$alte;
								}
							} else {
								if ($pos>$lst_ge) {
									$alts=($alts==0 || $lst_ge+1<$alts)?$lst_ge+1:$alts;
									$alte=($alte==0 || $pos>$alte)?$pos:$alte;
								} else {
									$alts=($alts==0 || $pos+1<$alts)?$pos+1:$alts;
									$alte=($alte==0 || $lst_ge>$alte)?$lst_ge:$alte;
								}
							}
						}
					} elsif ($B=~/MDS/) {
						$pos+=$A;
					}
				}
			}
			close IN;
			if ($jun1>$Cutoff/2 && $jun2>=$Cutoff/2) {
				if (!exists $ASevent{A3SS}{"$chr:$alts-$alte"}) {
					my ($S,$E)=(sort{$a<=>$b}($lst_gs,$lst_ge,$sec_gs,$sec_ge))[0,-1];
					$ASstat{A3SS}++;
					$ASevent{A3SS}{"$chr:$alts-$alte"}=1;
					$out.="$chr\tAS\tA3SS\t$S\t$E\t\.\t$strand\t\.\tID=A3SS$ASstat{A3SS}; $info; Junctions: 1^2=$jun1 1-^2=$jun2;\n";
					$out.="$chr\tAS\tconstitutive\t$sec_gs\t$sec_ge\t\.\t$strand\t\.\tID=A3SS$ASstat{A3SS};\n";
					$out.="$chr\tAS\talternative\t$alts\t$alte\t\.\t$strand\t\.\tID=A3SS$ASstat{A3SS};\n";
					$out.="$chr\tAS\tconstitutive\t$lst_gs\t$lst_ge\t\.\t$strand\t\.\tID=A3SS$ASstat{A3SS};\n";
				}
			}
		}
	}
	return $out;
}

sub detectAFE {
	my ($bam,$chr,$strand,$info,$gexon,$texon)=@_;
	$info=~s/[\;\s]+$//;
	my $out="";
	if (@$gexon>3 && @$texon>3) {
		my ($fst_gs,$fst_ge)=@{$gexon->[0]};
		my ($sec_gs,$sec_ge)=@{$gexon->[1]};
		my ($trd_gs,$trd_ge)=@{$gexon->[2]};
		my ($jun1,$jun2,$jun3)=(0,0,0);
		my ($S,$E)=(sort{$a<=>$b}($fst_gs,$fst_ge,$trd_gs,$trd_ge))[0,-1];
		open (IN,"samtools view $bam $chr:$S-$E|") || die $!;
		while(<IN>) {
			my @col=split /\t+/,$_;
			my ($chr,$pos,$CIAGR)=@col[2,3,5];
			next if ($CIAGR !~ /N/);
			my $loc=$pos;
			while($CIAGR=~/(\d+)(\D+)/g) {
				my ($A,$B)=($1,$2);
				if ($B eq "N") {
					$loc=$pos+$A;
					if ( ($strand eq "+" && $pos<=$fst_ge && $loc>=$trd_gs) ||
					($strand eq "-" && $pos<=$trd_ge && $loc>=$fst_gs) ){
						$jun1++;
					} elsif( ($strand eq "+" && $pos<=$fst_ge && $loc>=$sec_gs && $loc<=$sec_ge) ||
					($strand eq "-" && $pos>=$sec_gs && $pos<=$sec_ge && $loc>=$fst_gs)) {
						$jun2++;
					}elsif ( ($strand eq "+" && $pos>=$sec_gs && $pos<=$sec_ge && $loc>=$trd_gs) ||
					($strand eq "-" && $pos>=$sec_gs && $pos<=$sec_ge && $loc>=$fst_gs) ){
						$jun3++;
					}
				} elsif ($B =~ /[MDS]/) {
					$pos+=$A;
				}
			}
		}
		close IN;
		if ($jun1>=$Cutoff/2 && $jun2+$jun3>=$Cutoff/2) {
			if (!exists $ASevent{AFE}{"$chr:$sec_gs-$sec_ge"}) {
				$ASstat{AFE}++;
				$ASevent{AFE}{"$chr:$sec_gs-$sec_ge"}=1;
				$out.="$chr\tAS\tAFE\t$S\t$E\t\.\t$strand\t\.\tID=AFE$ASstat{AFE}; $info; Junctions: 1\^3=$jun1 1\^2=$jun2 2\^3=$jun3;\n";
				$out.="$chr\tAS\talternative\t$fst_gs\t$fst_ge\t\.\t$strand\t\.\tID=AFE$ASstat{AFE};\n";
				$out.="$chr\tAS\talternative\t$sec_gs\t$sec_ge\t\.\t$strand\t\.\tID=AFE$ASstat{AFE};\n";
				$out.="$chr\tAS\tconstitutive\t$trd_gs\t$trd_ge\t\.\t$strand\t\.\tID=AFE$ASstat{AFE};\n";
			}
		}
	}
	return $out;
}

sub detectALE {
	my ($bam,$chr,$strand,$info,$gexon,$texon)=@_;
	$info=~s/[\;\s]+$//;
	my $out="";
	if (@$gexon>3 && @$texon>3) {
		my $lsti=@$gexon;
		my $seci=$lsti-1;
		my $trdi=$seci-1;
		my ($lst_gs,$lst_ge)=@{$gexon->[$lsti-1]};
		my ($sec_gs,$sec_ge)=@{$gexon->[$seci-1]};
		my ($trd_gs,$trd_ge)=@{$gexon->[$trdi-1]};
		my ($jun1,$jun2,$jun3)=(0,0,0);
		my ($S,$E)=(sort{$a<=>$b}($lst_gs,$lst_ge,$trd_gs,$trd_ge))[0,-1];
		open (IN,"samtools view $bam $chr:$S-$E|") || die $!;
		while(<IN>) {
			my @col=split /\t+/,$_;
			my ($chr,$pos,$CIAGR)=@col[2,3,5];
			next if ($CIAGR !~ /N/);
			my $loc=$pos;
			while($CIAGR=~/(\d+)(\D+)/g) {
				my ($A,$B)=($1,$2);
				if ($B eq "N") {
					$loc=$pos+$A;
					if ( ($strand eq "+" && $pos<=$trd_ge && $loc>=$lst_gs) ||
					($strand eq "-" && $pos<=$lst_ge && $loc>=$trd_gs) ){
						$jun1++;
					} elsif( ($strand eq "+" && $pos<=$trd_ge && $loc>=$sec_gs && $loc>=$sec_ge) ||
					($strand eq "-" && $pos>=$sec_gs && $pos<=$sec_ge && $loc>=$trd_gs)) {
						$jun2++;
					}elsif ( ($strand eq "+" && $pos>=$sec_gs && $pos<=$sec_ge && $loc>=$lst_gs) ||
					($strand eq "-" && $pos>=$lst_gs && $pos<=$lst_ge && $loc>=$sec_gs) ){
						$jun3++;
					}
				} elsif ($B=~/[MDS]/) {
					$pos+=$A;
				}
			}
		}
		close IN;
		if ($jun1>=$Cutoff/2 && $jun2+$jun3>=$Cutoff/2) {
			if (!exists $ASevent{ALE}{"$chr:$sec_gs-$sec_ge"}) {
				$ASstat{ALE}++;
				$ASevent{ALE}{"$chr:$sec_gs-$sec_ge"}=1;
				$out.="$chr\tAS\tALE\t$S\t$E\t\.\t$strand\t\.\tID=ALE$ASstat{ALE}; $info; Junctions: $trdi\^$lsti=$jun1 $trdi\^$seci=$jun2 $seci\^$lsti=$jun3;\n";
				$out.="$chr\tAS\tconstitutive\t$trd_gs\t$trd_ge\t\.\t$strand\t\.\tID=ALE$ASstat{ALE};\n";
				$out.="$chr\tAS\talternative\t$sec_gs\t$sec_ge\t\.\t$strand\t\.\tID=ALE$ASstat{ALE};\n";
				$out.="$chr\tAS\talternative\t$lst_gs\t$lst_ge\t\.\t$strand\t\.\tID=ALE$ASstat{ALE};\n";
			}
		}
	}
	return $out;
}

sub detectMXSE {
	my ($bam,$chr,$strand,$info,$gexon,$texon)=@_;
	$info=~s/[\;\s]+$//;
	my $out="";
	if (@$gexon>2 && @$texon>=2) {
		my ($overlap,$skip,$from,$to)=(0,0,-1,-1);
		my $J=0;
		for (my $i=0;$i<@$gexon;$i++) {
			my ($gs,$ge)=@{$gexon->[$i]};
			#my ($pre_gs,$pre_ge)=@{$gexon->[$i-1]} if ($i>0);
			#my ($aft_gs,$aft_ge)=@{$gexon->[$i+1]} if ($i<@$gexon-1);
			for (my $j=$J;$j<@$texon;$j++) {
				my ($ts,$te)=@{$texon->[$j]};
				#my ($pre_ts,$pre_te)=@{$texon->[$j-1]} if ($j>0);
				#my ($aft_ts,$aft_te)=@{$texon->[$j+1]} if ($j<@$texon-1);
				if (($strand eq "+" && $ts>$ge) || ($strand eq "-" && $te<$gs)) {
					$skip++ if ($from>=0);
					last;
				} elsif ($ts<$ge && $te>$gs) {
					$overlap=overlapSize([$gs,$ge],[$ts,$te]);
					if ($overlap>=($ge-$gs+1)/2 || $overlap>=($te-$ts+1)/2) {
						if ($skip<=0) {
							$from=$i;
						} else {
							$to=$i;
						}
					} else {
						$skip++ if ($from>=0);
					}
				}
				$J=$j-1 if ($j>0);
			}
			if ($from>=0 && $to>0 && $to>$from+1) {
				if ($skip>2 && $to>$from+2) {
					my $fromi=$from+1;
					my $toi=$to+1;
					my ($fromS,$fromE)=@{$gexon->[$from]};
					my ($toS,$toE)=@{$gexon->[$to]};
					my ($S,$E)=(sort{$a<=>$b}($fromS,$fromE,$toS,$toE))[0,-1];
					my %jun=();
					open (IN,"samtools view $bam $chr:$S-$E|") || die $!;
					while(<IN>) {
						my @col=split/\t+/,$_;
						my ($chr,$pos,$CIAGR)=@col[2,3,5];
						next if ($CIAGR !~ /N/);
						my $loc=$pos;
						while($CIAGR=~/(\d+)(\D+)/g) {
							my ($A,$B)=($1,$2);
							if ($B eq "N") {
								$loc=$pos+$A;
								for (my $m=$from;$m<$to;$m++) {
									for (my $n=$m+1;$n<=$to;$n++) {
										my ($mS,$mE)=@{$gexon->[$m]};
										my ($nS,$nE)=@{$gexon->[$n]};
										if ($strand eq "+") {
											if ($pos>=$mS-3 && $pos<=$mE+3 && $loc>=$nS-3 && $loc<=$nE+3) {
												$jun{$m+1}{$n+1}++;
											}
										} else {
											if ($loc>=$mS-3 && $loc<=$mE+3 && $pos>=$nS-3 && $pos<=$nE+3) {
												$jun{$m+1}{$n+1}++;
											}
										}
									}
								}
							} elsif ($B =~ /[MDS]/) {
								$pos+=$A;
							}
						}
					}
					close IN;
					my $total=0;
					foreach my $m(keys %jun) {
						foreach my $n(keys %{$jun{$m}})
						{
							$total+=$jun{$m}{$n} if (abs($m-$n)>1);
						}
					}
					if (exists $jun{$fromi}{$toi} && $jun{$fromi}{$toi}>=$Cutoff && $total>=$Cutoff && !exists $ASevent{MXE}{"$chr:$S-$E"}) {
						$ASstat{MXE}++;
						$ASevent{MXE}{"$chr:$S-$E"}=1;
						$out.="$chr\tAS\tMXE\t$S\t$E\t\.\t$strand\t\.\tID=MXE$ASstat{MXE}; Junctions:";
						foreach my $m (sort{$a<=>$b}keys %jun) {
							foreach my $n (sort{$a<=>$b}keys %{$jun{$m}}) {
								$out.=" $m\^$n=$jun{$m}{$n}";
							}
						}
						$out.=";\n";
						$out.="$chr\tAS\talternative\t$fromS\t$fromE\t\.\t$strand\t\.\tID=MXE$ASstat{MXE};\n";
						for (my $x=$from+1;$x<$to;$x++)
						{
							my ($skipS,$skipE)=@{$gexon->[$x]};
							$out.="$chr\tAS\talternative\t$skipS\t$skipE\t\.\t$strand\t\.\tID=MXE$ASstat{MXE};\n";
							$ASevent{MXE}{"$chr:$skipS-$skipE"}=1;
						}
						$out.="$chr\tAS\tconstitutive\t$toE\t$toE\t\.\t$strand\t\.\tID=MXE$ASstat{MXE};\n";
					}
				} elsif ($skip>0) {
					my $fromi=$from+1;
					my $skipi=$fromi+1;
					my $toi=$to+1;
					my ($fromS,$fromE)=@{$gexon->[$from]};
					my ($skipS,$skipE)=@{$gexon->[$fromi]};
					my ($toS,$toE)=@{$gexon->[$to]};
					my ($S,$E)=(sort{$a<=>$b}($fromS,$fromE,$toS,$toE))[0,-1];
					my ($jun1,$jun2,$jun3)=(0,0,0);
					open (IN,"samtools view $bam $chr:$S-$E|") || die $!;
					while(<IN>) {
						my @col=split/\t+/,$_;
						my ($chr,$pos,$CIAGR)=@col[2,3,5];
						next if ($CIAGR !~ /N/);
						my $loc=$pos;
						while($CIAGR=~/(\d+)(\D+)/g) {
							my ($A,$B)=($1,$2);
							if ($B eq "N") {
								$loc=$pos+$A;
								if ( ($strand eq "+" && $pos<=$fromE+3 && $loc>=$toS-3) ||
								($strand eq "-" && $pos<=$toE+3 && $loc>=$fromS-3) ){
									$jun1++;
								} elsif ( ($strand eq "+" && $pos<=$fromE+3 && $loc>=$skipS-3 && $loc>=$skipE+3) ||
								($strand eq "-" && $loc>=$fromS-3 && $pos>=$skipS-3 && $pos>=$skipE+3) ){
									$jun2++;
								} elsif ( ($strand eq "+" && $pos>=$skipS-3 && $pos<=$skipE+3 && $loc>=$toS-3) ||
								($strand eq "-" && $loc>=$skipS-3 && $loc<=$skipE+3 && $pos<=$toE+3) ){
									$jun3++;
								}
							} elsif ($B=~/[MDS]/) {
								$pos+=$A;
							}
						}
					}
					close IN;
					if ($jun1>=$Cutoff && !exists $ASevent{SE}{"$chr:$skipS-$skipE"} && !exists $ASevent{MXE}{"$chr:$skipS-$skipE"}) {
						$ASstat{SE}++;
						$ASevent{SE}{"$chr:$S-$E"}=1;
						$out.="$chr\tAS\tSE\t$S\t$E\t\.\t$strand\t\.\tID=SE$ASstat{SE}; $info; Junctions: $fromi\^$toi=$jun1 $fromi\^$skipi=$jun2 $skipi\^$toi=$jun3;\n";
						$out.="$chr\tAS\tconstitutive\t$fromS\t$fromE\t\.\t$strand\t\.\tID=SE$ASstat{SE};\n";
						$out.="$chr\tAS\talternative\t$skipS\t$skipE\t\.\t$strand\t\.\tID=SE$ASstat{SE};\n";
						$out.="$chr\tAS\tconstitutive\t$toE\t$toE\t\.\t$strand\t\.\tID=SE$ASstat{SE};\n";
					}
				}
				($overlap,$skip,$from,$to)=(0,0,-1,-1);
			}
		}
	}
	return $out;
}


sub detectRI {
	my ($bam,$chr,$strand,$info,$gexon,$texon)=@_;
	$info=~s/[\;\s]+$//;
	my $out="";
	if (@$gexon>2 && @$texon>=2) {
		my ($overlap,$fromi,$toi)=(0,-1,-1);
		for (my $i=1;$i<@$gexon;$i++) {
			$fromi=$i-1;
			$toi=$i;
			my ($gs,$ge)=@{$gexon->[$i]};
			my ($pre_gs,$pre_ge)=@{$gexon->[$i-1]};
			my ($A,$B,$C,$D)=sort{$a<=>$b}($gs,$ge,$pre_gs,$pre_ge);
			for (my $j=0;$j<@$texon;$j++) {
				my ($ts,$te)=@{$texon->[$j]};
				my ($E,$F)=();
				if ($j>0)
				{
					my ($pre_ts,$pre_te)=@{$texon->[$j-1]};
					($E,$F)=sort{$a<=>$b}($ts,$te,$pre_ts,$pre_te)[1,2];
				} else {
					($E,$F)=sort{$a<=>$b}($ts,$te,$pre_gs,$pre_ge)[1,2];
				}
				my ($len1,$len2)=sort{$a<=>$b}($F-$E,$C-$B);
				if ($ts<=$C && $te>=$B && (($ts!=$gs && $te!=$ge)||($ts!=$pre_gs && $te!=$pre_ge)) ) {
					$overlap=overlapSize([$B,$C],[$ts,$te]);
					if ($overlap>0) {
						my ($jun1,$jun2)=(0,0);
						my ($alts,$alte)=(0,0);
						open (IN,"samtools view $bam $chr:$A-$D|") || die $!;
						while(<IN>) {
							my @col=split/\t+/,$_;
							my ($chr,$pos,$CIAGR)=@col[2,3,5];
							next if ($CIAGR !~ /N/);
							next if ($_ !~ /NH:i:1/ || $col[6] ne "=");
							my $loc=$pos;
							while($CIAGR=~/(\d+)(\D+)/g) {
								my ($len,$mark)=($1,$2);
								if ($mark eq "N") {
									$loc=$pos+$len;
									if ( $pos>=$A && $loc<=$D ) {
										if ($len<=$len2-3 && $len>=$len1-3) {
											$jun1++;
										} elsif ($len<=$len1 && $len1<$len2) {
											$jun2++;
											$alts=($pos>$B+1)?$pos:$B+1;
											$alte=($loc<$C-1)?$loc:$C-1;
										}
									}
								} elsif ($mark =~ /[MDS]/) {
									$pos+=$len;
								}
							}
						}
						close IN;
						if ($jun2>=$Cutoff && $alte>$alts) {
							if (!exists $ASevent{RI}{"$chr:$A-$D"}) {
								$ASstat{RI}++;
								$ASevent{RI}{"$chr:$A-$D"}=1;
								my $from=$fromi+1;
								my $to=$toi+1;
								$out.="$chr\tAS\tRI\t$A\t$D\t\.\t$strand\t\.\tID=RI$ASstat{RI}; $info; Junctions: $from\^$to=$jun1 $from\-$to=$jun2;\n";
								$out.="$chr\tAS\tconstitutive\t$pre_gs\t$pre_ge\t\.\t$strand\t\.\tID=RI$ASstat{RI};\n";
								$out.="$chr\tAS\talternative\t$alts\t$alte\t\.\t$strand\t\.\tID=RI$ASstat{RI};\n";
								$out.="$chr\tAS\tconstitutive\t$gs\t$ge\t\.\t$strand\t\.\tID=RI$ASstat{RI};\n";
							}
						}
					}
				}
			}
			
		}
	}
	return $out;
}

sub checkIndex {
	my $file = shift;
	if ($file=~/bam$/i) {
		if (-f "$file.bai") {
			return 1;
		} else {
			system "samtools index $file";
		}
	} elsif ($file=~/fa$/i) {
		if (-f "$file.fai") {
			return 1;
		} else {
			system "samtools faidx $file";
		}
	}
}

sub overlapSize
{
	my ($block1_p,$block2_p) = @_;
	my ($overlap_start,$overlap_end) = (sort {$a<=>$b} ($block1_p->[0],$block2_p->[0],$block1_p->[1],$block2_p->[1]))[1,2];
	return ($overlap_end - $overlap_start + 1);
}

__END__