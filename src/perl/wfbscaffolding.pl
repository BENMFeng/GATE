#!/usr/bin/perl -w
##########################################################################
#  Copyright (c) 2012 - BENM(Binxiao) Feng                               #
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
=head1 Name

wfbscaffolding -- WGS, FPC, BACends Scaffloding program.

=head1 Version

Author: BENM <binxiaofeng@gmail.com>
Date: Dec 26th, 2012, v0.1.2

=head1 Usage
  --wgs <infile>	Input WGS assembly fasta file
  --raw <infile>	Input WGS raw assembly fasta file if --wgs is RM.out file
  --repeatmask <infile>	Input RepeatMasker annoation gff file
  --fpc <infile>	Input FPC file
  --bac <infile>	Input BAC ends file
  --group <outfile>	Output FPC CTG group info file
  --pem <outfile>	Output BAEends mapping pairs file's prefix
  --roadmap <outfile>	Output Roadmap file
  --scaf <outfile>	Output Scaffolds result file
  --stat <outfile>	Output Stat of scafflods length report
  --formatdb <exc>	Input Excuted formatdb, default: /usr/local/bin/formatdb
  --blastall <exc>	Input Excuted blastall, default: /usr/local/bin/blastall
  --e_value <str>	Set E-value for blast, defalut: 1e-10
  --identity <float>	Set Identity value for blast, default: 90
  --aln_len <int>	Set LCS/HSP Length threshold, default: 50
  --overlap <float> 	Set Overlap rate of the minimum alignment blocks, default: 0.6
  --edge_gap <int>	Set maxmium edge gap size, default: 5
  --min <intk>		Set minimum insert size of BACEnds: 10k
  --max <intk>		Set maxmium insert size of BACEnds: 150k
  --verbose		print STDERR running progress information to screen
  --help		show this help
  

=head1 Example

  $ perl wfbscaffolding.pl --wgs MW_v2.fasta --fpc cassava2012.fpc --bac 20121112_BES.fa \
  --group cassava2012.fpc.group.xls --pem 20121112BES --roadmap CTG_SCF.txt --scaf superScaf.fa \
  --stat ctg_group.stat.txt --identity 97 --aln_len 100 --overlap 0.85 --min 40k --max 230k \
  --overlap 0.2 --edge_gap 5 --verbose > error.log 2>&1 &

=cut

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my %opts;
my ($FPC_file,$BAC_file,$WGS_file,$WGS_raw,$RepeatMasker_gff,$CTG_group,$PEM_file,$Roadmap_file,$Scaf_fa,$Stat_report,$formatdb,$blastall,$E_value,$Identity,$Aln_len,$Overlap,$Edge_gap,$Min,$Max,$Verbose,$Help);
GetOptions(
	\%opts,
	"wgs:s"=>\$WGS_file,
	"raw:s"=>\$WGS_raw,
	"repeatmask:s"=>\$RepeatMasker_gff,
	"fpc:s"=>\$FPC_file,
	"bac:s"=>\$BAC_file,
	"group:s"=>\$CTG_group,
	"pem:s"=>\$PEM_file,
	"roadmap:s"=>\$Roadmap_file,
	"scaf:s"=>\$Scaf_fa,
	"stat:s"=>\$Stat_report,
	"formatdb:s"=>\$formatdb,
	"blastall:s"=>\$blastall,
	"e_value:s"=>\$E_value,
	"identity:s"=>\$Identity,
	"aln_len:i"=>\$Aln_len,
	"overlap:s"=>\$Overlap,
	"edge_gap:i"=>\$Edge_gap,
	"min:s"=>\$Min,
	"max:s"=>\$Max,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if ( (!defined $WGS_file) || (!defined $FPC_file) || (!defined $BAC_file) || (!defined $PEM_file) || ($Help) );

## default setting
$formatdb ||= "/usr/local/bin/formatdb";
$blastall ||= "/usr/local/bin/blastall";
$E_value ||= "1e-10";
$Identity ||= 90;
$Aln_len ||= 50;
$Overlap ||= 0.6;
$Edge_gap ||= 5;
$Min = ((defined $Min)&&($Min=~m/(\d+)k/)) ? $1*1000 : 10000;
$Max = ((defined $Max)&&($Max=~m/(\d+)k/)) ? $1*1000 : 150000;
$WGS_raw ||=$WGS_file;

## main
my %Seq;
if ((defined $Roadmap_file)&&(-f $Roadmap_file))
{
	print STDERR "$Roadmap_file is existing! It will be overwrote!\n";
	system "rm $Roadmap_file\n";
}

my $c=check_name($BAC_file,$WGS_file);
if (($c==0)||($c==2))
{
	if (-f $BAC_file && !-f "$BAC_file.addleninfo")
	{
		rename_fa($BAC_file,"$BAC_file.addleninfo");
		$BAC_file="$BAC_file.addleninfo";
	}
	elsif (-f "$BAC_file.addleninfo")
	{
		$BAC_file="$BAC_file.addleninfo";
	}
	if (-f $FPC_file && !-f "$WGS_file.addfpcinfo")
	{
		fpc_name($FPC_file,$BAC_file);
		$BAC_file="$BAC_file.addfpcinfo";
	}
	elsif (-f "$WGS_file.addfpcinfo")
	{
		$BAC_file="$BAC_file.addfpcinfo";
	}
}
if ($c<2)
{
	if (defined $RepeatMasker_gff && -f $RepeatMasker_gff)
	{
		repeat_mask($RepeatMasker_gff,\$WGS_file,\%Seq);
	}
	rename_fa($WGS_file,"$WGS_file.addleninfo") if (-f $WGS_file && !-f "$WGS_file.addleninfo");
	$WGS_file="$WGS_file.addleninfo";
}
my @query_path=split /\//,$BAC_file;
my @subject_path=split /\//,$WGS_file;
my $blast_result=(@ARGV>0)?$ARGV[0]:"$query_path[-1]\_blastn_$subject_path[-1].m8";
blast($BAC_file,$WGS_file,$blast_result) unless (-f $blast_result);
if ($blast_result!~/besthit$/)
{
	get_best_blast($blast_result);
	$blast_result.=".besthit";
}

my %Scf_map=();
my %Wgs_map=();
blast2pe($blast_result,\%Scf_map,\%Wgs_map) if (-f $blast_result);
traceroad(\%Scf_map,\%Wgs_map);
outputScaf($Roadmap_file,$WGS_raw,\%Seq) if (defined $Roadmap_file && -s $Roadmap_file);
print STDERR (&get_time,": done!\n") if ($Verbose);

sub repeat_mask
{
	my ($gff,$wgs,$contig_seq) = @_;
	my %hash;
	print STDERR (&get_time,": parsing RepeatMasker anntation file: $gff\n") if ($Verbose);
	open (IN,$gff) || die $!;
	while(<IN>)
	{
		my @t=split /\t+/,$_;
		push @{$hash{$t[0]}},[$t[3],$t[4]];
	}
	close IN;
	open (FA,$$wgs) || die $!;
	print STDERR (&get_time,": hard mask for $$wgs: $$wgs.repeatmasker.out\n") if ($Verbose);
	open (OUT,">$$wgs.repeatmasker.out") || die $!;
	my ($name,$seq)=("","");
	while(<FA>)
	{
		if (/^\>(\S+)/ || eof)
		{
			my $id=$1 if (/^\>(\S+)/);
			if (!/^\>/ && eof && /[ACGTN]/i)
			{
				chomp;
				$seq.=$_;
			}
			if (length($seq)>0)
			{
				$$contig_seq{$name}=$seq if (defined $contig_seq);
				if (exists $hash{$name})
				{
					foreach my $frag(@{$hash{$name}})
					{
						if ($frag->[1]-1>length($seq))
						{
							print STDERR ("ERROR:\t$name\t",$frag->[0]-1,"\t",$frag->[1]-$frag->[0]+1,"\t",length($seq),"\n") if ($Verbose);
							substr($seq,$frag->[0]-1,length($seq)-$frag->[0]+1,"N"x(length($seq)-$frag->[0]+1)) if ($frag->[0]<length($seq));
							next;
						}
						substr($seq,$frag->[0]-1,$frag->[1]-$frag->[0]+1,"N"x($frag->[1]-$frag->[0]+1));
					}
				}
				Display(\$seq);
				print OUT ">$name\n";
				print OUT $seq;
			}
			$name=$id;
			$seq="";
		}
		else
		{
			chomp;
			$seq.=$_;
		}
	}
	close FA;
	close OUT;
	$$wgs = "$$wgs.repeatmasker.out";
}

sub blast
{
	my ($query,$subject,$out)=@_;
	print STDERR (&get_time,": $formatdb -i $subject -p F -o T\n") if ($Verbose);
	system "$formatdb -i $subject -p F -o T" unless (-f "$subject.nsq");
	print STDERR (&get_time,": $blastall -p blastn -i $query -d $subject -e $E_value -m8 -F F -W $Aln_len -o $out -a 12\n") if ($Verbose);
	system "$blastall -p blastn -i $query -d $subject -e $E_value -m8 -F F -W $Aln_len -o $out -a 12";
}

sub get_best_blast
{
	my $blast=shift;
	my $block=shift;
	my ($p,$q,$m,$n,$x,$y)=((defined $block)&&($block =~ m/subject/i)) ? (1,0,8,9,6,7) : (0,1,6,7,8,9);
	my $output=$blast.".besthit";
	my ($name,@info)=("",());
	open (IN,$blast) || die $!;
	print STDERR (&get_time,": searching for best blast hit: $blast\n") if ($Verbose);
	open (OUT,">$output") || die $!;
	while(<IN>)
	{
		chomp;
		my @t=split /\t+/,$_; #$QryId, $SubId, $Identity, $AlignLen, $MisMatch, $GapOpen, $QStart,$QEnd, $SStart, $SEnd, $EVal, $BitScore
		next if (($t[2]<$Identity)||($t[10]>$E_value)||($t[3]-$t[4]-$t[5]<$Aln_len));
		if ((($name ne "")&&($name ne $t[$p]))||(eof))
		{
			if (eof)
			{
				if ($name eq $t[$p])
				{
					if ($t[$m]>$t[$n])
					{
						($t[$m],$t[$n])=($t[$n],$t[$m]);
						($t[$x],$t[$y])=($t[$y],$t[$x]);
					}
					push @info,[@t];
				}
				else
				{
					print OUT "$_\n";
				}
			}
			my $out = "";
			if (@info>1)
			{
				@info = sort{$a->[$m]<=>$b->[$m]} @info;
				my @pre = @{$info[0]};
				my $OverlapSize=0;
				for (my $i=1;$i<@info;$i++)
				{
					my ($preM,$preN)=($pre[$m]<$pre[$n])?($pre[$m],$pre[$n]):($pre[$n],$pre[$m]);
					my ($preX,$preY)=($pre[$x]<$pre[$y])?($pre[$x],$pre[$y]):($pre[$y],$pre[$x]);
					my ($infoM,$infoN)=($info[$i][$m]<$info[$i][$n])?($info[$i][$m],$info[$i][$n]):($info[$i][$n],$info[$i][$m]);
					my ($infoX,$infoY)=($info[$i][$x]<$info[$i][$y])?($info[$i][$x],$info[$i][$y]):($info[$i][$y],$info[$i][$x]);
					$OverlapSize=get_overlap([$preM,$preN],[$infoM,$infoN]);
					my $alen=$1 if ($pre[$p]=~/length\=(\d+)/);
					my $blen=$1 if ($pre[$q]=~/length\=(\d+)/);
					my $aMatchRate=($pre[3]-$pre[4]-$pre[5])/$alen;
					my $bMatchRate=($pre[3]-$pre[4]-$pre[5])/$blen;
					my $preMatch=$pre[3]-$pre[4]-$pre[5];
					my $tmpMatch=$info[$i][3]-$info[$i][4]-$info[$i][5];
					if ($OverlapSize<=0)
					{
						if ( ( (($preM<=$Edge_gap)||($alen-$preN<=$Edge_gap)) && (($preX<=$Edge_gap)||($blen-$preY<=$Edge_gap)) )
						    || ($aMatchRate>=$Overlap) || ($bMatchRate>=$Overlap) )
						{
							$out.=(join "\t",@pre)."\n";
						}
					}
					else
					{
						if ( ($preMatch<$tmpMatch) || (($preMatch==$tmpMatch)&&(($info[$i][10]<$pre[10])||($info[$i][11]>$pre[11]))) )
						{
							$out=(join "\t",@{$info[$i]})."\n";
						}
						if (($preMatch==$tmpMatch)&&($info[$i][10]==$pre[10])&&($info[$i][11]==$pre[11]))
						{
							$out.=(join "\t",@{$info[$i]})."\n";
						}
					}
					if ( ($preMatch<$tmpMatch) || (($preMatch==$tmpMatch)&&(($info[$i][10]<$pre[10])||($info[$i][11]>$pre[11]))) )
					{
						@pre=@{$info[$i]};
					}
				}
				if (($OverlapSize<=0)||(@pre>0))
				{
					my $alen=$1 if ($pre[$p]=~/length\=(\d+)/);
					my $blen=$1 if ($pre[$q]=~/length\=(\d+)/);
					my ($infoM,$infoN)=($pre[$m]<$pre[$n])?($pre[$m],$pre[$n]):($pre[$n],$pre[$m]);
					my ($infoX,$infoY)=($pre[$x]<$pre[$y])?($pre[$x],$pre[$y]):($pre[$y],$pre[$x]);
					my $aMatch=($pre[3]-$pre[4]-$pre[5])/$alen;
					my $bMatch=($pre[3]-$pre[4]-$pre[5])/$blen;
					if ( ((($infoM<=$Edge_gap)||($alen-$infoN<=$Edge_gap))&&(($infoX<=$Edge_gap)||($blen-$infoY<=$Edge_gap)))||($aMatch>=$Overlap)||($bMatch>=$Overlap) )
					{
						$out.=(join "\t",@pre)."\n";
					}
				}
			}
			else
			{
				$out=(join "\t",@{$info[0]})."\n" if (@info>0);
			}
			print OUT $out;
			@info=();
		}
		if ($t[$m]>$t[$n])
		{
			($t[$m],$t[$n])=($t[$n],$t[$m]);
			($t[$x],$t[$y])=($t[$y],$t[$x]);
		}
		push @info,[@t];
		$name=$t[$p];
	}
	close IN;
	close OUT;
}

sub get_overlap
{
	my ($block1_p,$block2_p) = @_;
	my $combine_start = ($block1_p->[0] < $block2_p->[0]) ?  $block1_p->[0] : $block2_p->[0];
	my $combine_end   = ($block1_p->[1] > $block2_p->[1]) ?  $block1_p->[1] : $block2_p->[1];
	my $overlap_size = ($block1_p->[1]-$block1_p->[0]+1) + ($block2_p->[1]-$block2_p->[0]+1) - ($combine_end-$combine_start+1);
	return $overlap_size;
}

sub blast2pe
{
	my ($blast_result,$SCF,$MAP)=@_;
	open (IN,$blast_result) || die "failed to read $blast_result\n";
	print STDERR (&get_time,": reading blast.m8: $blast_result\n") if ($Verbose);
	my %BAC=();
	while(<IN>)
	{
		chomp;
		my @t=split /\t+/,$_;
		next if (($t[2]<$Identity)||($t[3]-$t[4]-$t[5]<$Aln_len));
		if ($t[0]=~/(\S+)\_([F|R])\_\S+/)
		{
			my ($clone,$p1)=($1,$2);
			my ($qS,$qE)=($t[6],$t[7]);
			my ($ss,$sS,$sE)=($t[8]<$t[9])?("+",$t[8],$t[9]):("-",$t[9],$t[8]);
			my $q_len=$1 if ($t[0]=~/length\=(\d+)/i);
			my $s_len=$1 if ($t[1]=~/length\=(\d+)/i);
			my $edge=0;
			my $ctginfo="";
			if ($t[0]=~/\_(ctg[^\_]+)\_(\d+)\-(ctg[^\_]+)\_(\d+)/)
			{
				next if ($1 ne $3);
				$ctginfo="$1:$2-$4";
				my @A=($t[0],$ss,$q_len,$t[6],$t[7],$t[1],$s_len,$sS,$sE,$ctginfo);
				push @{$BAC{$clone}{$p1}},[@A];
			}
		}
	}
	close IN;
	open (PE,">$PEM_file.pair") || die "failed to write to $PEM_file.pair\n";
	open (SE,">$PEM_file.single") || die "failed to write to $PEM_file.single\n";
	print STDERR (&get_time,": getting blast_pe...\n") if ($Verbose);
	foreach my $clone(keys %BAC)
	{
		if (exists $BAC{$clone}{"F"} && exists $BAC{$clone}{"R"})
		{
			my $p1="F";
			for (my $i=0;$i<@{$BAC{$clone}{$p1}};$i++)
			{
				my @A=@{${$BAC{$clone}{$p1}}[$i]};
				my $p2 = ($p1 eq "F") ? "R" : "F";
				for (my $j=0;$j<@{$BAC{$clone}{$p2}};$j++)
				{
					my @B=@{${$BAC{$clone}{$p2}}[$j]};
					my ($ctgA,$ctgAS,$ctgAE)=split /[\:\-]/,$A[-1];
					my ($ctgB,$ctgBS,$ctgBE)=split /[\:\-]/,$B[-1];
					next if ($ctgA ne $ctgB);
					my ($Dmin,$Dmax)=(0,0);
					if ($A[5] eq $B[5])
					{
						($Dmin,$Dmax) = mm(abs($A[7]-$B[7]),abs($A[7]-$B[8]),abs($A[8]-$B[7]),abs($A[8]-$B[8]));
						next if (($Dmin<$Min)||($Dmax>$Max));
						next if ($p1 eq $p2);
						if (!exists $$MAP{$A[5]}{$ctgA} || $$MAP{$A[5]}{$ctgA} ne "$clone:$p1")
						{
							push @{$SCF->{$ctgA}},[$ctgAS,$ctgAE,$clone,$p1,@A[5,1,7,8]];
							$$MAP{$A[5]}{$ctgA}="$clone:$p1";
						}
						if (!exists $$MAP{$B[5]}{$ctgB} || $$MAP{$B[5]}{$ctgB} ne "$clone:$p2")
						{
							push @{$SCF->{$ctgB}},[$ctgBS,$ctgBE,$clone,$p2,@B[5,1,7,8]] unless ($ctgBE != $ctgAE);
							$$MAP{$B[5]}{$ctgB}="$clone:$p2";
						} 
					}
					else
					{
						($Dmin,$Dmax) = mm(($A[7]+$B[7]),($A[7]+$B[6]-$B[8]),($A[6]-$A[8]+$B[7]),($A[6]-$A[8]+$B[6]-$B[8]));
						next if ($Dmax>$Max);
						if (!exists $$MAP{$A[5]}{$ctgA} || $$MAP{$A[5]}{$ctgA} ne "$clone:$p1")
						{
							push @{$SCF->{$ctgA}},[$ctgAS,$ctgAE,$clone,$p1,@A[5,1,7,8]];
							$$MAP{$A[5]}{$ctgA}="$clone:$p1";
						}
						if (!exists $$MAP{$B[5]}{$ctgB} || $$MAP{$B[5]}{$ctgB} ne "$clone:$p2")
						{
							push @{$SCF->{$ctgB}},[$ctgBS,$ctgBE,$clone,$p2,@B[5,1,7,8]];
							$$MAP{$B[5]}{$ctgB}="$clone:$p2";
						}
					}
					print PE ((join "\t",(@A,@B)),"\n");
				}
			}
		}
		else
		{
			foreach my $p("F","R")
			{
				if (exists $BAC{$clone}{$p})
				{
					for (my $i=0;$i<@{$BAC{$clone}{$p}};$i++)
					{
						print SE ((join "\t",@{${$BAC{$clone}{$p}}[$i]}),"\n");
					}
				}
			}
		}
	}
	close PE;
	close SE;
}

sub traceroad
{
	my ($SCF,$MAP)=@_;
	my $i=1;
	open (RD,">$Roadmap_file") if (defined $Roadmap_file);
	print STDERR (&get_time,": traceing roadmap\n") if ($Verbose);
	my $out="";
	foreach my $ctg(keys %$SCF)
	{
		print RD "#$ctg\n";
		my %info=();
		my %wgs_count=();
		my %wgs_bac=();
		my %bac_wgs=();
		my $i=0;
		my $pre_bac="";
		foreach my $map(sort{$a->[1]<=>$b->[1]}@{$SCF->{$ctg}})
		{
			next if (keys %{$$MAP{$map->[4]}}>1);
			$i++ if ($map->[2] ne $pre_bac);
			$wgs_bac{"$map->[4]:$map->[3]"}{$i}="$map->[3]:$map->[5]:$map->[6]:$map->[7]";
			$bac_wgs{$i}{"$map->[4]:$map->[3]"}=1;
			$wgs_count{$map->[4]}++;
			push @{$info{$i}},[@$map];
			$pre_bac=$map->[2];
		}
		my %roadmap;
		my %map;
		my %RC;
		foreach my $j(sort{$a<=>$b}keys %bac_wgs)
		{
			foreach my $wgs(keys %{$bac_wgs{$j}})
			{
				my $wgs_name=$1 if ($wgs=~/([^\:]+)\:/);
				my ($P,$strand,$Start,$End)=split /\:/,$wgs_bac{$wgs}{$j};
				if ($wgs_count{$wgs_name}>1)
				{
					if (!exists $map{$wgs_name})
					{
						$roadmap{$j}{3}="$wgs:$strand";
					}
					else
					{
						$roadmap{$j}{5}="$wgs:$strand";
					}
					$roadmap{$j}{$P}="$wgs:$strand";
					$map{$wgs_name}=$j;
				}
				else
				{
					if (!exists $roadmap{$j}{$P})
					{
						if (exists $roadmap{$j}{5} && !exists $roadmap{$j}{3})
						{
							$roadmap{$j}{3}="$wgs:$strand";
						}
						if (exists $roadmap{$j}{3} && !exists $roadmap{$j}{5})
						{
							$roadmap{$j}{5}="$wgs:$strand";
						}
						$roadmap{$j}{$P}="$wgs:$strand";
					}
					else
					{
						my $pre_wgs=$1 if ($roadmap{$j}{$P}=~/(\S+)\:[+-]/);
						my $pre_wgs_name=$1 if ($pre_wgs=~/([^\:]+)\:/);
						next if ($wgs_count{$pre_wgs_name}>1);
						my ($preStart,$preEnd)=(split /\:/,$wgs_bac{$pre_wgs}{$j})[2,3];
						if (abs($preEnd-$preStart)<abs($End-$Start))
						{
							if (exists $roadmap{$j}{5} && !exists $roadmap{$j}{3})
							{
								$roadmap{$j}{3}="$wgs:$strand";
							}
							if (exists $roadmap{$j}{3} && !exists $roadmap{$j}{5})
							{
								$roadmap{$j}{5}="$wgs:$strand";
							}
							$roadmap{$j}{$P}="$wgs:$strand";
						}
					}
				}
			}
			if (exists $roadmap{$j}{5} && !exists $roadmap{$j}{3})
			{
				my ($wgs5,$P5,$strand5)=split /\:/,$roadmap{$j}{5};
				my $P3=($P5 eq "F") ? "R" : "F";
				$roadmap{$j}{3}=$roadmap{$j}{$P3} if (exists $roadmap{$j}{$P3} && defined $roadmap{$j}{$P3});
			}
			if (!exists $roadmap{$j}{5} && exists $roadmap{$j}{3})
			{
				my ($wgs3,$P3,$strand3)=split /\:/,$roadmap{$j}{3};
				my $P5=($P3 eq "F") ? "R" : "F";
				$roadmap{$j}{5}=$roadmap{$j}{$P5} if (exists $roadmap{$j}{$P5} && defined $roadmap{$j}{$P5});
			}
			if (exists $roadmap{$j}{5} && exists $roadmap{$j}{3})
			{
				my ($wgs5,$p5,$strand5)=split /\:/,$roadmap{$j}{5};
				my ($wgs3,$p3,$strand3)=split /\:/,$roadmap{$j}{3};
				if (exists $roadmap{$j-1}{3})
				{
					my ($pre_wgs,$pre_p,$pre_strand)=split /\:/,$roadmap{$j-1}{3};
					if ($pre_wgs eq $wgs5)
					{
						my $rc=judge_RC($pre_p,$pre_strand,$p5,$strand5);
						if ($rc==2 || (!exists $RC{$j-1}{3} && $rc==1) || (exists $RC{$j-1}{3} && $rc==0) )
						{
							$RC{$j}{5}=1;
						}
						elsif ($rc==3 && $wgs3 ne $wgs5)
						{
							my $tmp=$roadmap{$j}{5};
							$roadmap{$j}{5}=$roadmap{$j-1}{3};
							$roadmap{$j-1}{3}=$tmp;
						}
					}
				}
				elsif (exists $roadmap{$j+1}{5})
				{
					my ($aft_wgs,$aft_p,$aft_strand)=split /\:/,$roadmap{$j+1}{5};
					if ($aft_wgs eq $wgs3)
					{
						my $rc=judge_RC($p3,$strand3,$aft_p,$aft_strand);
						if ($rc==1)
						{
							$RC{$j}{3}=1;
						}
						elsif ($rc==3 && $wgs3 ne $wgs5)
						{
							my $tmp=$roadmap{$j}{3};
							$roadmap{$j}{3}=$roadmap{$j+1}{5};
							$roadmap{$j+1}{5}=$tmp;
						}
					}
				}
				my $rc=judge_RC($p5,$strand5,$p3,$strand3);
				if (!exists $RC{$j}{5} && !exists $RC{$j}{3})
				{
					if ($rc==1)
					{
						$RC{$j}{5}=1;
					}
					elsif ($rc==2)
					{
						$RC{$j}{3}=1;
					}
					elsif ($rc==3)
					{
						my $tmp=$roadmap{$j}{3};
						$roadmap{$j}{3}=$roadmap{$j}{5};
						$roadmap{$j}{5}=$tmp;
					}
				}
				($wgs5,$p5,$strand5)=split /\:/,$roadmap{$j}{5};
				($wgs3,$p3,$strand3)=split /\:/,$roadmap{$j}{3};
				if ($wgs5 eq $wgs3)
				{
					if (exists $RC{$j}{5} && !exists $RC{$j}{3})
					{
						$RC{$j}{3}=1;
					}
					elsif (!exists $RC{$j}{5} && exists $RC{$j}{3})
					{
						$RC{$j}{5}=1;
					}
				}
				for my $flag(5,3)
				{
					my ($wgs,$p,$strand)=split /\:/,$roadmap{$j}{$flag};
					foreach my $map(@{$info{$j}})
					{
						if ($wgs eq $map->[4])
						{
							print RD (join "\t",@$map);
							if (exists $RC{$j}{$flag})
							{
								print RD "\tRC\n";
							}
							else
							{
								print RD "\n";
							}
						}
					}
				}
			}
			else
			{
				if (!exists$roadmap{$j}{'F'} || !exists $roadmap{$j}{'R'})
				{
					print Dumper $roadmap{$j};
					next;
				}
				my ($wgsF,$pF,$strandF)=split /\:/,$roadmap{$j}{'F'};
				my ($wgsR,$pR,$strandR)=split /\:/,$roadmap{$j}{'R'};
				if (exists $roadmap{$j-1}{3})
				{
					my ($pre_wgs,$pre_p,$pre_strand)=split /\:/,$roadmap{$j-1}{3};
					if ($pre_wgs eq $wgsF)
					{
						my $rc=judge_RC($pre_p,$pre_strand,$pF,$strandF);
						if ($rc==2 || (!exists $RC{$j-1}{3} && $rc==1) || (exists $RC{$j-1}{3} && $rc==0) )
						{
							$RC{$j}{5}=1;
						}
						elsif ($rc==3)
						{
							my $tmp=$roadmap{$j}{'F'};
							$roadmap{$j}{'F'}=$roadmap{$j-1}{3};
							$roadmap{$j-1}{3}=$tmp;
						}
						$roadmap{$j}{5}=$roadmap{$j}{'F'};
						$roadmap{$j}{3}=$roadmap{$j}{'R'};
					}
					elsif($pre_wgs eq $wgsR)
					{
						my $rc=judge_RC($pre_p,$pre_strand,$pR,$strandR);
						if ($rc==2 || (!exists $RC{$j-1}{3} && $rc==1) || (exists $RC{$j-1}{3} && $rc==0) )
						{
							$RC{$j}{5}=1;
						}
						elsif ($rc==3)
						{
							my $tmp=$roadmap{$j}{'R'};
							$roadmap{$j}{'R'}=$roadmap{$j-1}{3};
							$roadmap{$j-1}{3}=$tmp;
						}
						$roadmap{$j}{5}=$roadmap{$j}{'R'};
						$roadmap{$j}{3}=$roadmap{$j}{'F'};
					}
					else
					{
						my $rc=judge_RC($pF,$strandF,$pR,$strandR);
						if ($rc==3)
						{
							$roadmap{$j}{5}=$roadmap{$j}{'R'};
							$roadmap{$j}{3}=$roadmap{$j}{'F'};
						}
						else
						{
							if ($rc==1)
							{
								$RC{$j}{5}=1;
							}
							elsif ($rc==2)
							{
								$RC{$j}{3}=1;
							}
							$roadmap{$j}{5}=$roadmap{$j}{'F'};
							$roadmap{$j}{3}=$roadmap{$j}{'R'};
						}
					}
				}
				elsif (exists $roadmap{$j+1}{5})
				{
					my ($aft_wgs,$aft_p,$aft_strand)=split /\:/,$roadmap{$j+1}{5};
					if ($aft_wgs eq $wgsF)
					{
						my $rc=judge_RC($pF,$strandF,$aft_p,$aft_strand);
						if ($rc==1)
						{
							$RC{$j}{3}=1;
						}
						elsif ($rc==3)
						{
							my $tmp=$roadmap{$j}{'F'};
							$roadmap{$j}{'F'}=$roadmap{$j+1}{5};
							$roadmap{$j+1}{5}=$tmp;
						}
						$roadmap{$j}{3}=$roadmap{$j}{'F'};
						$roadmap{$j}{5}=$roadmap{$j}{'R'};
					}
					elsif ($aft_wgs eq $wgsR)
					{
						my $rc=judge_RC($pR,$strandR,$aft_p,$aft_strand);
						if ($rc==1)
						{
							$RC{$j}{3}=1;
						}
						elsif ($rc==3)
						{
							my $tmp=$roadmap{$j}{'R'};
							$roadmap{$j}{'R'}=$roadmap{$j+1}{5};
							$roadmap{$j+1}{5}=$tmp;
						}
						$roadmap{$j}{3}=$roadmap{$j}{'R'};
						$roadmap{$j}{5}=$roadmap{$j}{'F'};
					}
					else
					{
						my $rc=judge_RC($pF,$strandF,$pR,$strandR);
						if ($rc==3)
						{
							$roadmap{$j}{5}=$roadmap{$j}{'R'};
							$roadmap{$j}{3}=$roadmap{$j}{'F'};
						}
						else
						{
							if ($rc==1)
							{
								$RC{$j}{5}=1;
							}
							elsif ($rc==2)
							{
								$RC{$j}{3}=1;
							}
							$roadmap{$j}{5}=$roadmap{$j}{'F'};
							$roadmap{$j}{3}=$roadmap{$j}{'R'};
						}
					}
				}
				else
				{
					my $rc=judge_RC($pF,$strandF,$pR,$strandR);
					if ($rc==3)
					{
						$roadmap{$j}{5}=$roadmap{$j}{'R'};
						$roadmap{$j}{3}=$roadmap{$j}{'F'};
					}
					else
					{
						if ($rc==1)
						{
							$RC{$j}{5}=1;
						}
						elsif ($rc==2)
						{
							$RC{$j}{3}=1;
						}
						$roadmap{$j}{5}=$roadmap{$j}{'F'};
						$roadmap{$j}{3}=$roadmap{$j}{'R'};
					}
				}
				for my $flag(5,3)
				{
					my ($wgs,$p,$strand)=split /\:/,$roadmap{$j}{$flag};
					foreach my $map(@{$info{$j}})
					{
						if ($wgs eq $map->[4])
						{
							print RD (join "\t",@$map);
							if (exists $RC{$j}{$flag})
							{
								print RD "\tRC\n";
							}
							else
							{
								print RD "\n";
							}
						}
					}
				}
			}
		}
		print RD "\n";
	}
	close RD;
}

##ctg164
#156000  287000  EM029A17        R       MW035431_length=1687    -       1308    1439
#156000  287000  EM029A17        F       MW043459_length=12348   -       1       108
#156000  287000  EM029A17        F       MW044473_length=46961   -       1729    2005
#321000  421000  EM003D23        R       MW015018_length=20108   +       652     990
#321000  421000  EM003D23        F       MW041838_length=35600   -       24075   24775

##ctg234
#0000    133000  EM034M10        R       MW014927_length=33280   +       10497   11151
#0000    133000  EM034M10        F       MW042441_length=43174   -       39146   39417
#106000  224000  EM043G10        F       MW042441_length=43174   +       19300   20000
#106000  224000  EM043G10        R       MW041260_length=83870   +       72922   73437

sub outputScaf
{
	my ($roadmap,$wgscontig,$Seqhash)=@_;
	print STDERR (&get_time,": output scaffolding sequence\n") if ($Verbose);
	if (keys %Seq==0)
	{
		open (IN,$wgscontig) || die $!;
		my $name="";
		while(<IN>)
		{
			if (/^\>(\S+)/)
			{
				$name=$1;
				$name=~s/\_length\=\S+$//;
			}
			else
			{
				s/\s+//g;
				$$Seqhash{$name}.=$_;
			}
		}
		close IN;
	}
	my %Group;
	my $ctg_num=0;
	my $scf_num=0;
	my $key="";
	open (RD,$roadmap) || die $!;
	my $clone="";
	while(<RD>)
	{
		if ($_=~/^\s+/)
		{
			next;
		}
		elsif (/^\#ctg(\d+)/)
		{
			$ctg_num=$1;
			$scf_num++;
			$key="$ctg_num.$scf_num";
		}
		else
		{
			chomp;
			my @t=split;
			if ($clone eq $t[2])
			{
				if ($t[-1] eq "RC")
				{
					push @{$Group{$key}},[$t[1],$t[4],$1,0,$t[6],"$t[2]:$t[3]\@$t[4]:$t[5]:$t[6]..$t[7](RC)"] if ($t[4]=~/length\=(\d+)/);
				}
				else
				{
					push @{$Group{$key}},[$t[1],$t[4],0,$1,$t[7],"$t[2]:$t[3]\@$t[4]:$t[5]:$t[6]..$t[7]"] if ($t[4]=~/length\=(\d+)/);
				}
			}
			else
			{
				if ($key eq "" || !exists $Group{$key} || (@{$Group{$key}}>0 && ${$Group{$key}}[-1][1] ne $t[4]) )
				{
					$scf_num++;
					$key="$ctg_num.$scf_num";
				}
				if ($t[-1] eq "RC")
				{
					push @{$Group{$key}},[$t[0],$t[4],$1,0,$t[7],"$t[2]:$t[3]\@$t[4]:$t[5]:$t[6]..$t[7](RC)"] if ($t[4]=~/length\=(\d+)/);
				}
				else
				{
					push @{$Group{$key}},[$t[0],$t[4],0,$1,$t[6],"$t[2]:$t[3]\@$t[4]:$t[5]:$t[6]..$t[7]"] if ($t[4]=~/length\=(\d+)/);
				}
			}
			$clone=$t[2];
		}
	}
	close RD;
	my %nr;
	open (SCF,">$Scaf_fa") || die $!;
	open (ST,">$Stat_report") if (defined $Stat_report);
	$scf_num=1;
	foreach my $key(sort{$a<=>$b}keys %Group)
	{
		my %nr_name=();
		@{$Group{$key}}=sort{$a->[0]<=>$b->[0]}@{$Group{$key}};
		my ($ctg,$scf)=split /\./,$key;
		$scf_num=sprintf("%06d",$scf_num);
		my $name=">scaffold$scf_num ctg$ctg ${$Group{$key}}[0][-1]";
		$nr_name{${$Group{$key}}[0][1]}=1;
		my $seq="";
		my %wgscontig=();
		if (${$Group{$key}}[0][2]<${$Group{$key}}[0][3])
		{
			my $seqname=${$Group{$key}}[0][1];
			$seqname=~s/\_length\=\S+$//;
			$seq=$$Seqhash{$seqname};
		}
		else
		{
			my $seqname=${$Group{$key}}[0][1];
			$seqname=~s/\_length\=\S+$//;
			$seq=reverse($$Seqhash{$seqname});
			$seq=~tr/ACGTacgt/TGCAtgca/;
		}
		$wgscontig{${$Group{$key}}[0][1]}=1;
		for (my $i=1;$i<@{$Group{$key}};$i++)
		{
			if (${$Group{$key}}[$i][1] ne ${$Group{$key}}[$i-1][1])
			{
				$name.=" -- ${$Group{$key}}[$i][-1]";
				$nr_name{${$Group{$key}}[$i][1]}=1;
				if (!exists $wgscontig{${$Group{$key}}[$i][1]})
				{
					my $gap=${$Group{$key}}[$i][0]-${$Group{$key}}[$i-1][0]-(${$Group{$key}}[$i-1][2]-${$Group{$key}}[$i-1][4])-${$Group{$key}}[$i][4];
					$seq.="N"x$gap;
					if (${$Group{$key}}[$i][2]<${$Group{$key}}[$i][3])
					{
						my $seqname=${$Group{$key}}[$i][1];
						$seqname=~s/\_length\=\S+$//;
						$seq.=$$Seqhash{$seqname};
					}
					else
					{
						my $seqname=${$Group{$key}}[$i][1];
						$seqname=~s/\_length\=\S+$//;
						my $temp=reverse($$Seqhash{$seqname});
						$temp=~tr/ACGTacgt/TGCAtgca/;
						$seq.=$temp;
					}
				}
			}
			else
			{
				$name.=" -- ${$Group{$key}}[$i][-1]";
				$nr_name{${$Group{$key}}[$i][1]}=1;
			}
			$wgscontig{${$Group{$key}}[$i][1]}=1;
		}
		if (!exists $nr{(join " ",keys %nr_name)})
		{
			print ST ("scaffold$scf_num\tlength=ctg$ctg\tscaffolding_ctg_num: ",scalar(keys %wgscontig),"\tLength: ",length($seq),"\n") if (defined $Stat_report);
			print SCF "$name\n";
			Display(\$seq);
			print SCF $seq;
			$scf_num++;
		}
		$nr{(join " ",keys %nr_name)}=1;
	}
	close SCF;
	close ST;
}

sub judge_RC
{
	my ($p1,$s1,$p2,$s2)=@_;
	my $rc=0;
	if ($p1 eq "F" && $s1 eq "+" && $p2 eq "R" && $s2 eq "+")
	{
		$rc=2;
	}
	elsif ($p1 eq "F" && $s1 eq "-" && $p2 eq "R" && $s2 eq "-")
	{
		$rc=1;
	}
	elsif ($p1 eq "R" && $s1 eq "+" && $p2 eq "F" && $s2 eq "+")
	{
		$rc=2;
	}
	elsif ($p1 eq "R" && $s1 eq "-" && $p2 eq "F" && $s2 eq "-")
	{
		$rc=1;
	}
	elsif ($p1 eq "F" && $s1 eq "-" && $p2 eq "F" && $s2 eq "+")
	{
		$rc=3;
	}
	elsif ($p1 eq "R" && $s1 eq "-" && $p2 eq "F" && $s2 eq "+")
	{
		$rc=3;
	}
}


sub check_name
{
	my @file=@_;
	my $i=0;
	for (my $j=0;$j<@file;$j++)
	{
		my $name=`grep ">" $file[$j] |head -1`;
		print STDERR (&get_time,": checking name of $file[$j]: $name") if ($Verbose);
		if ($name=~/^\S+length\=\d+/i)
		{
			$i+=2**$j;
		}
	}
	print STDERR (&get_time,": checked name status: $i\n") if ($Verbose);
	return $i;
}

sub rename_fa
{
	my ($fa,$out)=@_;
	my ($name,$seq)=("","");
	open (IN,$fa) || die "failed to read $fa\n";
	print STDERR (&get_time,": reading $fa\n") if ($Verbose);
	open (OUT,">$out") || die "failed to write to $out\n";
	print STDERR (&get_time,": outputing $out\n") if ($Verbose);
	my $SeqName="";
	while(<IN>)
	{
		s/\s*$//g;
		next if (($_ eq "")||($_=~/^\s/));
		if (/^\>(.*)/)
		{
			$SeqName=$1;
			if (($seq ne "")&&(length($seq)>0))
			{
				my $tmp=$seq;
				$tmp=~s/N//gi;
				my $len=length($tmp);
				$name="$SeqName\_length=".$len unless ($name =~ m/\_length\=\d+/);
				#$$contig_seq{$SeqName}=$seq if (defined $contig_seq);
				print OUT ">$name\n";
				Display(\$seq);
				print OUT $seq;
				$seq="";
			}
			if ($SeqName=~/(\w+)\-M([FR])/)
			{
				my ($clone,$S)=($1,$2);
				$name=(length($clone)>8)?substr($clone,0,8)."_".$S."_".substr($clone,8,length($clone)-8):"$clone\_$S";
			}
			elsif ($SeqName=~/(gi)\|(\d+)/)
			{
				$name=$1."_".$2;
			}
			elsif ($SeqName=~/(\S+)/)
			{
				$name=$1;
			}
		}
		else
		{
			$seq.=$_;
		}
	}
	if ($seq ne "")
	{
		my $tmp=$seq;
		$tmp=~s/N//gi;
		my $len=length($tmp);
		#$$contig_seq{$SeqName}=$seq  if (defined $contig_seq);
		$name.="_length=".$len unless ($name =~ m/\_length\=\d+/);
		print OUT ">$name\n";
		Display(\$seq);
		print OUT $seq;
	}
	close IN;
	close OUT;
}

sub fpc_name
{
	my ($fpc,$bacend)=@_;
	my %hash;
	open (IN1,$fpc) || die $!;
	print STDERR (&get_time,": reading $fpc\n") if ($Verbose);
	my %CG;
	my ($pre,$bac,$ctg,$left,$right)=("","","","","");
	while(<IN1>)
	{
		if (/^BAC[^\"]*\"([^\"\s]+)\"/)
		{
			if ($1 ne $pre)
			{
				($pre,$bac,$ctg,$left,$right)=("","","","","");
			}
			$bac=$1;
			$_=<IN1>;
			if ($_=~/Map[^\"]*"([^\"\s]+)\"/)
			{
				$ctg=$1;
				$hash{$bac}=$ctg;
				if($_=~/Left\s+(\d+\.\d+)$/)
				{
					$left=$1;$left=~s/\.//;
					$hash{$bac}.="_$left";
					$_=<IN1>;
					if ($_=~/Map\s+\"([^\"\s]+)\"/)
					{
						$hash{$bac}.="-$1";
					}
					if($_=~/Right\s+(\d+\.\d+)/)
					{
						my $right=$1;$right=~s/\.//;
						if (defined $bac && defined $ctg && defined $left && $left ne "")
						{
							push @{$CG{$ctg}},[$bac,$left,$right];
						}
						$hash{$bac}.="_$right";
						($left,$right)=("","");
					}
				}
				elsif($_=~/Right\s+(\d+\.\d+)/)
				{
					$right=$1;$right=~s/\.//;
					if (defined $bac && defined $ctg && defined $left && $left ne "")
					{
						push @{$CG{$ctg}},[$bac,$left,$right];
					}
					$hash{$bac}.="_$right";
					($left,$right)=("","");
				}
			}
			$pre=$bac;
		}
	}
	close IN1;
	open (CG,">$CTG_group") if (defined $CTG_group);
	print CG "#Contig Name\t# of Clones\tClone Name\tLeft Position\tRight Position\n";
	foreach my $ctg(sort keys %CG)
	{
		my $num=scalar(@{$CG{$ctg}});
		my @CGgroup=sort{$a->[1]<=>$b->[1]}@{$CG{$ctg}};
		for (my $i=0;$i<@CGgroup;$i++)
		{
			print CG ("$ctg\t$num\t",(join "\t",@{$CGgroup[$i]}),"\n");
		}
	}
	close CG;
	open (IN2,$bacend) || die $!;
	my $out=(split /\//,$bacend)[-1];
	open (OUT,">./$out.addfpcinfo") || die $!;
	print STDERR (&get_time,": changing name of $bacend according to $fpc\n") if ($Verbose);
	while(<IN2>)
	{
		if (/>([^\_]+)\_(.*)\_(length\=\d+)/)
		{
			my ($id,$other,$len)=($1,$2,$3);
			my $clone=(length($id)>8)?substr($id,0,8):$id;
			if (exists $hash{$clone})
			{
				print OUT ">$clone\_$other\_$hash{$clone}\_$len\n";
			}
			#else
			#{
			#	print OUT $_;
			#}
		}
		else
		{
			print OUT $_;
		}
	}
	close IN2;
}

sub mm
{
	my ($min,$max)=(sort{$a<=>$b}@_)[0,-1];
	return ($min,$max);
}

sub Display
{
	my $seq_p=shift;
	my $num ||= (@_) ? shift : 100;
	my $disp;
	for (my $i=0;$i<length($$seq_p);$i+=$num) {
		$disp .= substr($$seq_p,$i,$num)."\n";
	}
	$$seq_p = ($disp) ? $disp : "\n";
}

sub get_time
{
	my ($sec,$min,$hour,$mday,$mon,$year) = (localtime)[0..5];
	($sec,$min,$hour,$mday,$mon,$year) = (
		sprintf("%02d", $sec),
		sprintf("%02d", $min),
		sprintf("%02d", $hour),
		sprintf("%02d", $mday),
		sprintf("%02d", $mon + 1),
		$year + 1900
	);
	return "[$year-$mon-$mday $hour:$min:$sec]";
}

##*.fpc
#BAC : "EM001A01"
#Map "ctg84" Ends Left 15.000
#Map "ctg84" Ends Right 123.000 Oldctg 84
#Gel_number    EM001
#Bands  2337217 109
#Approximate_match_to_cosmid "EM010B24"
#Creation_date 112 1 10 10 2 
#Modified_date 112 1 10 10 15 
#  
#BAC : "EM001A02"
#Map "ctg1012" Ends Left 84.000
#Map "ctg1012" Ends Right 181.000 Oldctg 1012
#Gel_number    EM001
#Bands  1434193 98
#Creation_date 112 1 10 10 2 
#Modified_date 112 1 10 10 16 
#  
#BAC : "EM001A03"
#Map "ctg75" Ends Left 377.000
#Map "ctg75" Ends Right 478.000 Oldctg 75
#Gel_number    EM001
#Bands  2331515 102
#Creation_date 112 1 10 10 2 
#Modified_date 112 1 10 10 15

##*.pair
#EM001H09_R_10114699_Bac-RR_H06_ctg116_14000-ctg116_82000_length=465	+	465	51	464	0.89	MW003006_length=15876	15876	6664	7077	0.89	EM001H09_F_10114698_Bac-FF_G06_ctg116_14000-ctg116_82000_length=613	-	613	1	445	0.72	MW043037_length=18130	18130	4706	5150	0.72
#EM002I02_R_10114731_Bac-RR_H10_ctg1259_58000-ctg1259_152000_length=701	+	701	51	698	0.92	MW043699_length=28160	28160	23214	23861	0.92	EM002I02_F_10114730_Bac-FF_G10_ctg1259_58000-ctg1259_152000_length=701	+	701	63	701	0.91	MW041257_length=51738	51738	30909	31547	0.91
#EM003D23_R_10114739_Bac-RR_H11_ctg164_321000-ctg164_421000_length=701	+	701	40	378	0.48	MW015018_length=20108	20108	652	990	0.48	EM003D23_F_10114738_Bac-FF_G11_ctg164_321000-ctg164_421000_length=701	-	701	1	701	1	MW041838_length=35600	35600	24075	24775	1
#EM003G09_R_10114769_Bac-RR_B03_ctg2157_219000-ctg2157_349000_length=700	-	700	48	700	0.93	MW043716_length=18353	18353	10500	11152	0.93	EM003G09_F_10114768_Bac-FF_E03_ctg2157_219000-ctg2157_349000_length=384	+	384	1	382	0.99	MW042602_length=68464	68464	54226	54607	0.99
#EM003E02_R_10114771_Bac-RR_H03_ctg223_184000-ctg223_345000_length=465	-	465	28	191	0.35	MW008738_length=15645	15645	4664	4827	0.35	EM003E02_F_10114770_Bac-FF_C03_ctg223_184000-ctg223_345000_length=701	+	701	1	701	1	MW043327_length=74674	74674	53608	54304	0.99
#EM002M14_R_10114809_Bac-RR_F08_ctg4_738000-ctg4_803000_length=701	-	701	53	701	0.92	MW043894_length=86832	86832	66929	67577	0.92	EM002M14_F_10114808_Bac-FF_E08_ctg4_738000-ctg4_803000_length=701	+	701	1	690	0.98	MW041689_length=75759	75759	59777	60466	0.98
#EM002P10_R_10114855_Bac-RR_C04_ctg363_37000-ctg363_148000_length=701	+	701	55	701	0.92	MW041363_length=93301	93301	54481	55125	0.92	EM002P10_F_10114854_Bac-FF_C02_ctg363_37000-ctg363_148000_length=701	+	701	124	701	0.82	MW041789_length=79394	79394	18261	18838	0.82
#EM003F24_R_10114883_Bac-RR_H05_ctg900_85000-ctg900_195000_length=701	-	701	236	701	0.66	MW042092_length=43544	43544	8767	9232	0.66	EM003F24_F_10114882_Bac-FF_G05_ctg900_85000-ctg900_195000_length=701	-	701	1	701	1	MW043200_length=33268	33268	14697	15397	1
#EM004A04_R_10114937_Bac-RR_B04_ctg593_468000-ctg593_570000_length=701	+	701	28	701	0.96	MW041653_length=63502	63502	30437	31111	0.96	EM004A04_F_10114936_Bac-FF_A04_ctg593_468000-ctg593_570000_length=701	-	701	1	701	1	MW042211_length=131195	131195	39892	40592	1
#EM004C14_R_10114951_Bac-RR_H05_ctg12_173000-ctg12_246000_length=701	+	701	49	701	0.93	MW042975_length=65142	65142	57735	58387	0.93	EM004C14_F_10114950_Bac-FF_G05_ctg12_173000-ctg12_246000_length=700	+	700	378	700	0.46	MW043549_length=58681	58681	4506	4828	0.46

##*.single
#EM001D16_F_10114652_Bac-FF_A01_ctg797_122000-ctg797_249000_length=700	-	700	1	478	0.68	MW044095_length=146554	146554	24327	24804	0.68
#EM001C04_R_10114655_Bac-RR_D01_ctg262_0000-ctg262_96000_length=700	+	700	151	611	0.66	MW001586_length=63708	63708	12807	13266	0.66
#EM001L15_F_10114662_Bac-FF_C02_ctg636_0000-ctg636_107000_length=700	-	700	1	141	0.2	MW044884_length=12780	12780	12944	13084	0.2
#EM001G22_R_10114667_Bac-RR_H02_ctg445_123000-ctg445_247000_length=700	+	700	58	594	0.77	MW000019_length=12412	12412	1092	1629	0.77
#EM001E05_R_10114669_Bac-RR_B03_ctg1883_0000-ctg1883_24000_length=701	+	701	501	696	0.28	MW021834_length=5145	5145	5084	5280	0.28
#EM001L19_F_10114678_Bac-FF_C04_ctg1804_28000-ctg1804_123000_length=700	+	700	1	412	0.59	MW017549_length=6239	6239	102	514	0.59
#EM002C17_R_10114689_Bac-RR_F05_ctg662_405000-ctg662_526000_length=129	-	129	9	128	0.92	MW011319_length=7672	7672	7466	7585	0.92
#EM001N09_F_10114692_Bac-FF_A06_ctg53_1066000-ctg53_1207000_length=543	+	543	1	353	0.65	MW044006_length=77617	77617	78296	78648	0.65
#EM001K12_F_10114696_Bac-FF_E06_ctg937_0000-ctg937_99000_length=176	+	176	1	175	0.99	MW000847_length=43955	43955	8466	8640	0.99
#EM001H09_F_10114698_Bac-FF_G06_ctg116_14000-ctg116_82000_length=613	-	613	1	445	0.72	MW043037_length=18130	18130	4706	5150	0.72