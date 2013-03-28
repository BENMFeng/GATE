#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use lib "$FindBin::Bin/../lib";
use File::Basename qw(basename dirname);
use SVG;
use Font;

=head1 name

Copyright (c) 2012 BENM(Binxiao) Feng                            
All Rights Reserved.                                                   
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

 GenePlot.pl -- PERL script with SVG package to draw alternative
                splicing events in order to compare to refGene
                and transcripts
 Prviouse Version PlotGeneSplicing.pl

Note: Require PERL v5.10.1

=head1 Description

Alternative Splicing Events Type:

SE   -- skipped exons

MXE  -- mutually exclusive exons

A5SS -- alternative 5' splicing sites

A3SS -- alternative 3' splicing sites

AFE  -- alternative first exons

ALE  -- alternative last exons

RI   -- retained intron

=head1 Usage

 % perl GenePlot.pl				[options]

   <I/O>
    --list <infile>				input gene list file
    --color_scheme				input color scheme configure file
    --ref <infile>				input refGene file in GFF format
    --as <label1:infile1,label2:infile2,...>	input alternative splicing events diagnose result files with their label name
    --tr <label1:infile1,label2:infile2,...>	input transcripts GTF/GFF format files with their label name
    --bam <label1:bam1,label2:bam2,...>		input alignment file in BAM format with their label name
    --jun <label1:bed1,label2:bed2,...>		input AS junction file in BED format with their label name, need disposeOverlap.pl
    --outdir <ouput_folder>			ouput folder

   [Parameters]
    --fpkm <float>				FPKM lowest limitet of isform, default: 0.1
    --sp <int>					set splicing reads limited for ploting: 2
    --only <str>				only plot for one kind of AS, default: null
    --uniq					jus stat for uniq mapping reads, default: null
    --mapq					filter mapping qualtiy, default: 10
    --depthunit					set depth show scale axis unit: 'int' or 'logint', 'log', 'ln', default: 1

    --instruction				show file instruction
    --verbose					print running status
    --help					show this usage

=head1 author

	BENM <BinxiaoFeng@gmail.com>

=head1 version

Prvious Version

	PlotGeneSplicing.pl v1.0a, 2012-03-14
	PlotGeneSplicing.pl v1.1a, 2012-03-15 add coverage depth
	PlotGeneSplicing.pl v1.2a, 2012-03-18 revise the bug of coverage depth drawing
	PlotGeneSplicing.pl v1.3a, 2012-03-19 revise the bug of parse_coverage
	PlotGeneSplicing.pl v1.4a, 2012-03-22 add junction reads info for proving skipped exon
	PlotGeneSplicing.pl v1.5a, 2012-03-29 update
	PlotGeneSplicing.pl v1.6a, 2012-10-11 update

Current Version

	GenePlot.pl v0.1.0, 2012-11-14 revise the most of functions

=head1 example

	perl GenePlot.pl --list gene_list.xls --color_scheme custom_color.conf \
	--fpkm 0.1 --only SE:MXE --uniq --mapq 20 --ref refGene.gff \
	--tr TS1:S1_transctripts.gtf,TS2:S2_transcripts.gtf,TS3:S3_transcripts.gtf\
	--as AS1:AS1_result.gff,AS2:AS2_result.gff,AS3:AS3_result.gff \
	--bam BS1:BS1.bam,BS2:BS2.bam,BS3:BS3.bam \
	--jun JS1:JS1.bed,JS2:JS2.bed,JS3:JS3.bed

=cut

my ($GeneList,$ColorScheme,$RefGene,$Transcripts,$ASresult,$BAMfile,$ASjun,$Outdir,
$FPKM,$SPreads,$Only,$Uniq,$MapQ,$DepthUnit,$Instruction,$Verbose,$Help);
GetOptions(
	'list:s'=>\$GeneList,
	"color_scheme:s"=>\$ColorScheme,
	'ref:s'=>\$RefGene,
	'tr:s'=>\$Transcripts,
	'as:s'=>\$ASresult,
	"bam:s"=>\$BAMfile,
	"jun:s"=>\$ASjun,
	"outdir:s"=>\$Outdir,
	"fpkm:s"=>\$FPKM,
	"sp:i"=>\$SPreads,
	"only:s"=>\$Only,
	"mapq:i"=>\$MapQ,
	"uniq"=>\$Uniq,
	"depthunit"=>\$DepthUnit,
	"verbose"=>\$Verbose,
	'help'=>\$Help
);
die `pod2text $0` if (($Help)||(!defined $RefGene)||(!defined $GeneList));

$Outdir ||= "svg_out";
`mkdir -p $Outdir` if (!-d $Outdir);

# default parameters
$FPKM ||= 0.1;
$SPreads ||= 2;
$MapQ ||= 10;
$DepthUnit ||= 1;

# global parameters
my $Figure_resolution = 1/10;		## 1 pixel stands for 10 bp in displaying
my $Vertical_skip = 40;			## set vertical distance between neighbouring figure units
my $Left_margin = 150;
my $Right_margin = 100;
my $Top_margin = 50;
my $Bottom_margin = 50;
my $Font_family = "ArialNarrow";	##font family for all the texts
my $Font_size = 8;			##font size for gene name 
my $Font_size_title = 30;		##font size for chromosome title
my $Font_size_mark = 8;		##font size for the marks in the left border
my $GeneBlock_Height = 16;
my $Coverage_figure_height = 0;

#default color scheme
my %Color = (
	'gene' => {
		'exon' => 'blue',
		'CDS' => 'blue',
		'start_codon' => 'blue',
		'end_codon' => 'blue',
		'transcript_start_site' => '',
		'transcript_end_sit' => '',
		'gene' => 'blue',
		'mRNA' => 'blue',
	},
	'refGene' => {
		'color' => 'blue',
	},
	'0' => 'black',
	'1' => 'red',
	'2' => 'orange',
	'3' => 'green',
	'4' => 'yellow',
	'5' => 'cyan',
	'6' => 'grey',
	'7' => 'pink',
	'8' => 'chocolate',
	'9' => 'wheat',
	'10' => 'black',
	'transcript' => 'violet'
);

my %GeneLH;
my %GeneRH;
my %GeneBH;
my @Label;
my %lh;
my %Data;
my %Junction;
my %jun_gene;
my %BAM_hash;
my %OnlyDetect;

# main

setupcustomizeColor($ColorScheme,\%Color) if (defined $ColorScheme && -s $ColorScheme);

parseOnlyDetectAStypes($Only,\%OnlyDetect) if (defined $Only && $Only ne "");

parseList($GeneList,\%GeneLH);

parseGff($RefGene,"refGene",\%Data);
push @Label,"refGene";

my @AS=(split /\,/,$ASresult) if (defined $ASresult);
if (@AS>0)
{
	foreach my $as(@AS)
	{
		my ($label,$file)=split /\:/,$as;
		if (-f $file)
		{
			parseGff($file,$label,\%Data);
			push @Label,$label unless (exists $lh{$label});
			$lh{$label}=1;
		}
	}
}

my @Tr=(split /\,/,$Transcripts) if (defined $Transcripts);
if (@Tr>0)
{
	foreach my $tr(@Tr)
	{
		my ($label,$file)=split /\:/,$tr;
		if (-f $file)
		{
			parseGff($file,$label,\%Data);
			push @Label,$label unless (exists $lh{$label});
			$lh{$label}=1;
		}
	}
}

my @Junc=(split /\,/,$ASjun)  if (defined $ASjun);
if (@Junc>0)
{
	foreach my $aj(@Junc)
	{
		my ($label,$file)=split /\:/,$aj;
		if (-f $file)
		{
			parseBed($file,$RefGene,$label,\%Junction);
		}
	}
	$Left_margin = 300;
	$Top_margin = 300;
	$Bottom_margin = 100;
}

if ((defined $BAMfile)&&($BAMfile=~/\S+\:\S+/))
{
	while($BAMfile=~/([^\:\,]+)\:([^\:\,]+)/g)
	{
		$BAM_hash{$1}=$2;
		push @Label,$1 unless (exists $lh{$1} && @Label>0);
		$lh{$1}=1;
	}
	$Coverage_figure_height=50;
}
plotSVGfigure(\%Data);

sub plotSVGfigure
{
	my $data_p = shift;
	print STDERR readTime()."ploting SVG figure...\n" if ($Verbose);
	foreach my $gene(keys %$data_p)
	{
		#next if ($GeneBH{$gene}<=@{$GeneRH{$gene}});
		next if (!exists $GeneRH{$gene} || !defined $GeneRH{$gene} || @{$GeneRH{$gene}}<0);
		print STDERR readTime()."\t$gene\n" if ($Verbose);
		my $chr_start=int((sort{$a->[0]<=>$b->[0]}@{$GeneRH{$gene}})[0][0]);
		my $chr_end=int((sort{$b->[1]<=>$a->[1]}@{$GeneRH{$gene}})[0][1]);
		my $chr=${$GeneRH{$gene}}[0][3];
		my $Figure_width=0;
		my $Figure_height=$GeneBH{$gene} * 2 * ($Vertical_skip + $GeneBlock_Height + $Font_size);
		my ($ruler_mark,$scalesize);

		adjustResolution($chr_end-$chr_start,\$Figure_width,\$Figure_height,\$ruler_mark,\$scalesize);

		my $juc_figure_height=0;
		if (exists $Junction{$gene})
		{
			$juc_figure_height=scalar(keys %{$jun_gene{$gene}})*$GeneBlock_Height/scalar(@Label);
			$Figure_height+=$juc_figure_height;
		}
		my $jun_cov_width=$Left_margin*1/3;
		my $jun_cov_height=$Left_margin*1/6;
		my $junction_y_skip=0;
		if (exists $Junction{$gene})
		{
			$junction_y_skip=$juc_figure_height/(scalar(keys %{$jun_gene{$gene}}));
		}
		my $junction_reads_cov_skip=$jun_cov_height+25;

		my ($y_coor_up,$y_coor_down,$figure_y_coor,$coverage_y_coor);
		#if (${$GeneRH{$gene}}[0][2] eq "+")
		#{
		#	$figure_y_coor=($Figure_height-scalar(@Label) * $Coverage_figure_height) * 0.618;
		#	($y_coor_up,$y_coor_down)=( $figure_y_coor - scalar(@Label)*$Coverage_figure_height - $juc_figure_height,($Figure_height-scalar(@Label) * $Coverage_figure_height) * 0.618);
		#}
		#else
		#{
		#	$figure_y_coor=($Figure_height-scalar(@Label) * $Coverage_figure_height) * (1-0.618);
		#	($y_coor_up,$y_coor_down)=( $figure_y_coor - scalar(@Label)*$Coverage_figure_height - $juc_figure_height,($Figure_height-scalar(@Label) * $Coverage_figure_height) * (1-0.618));
		#}
		$figure_y_coor=($Figure_height-scalar(@Label) * $Coverage_figure_height) * (1-0.618);
		($y_coor_up,$y_coor_down)=( $figure_y_coor - scalar(@Label)*$Coverage_figure_height - $juc_figure_height,($Figure_height-scalar(@Label) * $Coverage_figure_height) * (1-0.618));

		$coverage_y_coor = $figure_y_coor - 2*$Vertical_skip;
		my $junction_y_coor=0;
		if (exists $Junction{$gene})
		{
			$junction_y_coor=$coverage_y_coor - (keys %BAM_hash) * $Coverage_figure_height - $Vertical_skip;
		}
		my $junction_reads_cov_y=$junction_y_coor;

		my $svg = SVG->new('width',$Figure_width,'height',$Figure_height);
		my $bg_color = "white";
		$svg->rect('x',0, 'y',0,'width',$Figure_width,'height',$Figure_height,'stroke',$bg_color,'fill',$bg_color);

		$svg->text('x',$Figure_width - $Right_margin + textWidth($Font_family,$Font_size_mark,$ruler_mark)*0.5,'y',$figure_y_coor,'fill','black','-cdata',$ruler_mark,"font-family",$Font_family,"font-size",$Font_size_mark);
		plotRuler("svg",$svg,"Y",$figure_y_coor, "X_start",$Left_margin,"X_end",$Figure_width-$Right_margin,"bp_start",$chr_start,"bp_end",$chr_end,"scaletype","Kb","scaletypepos","right","scalestart","force","rulerstyle",2,"bigscalesize",$scalesize);
		my $Coverage_hi=0;
		my %CoverageDepth=();
		my($coverage_chr_start,$coverage_chr_end,$coverage_x_coor);

		for (my $i=0;$i<@Label;$i++)
		{
			my $label=$Label[$i];
			my $color=(exists $Color{$label}{'color'})?$Color{$label}{'color'}:$Color{($i+1)%10};
			print STDERR readTime()."\t$label : $color\n" if ($Verbose);
			if (exists $$data_p{$gene}{$label})
			{
				my @assp=();
				for my $info_p(@{$$data_p{$gene}{$label}})
				{
					my ($id,$pos1,$pos2,$strand,$attributes)=@$info_p;
					$color = $Color{$label}{$id} if (exists $Color{$label}{$id});
					#my $y_coor = ($strand eq "+") ? $y_coor_up : $y_coor_down;
					my $y_coor = $y_coor_down;
					my $x1_coor = $Left_margin + ($pos1-$chr_start)*$Figure_resolution;
					my $x2_coor = $Left_margin + ($pos2-$chr_start)*$Figure_resolution;
					if (($id=~/gene/i)||($id=~/mRNA/i)||($id=~/transcript/i)||($attributes=~/mRNA/i)||($attributes=~/gene\=/i)||($attributes=~/name\=/i))
					{
						#if ($strand eq "+")
						#{
						#	$y_coor-=$Vertical_skip;
						#	$y_coor_up = $y_coor;
						#}
						#else
						#{
						#	$y_coor+=$Vertical_skip;
						#	$y_coor_down = $y_coor;
						#}
						$y_coor+=$Vertical_skip;
						$y_coor_down = $y_coor;
						if ($id=~/transcript/)
						{
							$coverage_x_coor=((!defined $coverage_x_coor)||($x1_coor<$coverage_x_coor))?$x1_coor:$coverage_x_coor;
							$coverage_x_coor=$Left_margin if ($Left_margin<$coverage_x_coor);
							$coverage_chr_start=((!defined $coverage_chr_start)||($pos1<$coverage_chr_start))?$pos1:$coverage_chr_start;
							$coverage_chr_start=$chr_start if ($chr_start<$coverage_chr_start);
							$coverage_chr_end=((!defined $coverage_chr_start)||($pos2>$coverage_chr_start))?$pos2:$coverage_chr_end;
							$coverage_chr_end=$chr_end if ($chr_end>$coverage_chr_start);
						}
						$svg->text('x',$x1_coor, 'y',$y_coor-4,'fill','#000000','-cdata',"$label $attributes","font-family",$Font_family,"font-size",$Font_size);
						$svg->rect('x',$x1_coor, 'y',$y_coor + $GeneBlock_Height*1/2,'width',($x2_coor-$x1_coor+1),'height',$GeneBlock_Height*1/4,'stroke',$color,'fill',$color);
						if ( ($id=~/gene/i)||($id=~/mRNA/i)||($id=~/transcript/i) )
						{
							for (my $x3_coor=$x1_coor+($chr_end-$chr_start+1)*$Figure_resolution/50;$x3_coor<$x2_coor-($chr_end-$chr_start+1)*$Figure_resolution/50;$x3_coor+=($chr_end-$chr_start+1)*$Figure_resolution/50)
							{
								if ($strand eq "+")
								{
									$svg->line('x1',$x3_coor-$GeneBlock_Height*1/4,'y1',$y_coor+$GeneBlock_Height*1/4, 'x2',$x3_coor+$GeneBlock_Height*1/4,'y2',$y_coor+$GeneBlock_Height*5/8,'stroke',$color,'stroke-width',1);
									$svg->line('x1',$x3_coor-$GeneBlock_Height*1/4,'y1',$y_coor+$GeneBlock_Height, 'x2',$x3_coor+$GeneBlock_Height*1/4,'y2',$y_coor+$GeneBlock_Height*5/8,'stroke',$color,'stroke-width',1);
								}
								else
								{
									$svg->line('x1',$x3_coor-$GeneBlock_Height*1/4,'y1',$y_coor+$GeneBlock_Height*5/8, 'x2',$x3_coor+$GeneBlock_Height*1/4,'y2',$y_coor+$GeneBlock_Height*1/4,'stroke',$color,'stroke-width',1);
									$svg->line('x1',$x3_coor-$GeneBlock_Height*1/4,'y1',$y_coor+$GeneBlock_Height*5/8, 'x2',$x3_coor+$GeneBlock_Height*1/4,'y2',$y_coor+$GeneBlock_Height,'stroke',$color,'stroke-width',1);
								}
							}
						}
					}
					else
					{
						if ($id eq "alternative")
						{
							if ($attributes=~/\=RI/)
							{
								$svg->rect('x',$x1_coor, 'y',$y_coor+$GeneBlock_Height/8,'width',($x2_coor-$x1_coor+1),'height',$GeneBlock_Height,'stroke','purple','fill','purple');
							}
							else
							{
								$svg->rect('x',$x1_coor, 'y',$y_coor+$GeneBlock_Height/8,'width',($x2_coor-$x1_coor+1),'height',$GeneBlock_Height,'stroke',$color,'fill',$color,"opacity",0.2);
							}
							push @assp,[$x1_coor,$x2_coor];
						}
						elsif ($id =~ /UTR/)
						{
							$svg->rect('x',$x1_coor, 'y',$y_coor+$GeneBlock_Height*5/24,'width',($x2_coor-$x1_coor+1),'height',$GeneBlock_Height*2/3,'stroke',$color,'fill',$color);
						}
						else
						{
							$svg->rect('x',$x1_coor, 'y',$y_coor+$GeneBlock_Height/8,'width',($x2_coor-$x1_coor+1),'height',$GeneBlock_Height,'stroke',$color,'fill',$color);
						}
					}
				}
				if (exists $Junction{$gene}{$label})
				{
					$junction_y_skip=$juc_figure_height/@{$Junction{$gene}{$label}} if (@{$Junction{$gene}{$label}}>0);
					foreach my $junc(@{$Junction{$gene}{$label}})
					{
						next if ($junc->[0]<$chr_start||$junc->[1]>$chr_end||@assp<=0);
						my $jun_left1_coor=$Left_margin + ($junc->[0]-$chr_start)*$Figure_resolution;
						my $jun_left2_coor=$Left_margin + ($junc->[0]+$junc->[2]-$chr_start)*$Figure_resolution;
						my $jun_right1_coor=$Left_margin + ($junc->[1]-$junc->[3]-$chr_start)*$Figure_resolution;
						my $jun_right2_coor=$Left_margin + ($junc->[1]-$chr_start)*$Figure_resolution;
						
						my $cross_overlap=0;
						foreach my $sp(@assp)
						{
							if ($jun_left1_coor<=$sp->[1] && $jun_right2_coor>=$sp->[0])
							{
								$cross_overlap=1;
								last;
							}
						}
						next if ($cross_overlap==0);
						print STDERR readTime()."ploting junctions reads of $label $junc->[0] - $junc->[1]\n" if ($Verbose);
						$svg->rect('x',$jun_left1_coor, 'y',$junction_y_coor,'width',($jun_left2_coor-$jun_left1_coor+1),'height',$GeneBlock_Height/6,'stroke',$color,'fill',$color);
						$svg->line('x1',$jun_left2_coor,'y1',$junction_y_coor+$GeneBlock_Height/10,'x2',$jun_right1_coor,'y2',$junction_y_coor+$GeneBlock_Height/10,'stroke',$color,'stroke-width',1);
						#$svg->rect('x',$jun_left2_coor, 'y',$junction_y_coor+$GeneBlock_Height/8,'width',($jun_right1_coor-$jun_left2_coor+1),'height',$GeneBlock_Height/8,'stroke',$color,'fill',$color);
						$svg->rect('x',$jun_right1_coor, 'y',$junction_y_coor,'width',($jun_right2_coor-$jun_right1_coor+1),'height',$GeneBlock_Height/6,'stroke',$color,'fill',$color);
						$svg->text('x',$jun_right2_coor+2, 'y',$junction_y_coor+$Font_size*2/5,'fill','#000000','-cdata',$junc->[4],"font-family",$Font_family,"font-size",$Font_size*2/5);
						my $get_junc_reads=0;
						if ((defined $BAMfile)&&(exists $BAM_hash{$label}))
						{
							my $mark=$junc->[4];
							$mark.=" reads";
							my $sp_color=$Color{$label}{'splicingreads'}?($Color{$label}{'splicingreads'}):$color;
							$get_junc_reads=plot_junction_reads_coverage($svg,$BAM_hash{$label},$chr,$junc->[0],$junc->[1],($junc->[1]-$junc->[0]-$junc->[2]-$junc->[3]),$jun_left1_coor,$jun_right2_coor,$junction_y_coor,$jun_cov_width,$junction_reads_cov_y,$jun_cov_width,$jun_cov_height,$sp_color,$mark);
							if ($get_junc_reads==1)
							{
								#$junction_reads_cov_y += ($Left_margin*1/4>$junction_reads_cov_skip) ? $Left_margin*1/4 : $junction_reads_cov_skip;
								$junction_reads_cov_y -= $junction_reads_cov_skip;
							}
						}
						$junction_y_coor-=$junction_y_skip;
					}
				}
			}
			if ((defined $BAMfile)&&(exists $BAM_hash{$label}))
			{
				$coverage_chr_start = $chr_start if (!defined $coverage_chr_start);
				$coverage_chr_end   = $chr_end if (!defined $coverage_chr_end);
				my $bin=1/$Figure_resolution;
				@{$CoverageDepth{$label}}=parseCoverage($BAM_hash{$label},$chr,$coverage_chr_start,$coverage_chr_end,$bin,\$Coverage_hi);
			}
		}
		for (my $i=0;$i<@Label;$i++)
		{
			my $label=$Label[$i];
			my $windows=1;
			my $color=(exists $Color{$label}{'depth'})?$Color{$label}{'depth'}:$Color{($i+1) % 10};
			$coverage_chr_start = $chr_start if (!defined $coverage_chr_start);
			$coverage_chr_end   = $chr_end if (!defined $coverage_chr_end);
			$coverage_x_coor=$Left_margin if (!defined $coverage_x_coor);
			if ((defined $BAMfile)&&(exists $BAM_hash{$label})&&(defined $coverage_chr_start)&&(defined $coverage_chr_end))
			{
				print STDERR readTime()."ploting $label coverage...\n" if ($Verbose);
				plotCoverage($svg,\@{$CoverageDepth{$label}},$coverage_x_coor,$coverage_y_coor,$Coverage_figure_height,$Coverage_hi,$windows,$color) if ($Coverage_hi>0);
				$coverage_y_coor-=$Coverage_figure_height;
			}
		}
		print STDERR readTime()."Outupt figure: $Outdir/$gene\_AS.svg\n" if ($Verbose);
		my $FigureName="$Outdir/$gene\_AS.svg";
		open OUT,">$FigureName" || die "fail creat $FigureName";
		print OUT $svg->xmlify();
		close OUT;
		print STDERR readTime()."End\n" if ($Verbose);
	}
}

sub plot_junction_reads_coverage
{
	print STDERR readTime()." ploting junction reads coverage...\n" if ($Verbose);
	my $svg = shift;
	my $bam = shift;
	my $chr = shift;
	my $start = shift;
	my $end = shift;
	my $gap = shift;
	my $junc_x1_coor = shift;
	my $junc_x2_coor = shift;
	my $junc_y_coor = shift;
	my $x_coor = shift;
	my $y_coor = shift;
	my $width = shift;
	my $height = shift;
	my $label_color = shift;
	my $mark = shift;
	my %jun_cov;
	#my @jun_cov;
	my ($color1,$color2)=('red','blue');
	($color1,$color2)=split /\:/,$label_color if ($label_color=~/\:/);
	open (IN,"samtools view -q $MapQ -F 4 $bam $chr:$start-$end|") || die $!;
	while(<IN>)
	{
		my @t=split;
		if ($_=~/NH\:i\:(\d+)/) {
			next if ($1>1 && defined $Uniq);
		}
		#next if ($t[6] ne "=");
		my $len=length($t[9]);
		if ($t[5]=~/(\d+)[^N\d]$gap(N)(\d+)[^N\d]/)
		{
			my ($left,$right)=($1,$3);
			map{$jun_cov{$_}++}((101-$left)..100);
			map{$jun_cov{$_}++}(101..(100+$right));
		}
	}
	close IN;
	if (keys %jun_cov>0) {
		my $max_jun_cov=(sort{$b<=>$a}values %jun_cov)[0];
		$max_jun_cov=($max_jun_cov>100)?int($max_jun_cov/100)*100+int(($max_jun_cov/int($max_jun_cov/100+0.5))/10)*10:100;
		$svg->line('x1',$x_coor-$width/2,'y1',$y_coor,'x2',$x_coor-$width/2,'y2',$y_coor-$height,'stroke','#000000','stroke-width',1);
		my $vertical_y=$y_coor;
		my $k=0;
		$svg->text('x',$x_coor-$width/2 - textWidth($Font_family,$Font_size_mark,$k),'y',$vertical_y,'fill','black','-cdata',"0","font-family",$Font_family,"font-size",$Font_size_mark);
		$k+=int($max_jun_cov/3+0.5);
		$vertical_y-=$height/3;
		$svg->text('x',$x_coor-$width/2 - textWidth($Font_family,$Font_size_mark,$k),'y',$vertical_y,'fill','black','-cdata',"50","font-family",$Font_family,"font-size",$Font_size_mark);
		$k+=int($max_jun_cov/3+0.5);
		$vertical_y-=$height/3;
		$svg->text('x',$x_coor-$width/2 - textWidth($Font_family,$Font_size_mark,$k),'y',$vertical_y,'fill','black','-cdata',"100","font-family",$Font_family,"font-size",$Font_size_mark);
		$k+=int($max_jun_cov/3+0.5);
		$vertical_y-=$height/3;
		$svg->text('x',$x_coor-$width/2 - textWidth($Font_family,$Font_size_mark,$k),'y',$vertical_y,'fill','black','-cdata',"150","font-family",$Font_family,"font-size",$Font_size_mark);
		my $jun_x_coor=$x_coor+2;
		my $bin_width=$width/202;
		for my $i(1..202){
			if (exists $jun_cov{$i})
			{
				my $cov=$jun_cov{$i};
				if ($i<101)
				{
					$svg->rect('x',$jun_x_coor, 'y',$y_coor-$cov*$height/$max_jun_cov,'width',$bin_width,'height',$cov*$height/$max_jun_cov,'stroke',$color1,'fill',$color1);
				}
				else
				{
					$svg->rect('x',$jun_x_coor, 'y',$y_coor-$cov*$height/$max_jun_cov,'width',$bin_width,'height',$cov*$height/$max_jun_cov,'stroke',$color2,'fill',$color2);
				}
			}
			$jun_x_coor+=$bin_width;
		}
		$svg->text('x',$x_coor+$width/2 - textWidth($Font_family,$Font_size_mark,$mark)/2,'y',$vertical_y,'fill','black','-cdata',$mark,"font-family",$Font_family,"font-size",$Font_size_mark);
		return 1;
	}
	else
	{
		return 0;
	}
}

sub adjustResolution
{
	my ($scaleWidth,$figure_width,$figure_height,$ruler_mark,$scalesize)=@_;
	if ($scaleWidth<20000)
	{
		$Figure_resolution=1/10;
	}
	elsif ($scaleWidth<200000)
	{
		$Figure_resolution=1/100;
	}
	elsif ($scaleWidth<20000000)
	{
		$Figure_resolution=1/1000;
	}
	elsif ($scaleWidth<200000000)
	{
		$Figure_resolution=1/10000;
	}
	else
	{
		$Figure_resolution=1/100000;
	}
	$$figure_width += $Left_margin + $scaleWidth * $Figure_resolution + $Right_margin;
	$$figure_height += $Top_margin + scalar(@Label)*($Coverage_figure_height+$Vertical_skip) + $Bottom_margin;
	
	$$ruler_mark = "10bp" if ($Figure_resolution==1);
	$$ruler_mark = "100bp" if ($Figure_resolution==1/10);
	$$ruler_mark = "Kb" if ($Figure_resolution==1/100);
	$$ruler_mark = "10Kb" if ($Figure_resolution==1/1000);
	$$ruler_mark = "100Kb" if ($Figure_resolution==1/10000);
	$$ruler_mark = "Mb" if ($Figure_resolution==1/100000);
	
	$$scalesize = 10 if ($$ruler_mark eq "10bp");
	$$scalesize = 100 if ($$ruler_mark eq "100bp");
	$$scalesize = 1000 if ($$ruler_mark eq "Kb");
	$$scalesize = 10000 if ($$ruler_mark eq "10Kb");
	$$scalesize = 100000 if ($$ruler_mark eq "100Kb");
	$$scalesize = 1000000 if ($$ruler_mark eq "Mb");
}


sub plotCoverage
{
	my $svg = shift;
	my $ary_p = shift;
	my $X_coor = shift;
	my $Y_coor = shift;
	my $coverage_figure_heigth = shift;
	my $coverage_max = shift;
	my $win = shift;
	my $color = shift;
	my $coverage_resolution = $coverage_figure_heigth / $coverage_max;

	my $x_coor = $X_coor;
	for (my $i=0; $i<@$ary_p; $i++)
	{
		$ary_p->[$i] = $coverage_max if($ary_p->[$i] > $coverage_max); 
		if ($ary_p->[$i]==0)
		{
			$x_coor+=$win;
			next;
		}
		my $y_coor = $Y_coor + $coverage_figure_heigth - $ary_p->[$i] * $coverage_resolution;
		if ($DepthUnit eq "int")
		{
			$svg->rect('x',$x_coor, 'y',$y_coor,'width',$win,'height',$ary_p->[$i] * $coverage_resolution,'stroke',$color,'fill',$color);
		}
		elsif ($DepthUnit eq "log")
		{
			$coverage_resolution = $coverage_figure_heigth / log($coverage_max);
			$svg->rect('x',$x_coor, 'y',$y_coor,'width',$win,'height',log($ary_p->[$i]) * $coverage_resolution,'stroke',$color,'fill',$color);
		}
		elsif ($DepthUnit eq "ln")
		{
			$coverage_resolution = $coverage_figure_heigth * log(10) / log($coverage_max);
			$svg->rect('x',$x_coor, 'y',$y_coor,'width',$win,'height',log($ary_p->[$i])/log(10) * $coverage_resolution,'stroke',$color,'fill',$color);
		}
		elsif ($DepthUnit =~ /log(\d+)/)
		{
			$coverage_resolution = $coverage_figure_heigth * log($1) / log($coverage_max);
			$svg->rect('x',$x_coor, 'y',$y_coor,'width',$win,'height',log($ary_p->[$i])/log($1) * $coverage_resolution,'stroke',$color,'fill',$color);
		}
		elsif ($DepthUnit =~ /^(\d+)$/)
		{
			$coverage_resolution = $coverage_figure_heigth * $1 / $coverage_max;
			$svg->rect('x',$x_coor, 'y',$y_coor,'width',$win,'height',$ary_p->[$i]/$1 * $coverage_resolution,'stroke',$color,'fill',$color);
		}
		elsif ($DepthUnit =~ /^(\d+)max/)
		{
			$ary_p->[$i]=($ary_p->[$i]>$1)?$1:$ary_p->[$i];
			$svg->rect('x',$x_coor, 'y',$y_coor,'width',$win,'height',$ary_p->[$i] * $coverage_resolution,'stroke',$color,'fill',$color);
		}
		$x_coor+=$win;
	}
}

sub parseCoverage
{
	my ($bam_file,$chr,$chr_start,$chr_end,$box_bin,$max_value)=@_;
	print STDERR readTime()."parse coverage of region:$chr:$chr_start-$chr_end from $bam_file\n" if ($Verbose);
	if (!-e "$bam_file.bai")
	{
		print STDERR readTime()."No bam index! Try:\n\tsamtools index $bam_file\n" if ($Verbose);
		`samtools index $bam_file`;
	}
	open (IN,"samtools view -q $MapQ -F 4 $bam_file $chr:$chr_start-$chr_end|");
	my %Depth=();
	#$#Depth=$chr_end-$chr_start+1;
	my $lines=0;
	while (<IN>)
	{
		my @t=split /\s+/,$_;
		if ($_=~/NH\:i\:(\d+)/) {
			next if ($Uniq && $1>1);
		}
		#next if ($t[6] ne "=");
		if ($t[5]=~/^(\d+)M$/)
		{
			my $match=$1;
			map{$Depth{$_-$chr_start}+=1}($t[3]..($t[3]+$match));
		}
		else
		{
			my $loc=$t[3];
			my $pos=$t[3];
			if ($t[5] eq "=")
			{
				map{$Depth{$_-$chr_start}+=1}($pos..($pos+length($t[9])-1));
			}
			else
			{
				while($t[5]=~/(\d+)(\D)/g)
				{
					my ($match,$cigar)=($1,$2);
					if (($cigar eq "M")||($cigar eq "S"))
					{
						$pos+=$match;
					}
					elsif ($cigar ne "I")
					{
						if ($pos>$loc){map{$Depth{$_-$chr_start}+=1}($loc..($pos-1));}
						$pos+=$match;
						$loc=$pos+1;
					}
					else
					{
						if ($pos>$loc){map{$Depth{$_-$chr_start}+=1}($loc..($pos-1));}
						if ($pos>$loc){$loc=$pos}
					}
				}
				if ($pos>$loc){map{$Depth{$_-$chr_start}+=1}($loc..($pos-1));}
			}
		}
	}
	close IN;
	
	my @coverage_array=();
	
	#my $box_bin = int(1/$Figure_resolution);
	
	for (my $i = 0; $i < $chr_end-$chr_start+1 ; $i += $box_bin)
	{
		my $coverage_median;
		my @box=();
		if ($i>$box_bin && $i<$chr_end-$chr_start+1-$box_bin)
		{
			for (my $j = $i-$box_bin/2; $j < $i + $box_bin/2; $j++)
			{
				if (exists $Depth{$j})
				{
					push @box,$Depth{$j};
				}
				else
				{
					push @box,0;
				}
			}
		}
		else
		{
			for (my $j = $i; $j < $i + $box_bin; $j++)
			{
				if (exists $Depth{$j})
				{
					push @box,$Depth{$j};
				}
				else
				{
					push @box,0;
				}
			}
		}
		@box=sort{$a<=>$b}@box;
		my $number=scalar @box;
		$coverage_median=$box[int $number/2];
		push @coverage_array, $coverage_median;
		$$max_value = ($$max_value>$coverage_median) ? $$max_value : $coverage_median;
	}

	return @coverage_array;
}

#ID	FirstTableID	Start	End	Length	OverlapNumber	OverlapSize	OverlapRate	SecondTableID:Start,End,Length,OverlapSize...
#chr10|+	JUNC00022798|65|98,57	6148103	6150705	2603	2	5206	2	mRNA|ID=NM_032905; name=RBM17;:6130949,6159422,28474,2603	mRNA|ID=NM_001145547; name=RBM17;:6131309,6159422,28114,2603
sub parseBed
{
	my ($bed,$gff,$label,$hash)=@_;
	print STDERR readTime()."parsing BED file: $bed\n" if ($Verbose);
	my $overlap=join "_overlap_",((split /\//,$bed)[-1],(split /\//,$gff)[-1]);
	print STDERR readTime()."disposing overlap between file: $bed and file: $gff\n" if ($Verbose);
	`disposeOverlap.pl --E O --i1 $bed --f1 0,5-3,4,10-1-2  --i2 $gff --f2 0,6-2,8-3-4 > $overlap` unless (-s $overlap);
	print STDERR "disposeOverlap.pl --E O --i1 $bed --f1 0,5-3,4,10-1-2  --i2 $gff --f2 0,6-2,8-3-4 > $overlap\n" if ($Verbose);

	my @sp_jun=();
	my ($pre_chr,$pre_strand)=("","");
	my ($skip_s,$skip_e)=(0,0);
	open (IN,$overlap) || die $!;
	<IN>;
	while(<IN>)
	{
		chomp;
		my @t=split /\t+/,$_;
		my ($seqid,$strand)=($1,$2) if ($t[0]=~/(\S+)\|([\+\-])$/);
		my ($juncid,$num,$len1,$len2)=($1,$2,$3,$4) if ($t[1]=~/(\S+)\|(\d+)\|(\d+)\,(\d+)$/);
		next if ($num<$SPreads);
		my %Skip;
		for (my $i=8;$i<@t;$i++)
		{
			if ($t[$i]=~/exon.*\:(\d+)\,(\d+)\,\d+\,\d+/)
			{
				$Skip{"$1 $2"}++;
			}
		}
		if ((scalar(keys %Skip)>=3)||($t[2]>=$skip_s&&$t[3]<=$skip_e))
		{
			($skip_s,$skip_e)=($t[2],$t[3]) if (scalar(keys %Skip)>=3);
			if ((@sp_jun>0)&&($pre_chr eq $seqid)&&($pre_strand eq $strand))
			{
				foreach my $sp(@sp_jun)
				{
					if ((exists $jun_gene{$sp->[5]})&&($sp->[0]>=$skip_s)&&($sp->[1]<=$skip_e))
					{
						my $sp_juncid=(split /\:/,$sp->[4])[0];
						if (!exists $jun_gene{$sp->[5]}{$sp_juncid})
						{
							push @{$$hash{$sp->[5]}{$label}},[@$sp[0..4]];
							$jun_gene{$sp->[5]}{$sp_juncid}++;
						}
					}
				}
			}
			else
			{
				@sp_jun=();
			}
			for (my $i=8;$i<@t;$i++)
			{
				if ($t[$i]=~/ID\=([^\;]+)\;/)
				{
					if (exists $GeneLH{$1})
					{
						my $getGene=$1;
						#print STDERR "$getGene\n";
						if (!exists $jun_gene{$getGene}{$juncid})
						{
							push @{$$hash{$getGene}{$label}},[$t[2],$t[3],$len1,$len2,"$juncid: $num"];
							$jun_gene{$getGene}{$juncid}++;
						}
					}
				}
				if ($t[$i]=~/name\=([^\;]+)\;/)
				{
					if (exists $GeneLH{$1})
					{
						my $getGene=$1;
						#print STDERR "$getGene\n";
						if (!exists $jun_gene{$getGene}{$juncid})
						{
							push @{$$hash{$getGene}{$label}},[$t[2],$t[3],$len1,$len2,"$juncid: $num"];
							$jun_gene{$getGene}{$juncid}++;
						}
					}
				}
			}
		}
		else
		{
			($skip_s,$skip_e)=(0,0);
			for (my $i=8;$i<@t;$i++)
			{
				if ($t[$i]=~/ID\=([^\;]+)\;/)
				{
					if (exists $GeneLH{$1})
					{
						my $getGene=$1;
						push @sp_jun,[$t[2],$t[3],$len1,$len2,"$juncid: $num",$getGene] unless (exists $jun_gene{$getGene}{$juncid});
					}
				}
				if ($t[$i]=~/name\=([^\;]+)\;/)
				{
					if (exists $GeneLH{$1})
					{
						my $getGene=$1;
						push @sp_jun,[$t[2],$t[3],$len1,$len2,"$juncid: $num",$getGene] unless (exists $jun_gene{$getGene}{$juncid});
					}
				}
			}
		}
		($pre_chr,$pre_strand)=($seqid,$strand);
	}
	close IN;
}

## GFF format
## seqid  source  type  start  end  score  strand  phase  attributes
sub parseGff
{
	my ($file,$label,$hash)=@_;
	my $getGene="";
	print STDERR readTime()."parsing Gene Gff file: $file\n" if ($Verbose);
	open (IN,$file) || die $!;
	while(<IN>)
	{
		s/\s*$//;
		next if (($_=~/\#/)||($_ eq ""));
		my ($seqid,$source,$type,$start,$end,$score,$strand,$phas,$attributes)=split /\t+/,$_;
		next if ($strand !~ /[\+\-]/);
		if ((($type =~ /mRNA/i)||($type =~ /gene/i))&&($label=~/refGene/i))
		{
			$getGene="";
			while ($attributes=~/\=([^\;]+)/gi)
			{
				if ((exists $GeneLH{"$file:$1"})||(exists $GeneLH{"$1"}))
				{
					print STDERR "\tfind gene: $1\n" if ($Verbose);
					$getGene=$1;
					push @{$GeneRH{$seqid}},[$start,$end,$strand,$getGene];
					push @{$GeneRH{$getGene}},[$start,$end,$strand,$seqid];
					$GeneBH{$getGene}++;
				}
			}
		}
		elsif (($label!~/refGene/i)&&(($type=~/transcript/i)||($attributes=~/gene\=/i)||($attributes=~/name\=/i)))
		{
			$getGene="";
			if ($source eq "AS")
			{
				next if (defined $Only && !exists $OnlyDetect{$type});
			}
			if ($attributes=~/FPKM\s+\"([^\"\s]+)\"/)
			{
				next if ($1<$FPKM);
			}
			if ($type=~/transcript/i)
			{
				if (exists $GeneRH{$seqid})
				{
					my $skip=0;
					if ($attributes=~/gene\_id\s+\"([^\"\s]+)\"/ || $attributes=~/transcript\_id\s+\"([^\"\s]+)\"/)
					{
						if (exists $GeneBH{$1})
						{
							$getGene=$1;
							$GeneBH{$1}++;
							print STDERR readTime()."\t$label: get gene's transcript: $attributes\n" if ($Verbose);
							$skip=1;
						}
					}
					if ($skip==0)
					{
						foreach my $gene_region(@{$GeneRH{$seqid}})
						{
							if (($start<=$gene_region->[1])&&($end>=$gene_region->[0])&&($strand eq $gene_region->[2]))
							{
								$getGene=$gene_region->[3];
								print STDERR readTime()."\t$label: get gene's transcript: $attributes\n" if ($Verbose);
								$GeneBH{$getGene}++;
							}
						}
					}
				}
			}
			elsif ($source=~/AS/i)
			{
				while ($attributes=~/\=([^\;]+)\;/gi)
				{
					my $gene=$1;
					if ((exists $GeneLH{"$file:$gene"})||(exists $GeneLH{"$gene"}))
					{
						$getGene=$gene;
						$GeneBH{$getGene}++;
						print STDERR readTime()."\t$label: get gene's alternativew isform: $attributes\n" if ($Verbose);
					}
					while($gene=~/([^\,]+)\,/gi)
					{
						if ((exists $GeneLH{"$file:$1"})||(exists $GeneLH{"$1"}))
						{
							$getGene=$1;
							$GeneBH{$getGene}++;
							print STDERR readTime()."\t$label: get gene's alternativew isform: $attributes\n" if ($Verbose);
						}
					}
				}
			}
		}
		push @{$$hash{$getGene}{$label}},[$type,$start,$end,$strand,$attributes] if ($getGene ne "");
	}
	close IN;
}

sub parseList
{
	my ($list,$hash)=@_;
	print STDERR readTime()."parsing Gene list file: $list\n" if ($Verbose);
	open (IN,$list) || die $!;
	while(<IN>)
	{
		s/\s*$//;
		if (/([^\:]+)[\:^\s+](\S+)/)
		{
			$$hash{"$1:$2"}=1;
			print STDERR "\t$1:$2\n" if ($Verbose);
		}
		elsif (/(\S+)\s+(\S+)/)
		{
			$$hash{"$1:$2"}=1;
			print STDERR "\t$1:$2\n" if ($Verbose);
		}
		elsif (/^(\S+)$/)
		{
			$$hash{$1}=1;
			print STDERR "\t$1\n" if ($Verbose);
		}
		else
		{
			die "can't be recognized format of $list!\n";
		}
	}
	close IN;
}

sub parseOnlyDetectAStypes
{
	my ($only,$asdettype)=@_;
	my @types=split /\:/,$Only;
	map{$$asdettype{$_}=1;print STDERR readTime()."only detect AS type: $_\n" if ($Verbose);}@types;
}

sub setupcustomizeColor
{
	my ($file,$color)=@_;
	print STDERR readTime()."setup customize color scheme\n" if ($Verbose);
	my $label="";
	open (CS,$file) || die "Can't open $file for reading\n";
	while(<CS>)
	{
		chomp;
		next if ($_ eq "" || $_=~/^#/);
		if (/Label\=(\S+)/)
		{
			$label=$1;
		}
		elsif (/([^\=\s]\=(\S+))/)
		{
			my ($L,$R)=($1,$2);
			$$color{$label}{$L}=$R if ($label ne "");
			if ($L eq 'color')
			{
				print STDERR "\t$label color is set as $R\n" if ($Verbose);
			}
			else
			{
				
				print STDERR "\t$label $L color is set as $R\n" if ($Verbose);
			}
		}
	}
	close CS;
}

sub parseConfig
{
	my $file=shift;
	my %hash;
	open (CF,$file) || die "Can't open $file for reading!\n";
	while(<CF>)
	{
		$hash{$1}=$2 if (/([^\=\s]+)\=([^\=\s]+)/);
	}
	close CF;
	$Figure_resolution ||= $hash{'Figure_resolution'};
	$Vertical_skip ||= $hash{'Vertical_skip'};
	$Left_margin ||= $hash{'Left_margin'};
	$Right_margin ||= $hash{'Right_margin'};
	$Top_margin ||= $hash{'Top_margin'};
	$Bottom_margin ||= $hash{'Bottom_margin'};
	$Font_family ||= $hash{'Font_family'};
	$Font_size ||= $hash{'Font_size'};
	$Font_size_title ||= $hash{'Font_size_title'};
	$Font_size_mark ||= $hash{'Font_size_mark'};
	$GeneBlock_Height ||= $hash{'GeneBlock_Height'};
	$Coverage_figure_height ||= $hash{'Coverage_figure_height'};
}

sub plotRuler
{
	my %rulcfg = @_;
	print STDERR readTime()." plot ruler...\n" if ($Verbose);
	$rulcfg{scaletype} ||= "bp";
	$rulcfg{scaletypepos} ||= "right";
	$rulcfg{scalestart} ||= "auto";
	$rulcfg{rulerstyle} ||= "1";
	$rulcfg{font_family} ||= "ArialNarrow";
	$rulcfg{font_size} ||= 8;
	
	my $scale_size = 10;
	my $divid = 100;
	my $unit_start;

	my $bp_len = $rulcfg{bp_end} - $rulcfg{bp_start};
	my $X_len = $rulcfg{X_end} - $rulcfg{X_start};

	##caculate the length of smallest unit 
	my ($str,$str1,$str2,$unit);
	$str = $bp_len / $divid;
	$str = sprintf("%e",$str);
	if ($str=~/([\d\.]+)e([\+\-\d]+)/) {
		$str1 = $1;
		$str2 = $2;
	}
	$str1 = int ( $str1 + 0.5 );
	$unit = $str1 * 10 ** $str2;
	$unit = $rulcfg{bigscalesize} if(defined $rulcfg{bigscalesize});

	my $g = $rulcfg{svg}->group('id'=>times().rand());
	
	## draw the main axis
	$g->line('x1',$rulcfg{X_start},'y1',$rulcfg{Y},'x2',$rulcfg{X_end},'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1);
	return if($rulcfg{bp_end} - $rulcfg{bp_start}  == 0);

	##decide unit start
	if ($rulcfg{scalestart} eq "auto")
	{
		$unit_start = $rulcfg{bp_start} + ($unit - $rulcfg{bp_start} % $unit);
	}
	if ($rulcfg{scalestart} eq "force")
	{
		$unit_start = int($rulcfg{bp_start} / 10 + 0.5) * 10;
	}
	
	## draw small scale lines
	for (my $i=$unit_start; $i<=$rulcfg{bp_end}; $i+=$unit)
	{
		my $X = $rulcfg{X_start} + ($i - $rulcfg{bp_start}) / $bp_len * $X_len;
		$g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size/2,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '1');
		$g->line('x1',$X,'y1',$rulcfg{Y} + $scale_size/2,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '2');
		$g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size/2,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '3');
	}
	## draw big scale lines and text scales
	for (my $i=$unit_start; $i<=$rulcfg{bp_end}; $i+=$unit*10)
	{
		my $X = $rulcfg{X_start} + ($i - $rulcfg{bp_start}) / $bp_len * $X_len;
		$g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '1');
		$g->line('x1',$X,'y1',$rulcfg{Y} + $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1)  if ($rulcfg{rulerstyle} eq '2');
		$g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '3');
		my $disp_scale_text = $i / $rulcfg{bigscalesize}  if(defined $rulcfg{bigscalesize});
		$rulcfg{svg}->text('x',$X - textWidth($rulcfg{font_family},$rulcfg{font_size}/2,$disp_scale_text) / 2,'y',$rulcfg{Y}+textHeight($rulcfg{font_size}),'fill','#000000','-cdata',$disp_scale_text,'font-size',$rulcfg{font_size}, 'font-family',$rulcfg{font_family}) if ($rulcfg{rulerstyle} eq '1');
		$rulcfg{svg}->text('x',$X - textWidth($rulcfg{font_family},$rulcfg{font_size}/2,$disp_scale_text) / 2,'y',$rulcfg{Y}+$scale_size+textHeight($rulcfg{font_size})+2,'fill','#000000','-cdata',$disp_scale_text,'font-size',$rulcfg{font_size}, 'font-family',$rulcfg{font_family})  if ($rulcfg{rulerstyle} eq '2');
		$rulcfg{svg}->text('x',$X - textWidth($rulcfg{font_family},$rulcfg{font_size}/2,$disp_scale_text) / 2,'y',$rulcfg{Y}-$scale_size-textHeight($rulcfg{font_size})+6,'fill','#000000','-cdata',$disp_scale_text,'font-size',$rulcfg{font_size}, 'font-family',$rulcfg{font_family}) if ($rulcfg{rulerstyle} eq '3');
	}

}

sub readTime () {
	my  ($sec,$min,$hour,$mday,$mon,$year) = (localtime)[0..5];
	($sec,$min,$hour,$mday,$mon,$year) = (
		sprintf("%02d", $sec),
		sprintf("%02d", $min),
		sprintf("%02d", $hour),
		sprintf("%02d", $mday),
		sprintf("%02d", $mon + 1),
		$year + 1900
	);
	return "[ $year-$mon-$mday $hour:$min:$sec ] ";
}

sub instruction
{
	print qq(
# gene.list
DDX11L1

# color_scheme.configure
Label=refGene
gene=blue
exon=blue
Label=UF1
transcripts=cyan
alternative=cyan
junctions=cyan
depth=cyan
splicingreads=green:red
Label=BF1
transcripts=orange
alternative=orange
junctions=orange
depth=orange
splicingreads=violet:yellow

# refGene.gff
chr1    refGene mRNA    11874   14408   .       +       .       ID=NR_046018; name=DDX11L1;
chr1    refGene exon    11874   12227   .       +       .       Parent=NR_046018;
chr1    refGene exon    12613   12721   .       +       .       Parent=NR_046018;

# transcripts.gtf  (Cufflink)
n790	Cufflinks	transcript	60	874	1000	.	.	gene_id "CUFF.1"; transcript_id "CUFF.1.1"; FPKM "5.4575354634"; frac "1.000000"; conf_lo "3.394781"; conf_hi "7.520290"; cov "9.308451";
n790	Cufflinks	exon	60	874	1000	.	.	gene_id "CUFF.1"; transcript_id "CUFF.1.1"; exon_number "1"; FPKM "5.4575354634"; frac "1.000000"; conf_lo "3.394781"; conf_hi "7.520290"; cov "9.308451";

# AS file format
chr1    AS      RI      661139  665335  .       -       .       ID=RI1; mRNA=NR_028327; gene=LOC100133331;
chr1    AS      constitutive    661139  665184  .       -       .       Parent=RI1;
chr1    AS      alternative     665185  665277  .       -       .       Parent=RI1;
chr1    AS      constitutive    665278  665335  .       -       .       Parent=RI1;

# junctions.bed (TopHat)
track name=junctions description="TopHat junctions"
chr1	14739	15034	JUNC00000001	3	-	14739	15034	255,0,0	2	90,65	0,230
chr1	14985	15842	JUNC00000002	1	-	14985	15842	255,0,0	2	53,47	0,810
chr1	16956	17331	JUNC00000003	21	-	16956	17331	255,0,0	2	99,99	0,276
chr1	16984	17291	JUNC00000004	2	-	16984	17291	255,0,0	2	71,33	0,274
chr1	17312	17673	JUNC00000005	2	-	17312	17673	255,0,0	2	56,68	0,293
chr1	17689	17961	JUNC00000006	1	-	17689	17961	255,0,0	2	53,47	0,225

# bam format refers to SAMtools
);
}