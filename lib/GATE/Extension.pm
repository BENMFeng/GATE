=pod 

=head1 NAME

GATE::Extension - Generate the extension bits for GATE.pm

=head1 AUTHOR

BENM(Binxiao) Feng, binxiaofeng@gmail.com

=head1 SEE ALSO
##########################################################################
#  Copyright (c) 2012 - 2013 - BENM(Binxiao) Feng                        #
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


=cut

package GATE::Extension;
use strict;
use FindBin qw($Bin $Script);
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../../lib";
use File::Basename qw(basename dirname);
use GATE::DO;

use vars qw(@ISA $VERSION @TYPES %TYPES);

$VERSION = "0.3, 2013-09-12";

#@ISA=qw(GATE::Element GATE::Extension);

# DTD declarations handled in this modules
use constant ELEMENT => "ELEMENT";
use constant ATTLIST => "ATTLIST";
use constant NOTATION => "NOTATION";
use constant ENTITY => "ENTITY";

@TYPES=(ELEMENT,ATTLIST,NOTATION,ENTITY);
%TYPES=map { $_ => 1 } @TYPES;

#-----------------

sub make_soapdenovo_config ($) {
	my $self = shift;
}

sub stat_reads ($) {
	my $self = shift;
	my @Reads=@{$self->{'reads'}};
	foreach my $reads(@Reads) {
		my $type=check_format($reads);
		if ($type==0)
		{
			
		} elsif ($type==1)
		{
			
		}
	}
}

sub stat_mappedreads ($) {
	my ($self,%attrs) = @_;
	my $cmd="";
	my $samtools = checkPath($self->{"software:samtools"});
	my $sam2bed = checkPath($self->{"software:sam2bed"});
	my $overlap = checkPath($self->{"software:overlap"});
	my $msort = checkPath($self->{"software:msor"});
	$cmd .= qq(alias sortbychr="$msort -k '1{1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM}'"\n);
	if (exists $attrs{'bam'}) {
		my $bam=$attrs{'bam'};
		my $lib=$attrs{'lib'};
		my $prefix=$1 if ($bam=~/(\S+)\.bai/);
		if (!-f "$bam.bai" && !-f "$prefix.bai") {
			$cmd.=qq($samtools $bam\n);
		}
		$cmd .= qq($samtools -F 4 $bam | $sam2bed -n > $lib.bed\n);
		$cmd .= qq($samtools -F 4 $bam | $sam2bed -u -n > $lib.uniq.bed\n);
		$cmd .= qq(cut -f 1,4 $lib.bed |sort -u |cut -f 1 | \${sortbychr} |uniq -c |awk '{print \$2"\\t"\$1}' > $lib.chr.mappedreads.txt\n);
		$cmd .= qq(cut -f 1,4 $lib.uniq.bed |sort -u |cut -f 1 | \${sortbychr} |uniq -c |awk '{print \$2"\\t"\$1}' > $lib.uniq.chr.mappedreads.txt\n);
		$cmd .= qq($overlap --C --i1 $lib.bed --f1 0-1-2 |perl -e 'my \%hash;<>;while(<>){my \@t=split;\$hash{\$t[0]}+=\$t[3]}foreach my \$chr(keys %hash){print "\$chr\\t\$hash{\$chr}\\n"}' |\${sortbychr} > $lib.chr.cov-len.txt\n);
		$cmd .= qq($overlap --C --i1 $lib.uniq.bed --f1 0-1-2 |perl -e 'my \%hash;<>;while(<>){my \@t=split;\$hash{\$t[0]}+=\$t[3]}foreach my \$chr(keys %hash){print "\$chr\\t\$hash{\$chr}\\n"}' |\${sortbychr} > $lib.uniq.chr.cov-len.txt\n);
		$cmd .= qq($overlap --C --i1 $lib.bed --f1 0-3-1-2 |perl -e 'my \%hash;<>;while(<>){my \@t=split;my \$depth=0;foreach my \$i\(5..\(\@t-1\)\){\$depth+=\$1 if (\$t[\$i]=~\/\\,\(\\d+\)\$\/);}\$hash{\$t[0]}{0}+=\$depth\/\$t[3];\$hash{\$t[0]}{1}++;}foreach my \$chr(keys \%hash){print \("\$chr\\t",\$hash{\$chr}{0}/\$hash{\$chr}{1},"\\n"\)}'|\${sortbychr} > $lib.cov-depth.txt\n);
		$cmd .= qq($overlap --C --i1 $lib.uniq.bed --f1 0-3-1-2 |perl -e 'my \%hash;<>;while(<>){my \@t=split;my \$depth=0;foreach my \$i\(5..\(\@t-1\)\){\$depth+=\$1 if (\$t[\$i]=~\/\\,\(\\d+\)\$\/);}\$hash{\$t[0]}{0}+=\$depth\/\$t[3];\$hash{\$t[0]}{1}++;}foreach my \$chr(keys \%hash){print \("\$chr\\t",\$hash{\$chr}{0}/\$hash{\$chr}{1},"\\n"\)}'|\${sortbychr} > $lib.uniq.cov-depth.txt\n);
	}
	return $cmd;
}

sub stat_aln ($) {
	my $self = shift;
}

sub check_format ($) {
	my $self=shift;
	if (-B $self) {
		if ($self=~/bam$/i)
		{
			my $index=checkIndex('samtools',$self);
			if ($index==1)
			{
				return "bam";
			}
		} elsif ($self=~/fasta.gz$/i || $self=~/fa.gz$/i) {
			return "fasta.gz";
		} elsif ($self=~/fastq.gz$/i || $self=~/fq.gz$/i) {
			return "fastq.gz";
		} elsif ($self=~/gz$/i)
		{
			open (IN,"zcat $self|");
			my $line=<IN>;
			if ($line=~/^\>/) {
				return "fasta.gz";
			} elsif ($line=~/^\@/) {
				return "fastq.gz";
			}
			close IN;
		}
	} elsif (-T $self) {
		if ($self=~/sam/i) {
			return "sam";
		} elsif ($self=~/fasta$/i || $self=~/fa$/i || $self=~/seq/i || $self=~/fna/i || $self=~/fas/i) {
			return "fasta";
		} elsif ($self=~/fastq$/i || $self=~/fq/i) {
			return "fastq";
		} else {
			open (IN,$self);
			my $line=<IN>;
			if ($line=~/^\>/) {
				return "fasta";
			} elsif ($line=~/^\@/) {
				return "fastq";
			}
			close IN;
		}
	}
}

sub check_version ($)
{
	my $selft=shift;
}

sub plottingLociVsCoverage ($)
{
	my ($infile,$outfile,$title);
	my $R_cmd=qq(
##Author: Johan Dahlberg
library(ggplot2)

#Function to put the data in the correct format for plotting
preprocessCoverageFile <- function\(file, timeSeries\) {
  # Load coverage data
  data = read.delim\(file, header=T\)
  # Remove everything before ":"
  data\$Locus<-gsub\(".*:","",data\$Locus\)
  data\$Locus<-as.numeric\(data\$Locus\)
  data<-subset\(data, select=c\(Locus, Total_Depth\)\)
  data\$Time<-c\(as.character\(timeSeries\)\)
  data
}

#Function to plot the coverage
plotCoverage <- function\(coverage, title\) {
  ggplot\(data=coverage, aes\(x=Locus, y=Total_Depth, colour=Time\)\) +
    geom_line\(\) +
    ylab\("Depth of coverage"\) +
    xlab\("Loci"\) +
    labs\(title=title\)
}


# Example of using the above functions.

outputPDFFile = "$outfile"
# Load coverage data
coverageSample = preprocessCoverageFile($infile, $title)

# Create dataset with all data together.

# Plot the data to a pdf
pdf\(outputPDFFile, paper="a4r", width=10\)

plotCoverage\(coverageSample, "Coverage plot for $title"\)

dev.off\(\)
);
}

sub plottingCumulativeCoverage($){
	my $R_cmd=qq(
##Author: Johan Dahlberg
library\(ggplot2\)

# Load coverage data \(Note that I have removed all text from the original output files from gatk\)
lib1 = read.delim\("<file1>", header=F\)
lib2 = read.delim\("<file2>", header=F\)

# Remove extra line
lib1<-lib1[,2:length\(lib1\)]
lib2<-lib2[,2:length\(lib2\)]

# Combining the data sets
data<-data.frame\(t\(lib1\),rep\("manual"\)\)
colnames\(data\)<-c\("coverage","loci","prep_type"\)
lib2<-cbind\(t\(lib2\),rep\("robot"\)\)
colnames\(lib2\)<-c\("coverage","loci","prep_type"\)
data<-rbind\(data, lib2\)

# Remove the "gte_" part from the coverage function
data\$coverage<-sapply\(data\$coverage, function\(x\) sub\("gte_", "",x\)\)

# Convert to reasonable data types...
data\$coverage<-as.numeric\(data\$coverage\)
data\$loci<-as.numeric\(as.character\(data\$loci\)\)

sumManualLociCovered <- sum\(as.numeric\(as.character\(data[data\$prep_type == "manual",]\$loci\)\)\)
sumRobotLociCovered <- sum\(as.numeric\(as.character\(data[data\$prep_type == "robot",]\$loci\)\)\)


totalNumberOfLociCovered <- data\$loci[1]

data\$loci <- data\$loci / totalNumberOfLociCovered

# Normalized over the total number of reads in each file.
data[data\$prep_type == "manual",]\$coverage <-data[data\$prep_type == "manual",]\$coverage/764574276
data[data\$prep_type == "robot",]\$coverage <-data[data\$prep_type == "robot",]\$coverage/208914034

ggplot\(data=data, aes\(x=log\(coverage\), y=loci, color=prep_type\)\) +
  geom_point\(\) +
  ylab\("% of bases covered"\) +
  xlab\("log\(1/base\)"\) +
  opts\(title = "Cumulative coverage compassion"\)
);
}


sub check_run ($;@) {
	my $self = shift;
	my ($shell,$prefix)=@_;
	$prefix ||= "new";
	my $i=0;
	my @filehandle=();
	my @outhandle=();
	parse_sh($shell,\@filehandle,\@outhandle,$i,$prefix);
}

sub parse_sh {
	my ($shfile,$inhandle,$outhandle,$i,$prefix)=@_;
	my $glob;
	my $skip=0;
	my @path=();
	my %hash=();
	my $in=$$inhandle[$i];
	my $out=$$outhandle[$i];
	open($in,$shfile) || die "Can't open $shfile for reading!\n";
	open($out, ">$prefix.$shfile") || die "Can't write to new.$shfile\n";
	while (<$in>) {
		chomp;
		my $line=$_;
		next if (!defined $_ || $_ eq "" || $_=~/^\s+$/);
		while ($line=~/([^\"\s]+\.sh)/g) {
			$i++;
			parse_sh($1,$inhandle,$outhandle,$i,$prefix);
		}
		if ($line=~/export\s+([^\=\s]+)\=(\S+)/) {
			my ($alias,$source)=($1,$2);
			$source=~s/\"//g;
			$hash{$alias}=$source;
		}
		if ($line=~/cd\s+(\S+)/) {
			my $dir=$1;
			$dir=~s/\"//g;
			while ($dir =~ /\.\./g) {
				pop @path;
			}
			if ($dir ne ".." && $dir !~ /^\//) {
				if ($dir=~/\$\{([^\{\}\s]+)\}/ || $dir=~/\$(\S+)/) {
					my $alias=$1;
					if (exists $hash{$alias}) {
						$dir=$hash{$alias};
					} else {
						die "$shfile:\n\t$line\n";
					}
				}
				push @path,$dir;
			} elsif ($dir=~/^\//) {
				@path=();
				push @path,$dir;
			}
			my $fullpath=join "/",@path;
			my @finished=glob("$fullpath/finished*");
			if (@finished>0) {
				$skip=1;
			} else {
				$skip=0;
			}
		}
		if ($skip==0 || $line=~/^cd\s+\S+$/) {
			print $out "$line\n";
		} else {
			print $out "\#$line\n";
		}
	}
	close $in;
	close $out;
}
1;