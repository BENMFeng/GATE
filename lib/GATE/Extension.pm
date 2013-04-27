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

$VERSION = "0.1";

#@ISA=qw(GATE::Element GATE::Extension);

# DTD declarations handled in this modules
use constant ELEMENT => "ELEMENT";
use constant ATTLIST => "ATTLIST";
use constant NOTATION => "NOTATION";
use constant ENTITY => "ENTITY";

@TYPES=(ELEMENT,ATTLIST,NOTATION,ENTITY);
%TYPES=map { $_ => 1 } @TYPES;

#-----------------
sub new {
	return shift->SUPER::new(@_);
}

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

1;