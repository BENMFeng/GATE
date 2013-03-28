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


1;