=pod 

=head1 NAME

GATE - Perl extension for Genomics Applications, also for Transcriptomics & Epgenetics analysis pipeline

=cut

package GATE;

use strict;
use vars qw($VERSION @ISA $AUTOLOAD);
use Exporter;
use FindBin qw($Bin $Script);
use lib "$FindBin::Bin/../lib";
use File::Basename qw(basename dirname);
use GATE::Element;
#use GATE::Extension;
@ISA = qw(GATE::Element GATE::DO);
#@ISA = qw(GATE::Element GATE::Extension GATE::DO);


$VERSION = "1.2b";

#-------------------------------------------------------------------------------

=pod 

=head2 VERSION

Version 1.0a, 12 September, 2012

Refer to L<RST::Manual> for the complete manual

=head1 DESCRIPTION


Refer to L<RST::Manual> for the complete manual

=head1 AUTHOR

BENM(Binxiao) Feng <BinxiaoFeng@gmail.com>

=head1 CREDITS


=head1 EXAMPLES


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


#-------------------------------------------------------------------------------

my %default_attrs = (
    # processing options
    -config     => 'config.example',
    -workdir    => "./",
    -auto       => 0,       # permit arbitrary autoloads (only at import)
    -verbose    => 1,

    # default options
    "CustomSetting:MergeSamFiles" => 'USE_THREADING=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT',
    "CustomSetting:bwaaln"  => '-t 12 -n 2 -q 10',
    "CustomSetting:tophat"  => '-p 12 -r 200',
    "CustomSetting:cufflinks" => '-p 12 -u',
    "CustomSetting:cuffmerge" => '-p 12',
    "CustomSetting:cuffcompare" => '-R -T -V -CG',
    "CustomSetting:cuffdiff" => '-p 12 -u --geometric-norm',
    "CustomSetting:PATH" => '/usr/local/bin',
    "CustomSetting:PYTHONPATH" => '/usr/local/lib/python2.7/site-packages',
    "CustomSetting:heap" => "4g",
    "CustomSetting:multithreads" => 4,
    "CustomSetting:qsub" => "-cwd -l cpu=12,vf=16G",
    "CustomSetting:mergeOverlapPE" => "-seed 11",
    "CustomSetting:filterPCRdup" => "-cmp 100",
    "CustomSetting:scanAP" => "-k 10 -e 3 -m 10 -l",
    "CustomSetting:fastqcut" => "-q 20 -f L -n 5 -len 35 -sd 10 -trim",
    "CustomSetting:trim_seq" => "-edge 3 -len_t 25 -len_p 0.25 --trim_mode both --ap_mode -verbose",
    "CustomSetting:qc_outdir" => "qc",
    "CustomSetting:SolexaQA" => "-v -m -s 10000 -b -sanger",

    "CustomSetting:aln_outdir" => "aln",
    "CustomSetting:var_outdir" => "var",
    "CustomSetting:as_outdir" => "as",

    "CustomSetting:PBS-header" => qq(#!/bin/bash
#PBS -N  mpi_template
#PBS -l nodes=1:ppn=10
#PBS -l walltime=01:40:00
#PBS -j oe
#PBS -q mem.q
NSLOTS=`cat \${PBS_NODEFILE} | wc -l`
cd  \$PBS_O_WORKDIR\n\n),
);

#sub import {
#	my $package=shift;
#
#	my $attr=undef;
#	foreach (@_) {
#	if ($attr) {
#		$default_attrs{$attr}=$_;
#		undef $attr;
#	} elsif (exists $default_attrs{$_}) {
#		$attr=$_;
#	} else {
#		/^-/ and die "Unknown attribute '$_' in import list\n";
#		$RST::Element::autosubs{$_}=1; # add to list of autoloadable tags
#	}
#	}
#	return ();
#}

sub new ($;@) {
	my ($proto,%attrs)=@_;
	my $class=ref $proto || $proto;
	my $self;
	
	# establish defaults for unspecified attributes
	foreach my $attr (keys %default_attrs) {
		$attrs{$attr}=$default_attrs{$attr} unless exists $attrs{$attr}
	}

	$self->{$_} = $attrs{$_} foreach keys %default_attrs;
	bless $self;
	return $self;
}

sub set_attributes { 
	my $object = shift; 
	my $attribute_name; 
	if ( ref ($_[0] ) ) { 
		my ($attribute_name_list, $attribute_value_list) = @_; 
		my $i = 0; 
		foreach $attribute_name (@{$attribute_name_list}) { 
			my $set_method_name = "set_" . $attribute_name; 
			$object->$set_method_name ($attribute_value_list->[$i++]); 
		} 
	} else { 
		my ($attribute_name, $attribute_value); 
		while (@_) { 
			$attribute_name = shift; 
			$attribute_value = shift; 
			my $set_method_name = "set_" . $attribute_name; 
			$object->$set_method_name ($attribute_value); 
		} 
	} 
} 

sub get_attributes { 
	my $object = shift; 
	my (@retval); 
	foreach my $attribute_name (@_) { 
		my $get_method_name = "get_" . $attribute_name; 
		push ( @retval, $object->$get_method_name() ); 
	} 
	return @retval; 
}

sub get_attribute_names {
	my $package = shift;
	if ( ref ($package) ) {
		$package = ref ($package); 
	}
	my @result = @{"${package}::_ATTRIBUTES_"};
	if ( @{"${package}::ISA"} > 0 ) {
		foreach my $base_package (@{"${package}::ISA"}) {
			push ( @result, get_attribute_names ($base_package) );
		}
	}
	return @result;
}

1;