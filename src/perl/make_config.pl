#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use FindBin qw($Bin $Script);

my ($path,$clone,$pattern,$sm,$index,$barcode,$Help);
my %opts;
GetOptions(\%opts,"path:s"=>\$path,"clone:s"=>\$clone,"pattern:s"=>\$pattern,"sm:s"=>\$sm,"index:s"=>\$index,"barcode:s"=>\$barcode,"help"=>\$Help);
die qq(perl $0 <-path input_data_dir> [-clone example_config] [-pattern "_R1_:_R2_"] [-sm sample_list.txt] [-index index.list] [-barcode C]\n) if (!defined $path || $Help);

$path=~s/\/$//;
$pattern ||= "_R1_:_R2_";
my ($R1,$R2)=split /[\:\-]/,$pattern;
my @SM=parselist($sm) if (defined $sm);
my @Index=parselist($index) if (defined $index);
my @Barcode=parselist($barcode) if (defined $barcode);

my @filenames;
my %get_file;
&recur_read_dir ($path, \@filenames);
my $prelb="";
my $sm_n=0;
my $index_n=0;
my $barcode_n=0;
my $mergePE=0;
if (defined $clone && -f $clone) {
	my $print=0;
	open(IN, $clone) or die $!;
	while (<IN>) {
		if (/\[(\S+)\]/) {
			if ($1 ne "LIB") {
				$print=1;
			} else {
				$print=0;
			}
		}
		$mergePE=1 if (/mergeOverlapPE/);
		print if ($print==1)
	}
	close IN;
}

foreach my $myfile(@filenames)
{
	if ($myfile=~/$R1/ && -f $myfile && ($myfile=~/fastq$/i || $myfile=~/fq$/i || $myfile=~/gz$/i) )
	{
		my $otherfile=$myfile;
		$otherfile=~s/$R1/$R2/;
		my $lb=$2 if ($myfile=~/$path([^\/].*\/|\/)([^\/\s]+)\//);
		if ($lb ne $prelb)
		{
			print "[LIB]\nLB=$lb\nID=$lb\n";
		}
		if (defined $sm)
		{
			print "SM=$SM[$sm_n]\n";
			$sm_n++ unless ($sm_n==@SM-1);
		}
		else
		{
			print "SM=$lb" if ($lb ne $prelb);
		}
		if (defined $index)
		{
			print "Index=$Index[$index_n]\n";
			$index_n++ unless ($index_n==@Index-1);
		}
		if (defined $barcode) {
			print "Barcode=$Barcode[$barcode_n]\n";
			$barcode_n++ unless ($barcode_n==@Barcode-1);
		}
		my $PL=checkPL($myfile);
		print $PL if (defined $PL);
		if (-f $otherfile)
		{
			print "MergePE=TRUE\n" if ($mergePE==1);
			print "fq1=$myfile\nfq2=$otherfile\n";
		}
		else
		{
			print "fq=$myfile\n";
		}
		$prelb=$lb;
	}
}

sub parselist {
	my $self=shift;
	my @ary=();
	if (-f $self) {
		open(IN,$self) or die $!;
		while (<IN>) {
			push @ary,(split /\s+/,$_);
		}
		close IN;
	} else {
		@ary=split/\:/,$self;
	}
	return @ary;
}

sub recur_read_dir {
	my ($path, $r_filename_list) =@_;
	my $h_dir;
	my @all_path=sort glob("$path*");
	foreach my $cur_path(@all_path)
	{
		opendir($h_dir, $cur_path) or die "serious dainbramage: $!";
		my @allfiles = grep { not /^\.{1,2}\z/ } readdir $h_dir;
		my $filename;
		for(my $i =0; $i< @allfiles; $i++ ) {
			$filename = $allfiles[$i];
			my $absolute_name = get_absolute_name($cur_path, $filename);
			if( -d $absolute_name ) {
				&recur_read_dir($absolute_name, $r_filename_list);
				push(@$r_filename_list, $absolute_name) if (!exists $get_file{$absolute_name});
				$get_file{$absolute_name}=1;
			}
			else {
				push(@$r_filename_list, $absolute_name) if (!exists $get_file{$absolute_name});
				$get_file{$absolute_name}=1;
			}
		}
		closedir $h_dir;
	}
}

sub get_absolute_name {
	my ($dir_name, $file_name) =@_;
	my $absolute_name =$dir_name."/".$file_name;
	return $absolute_name;
}

sub checkPL {
	my $file=shift;
	open (FILEHANDLE,$file) if ($file!~/gz$/);
	open(FILEHANDLE,"gzip -cd $file|") if ($file=~/gz$/);
	my $firstline=<FILEHANDLE>;
	close FILEHANDLE;
	if ($firstline=~/^\@H/) {
		return "PL=ILLUMINA\nPU=HISEQ\n";
	} elsif ($firstline=~/^\@M/) {
		return "PL=ILLUMINA\nPU=MISEQ\n";
	}
}
