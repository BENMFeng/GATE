#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my ($Ref,$k1,$k2,$J,$Interpolation,$Replace,$Ignore,$Help);
my %opts;
GetOptions(\%opts,"ref:s"=>\$Ref,"key:s"=>\$k1,"ink:s"=>\$k2,"j:s"=>\$J,"split:s"=>\$Interpolation, "e:s"=>\$Replace, "ignore"=>\$Ignore, "help"=>\$Help);
my $usage = qq(
join.pl -- join tables according to the same columns(hash tables)
Usage: perl $0 [option] <Input files> [> Output file]
Option: -ref <file>		reference file,if undefined, it will join all the input files together by hash table
        -key FIELD		reference hash key columns field of ref, delmit by ",-.;:"
        -ink FIELD		input hash key columns field of Input files, delmit by ",-.;:"; ":" is the first-degree layout, "all" is set for all columns expcet "ink" column
        -j FIELD		joined columns field of Input files, delmit by ",-.;:"; ":" is the first-degree layout, "all" is set for all columns expcet "ink" column
        -split <str> 		set the interpolation for splitting data columns, such as [s+,s,t+,t] and some special symbol, default: t (i.e. \@col=split/\\t/)
        -e <str>		replace missing input fields with a string, default: "-"
        -ignore			ignore differences in case when comparing fields
        -help			help information

Example: perl join.pl -ref ref.txt -key 1-3 -ink all:1,2,4-6 -j all:7-10 -split "s+" Input1.txt Input2.txt Input3.txt >Output.txt

Notice: "-" is for successive number, "4-6" mean 4,5,6
[Windows] It just need to run join.pl in ActivePerl with the files of Input1.txt Input2.txt..., Parameter.txt
in the same directory, and it would output a file named "Output.txt".

After v2.4 alpha it will append the redundantly matching contents

Author : BENM <BinxiaoFeng\@gmail.com>
Version: 2.5 alpha
Data: Dec 20th, 2011
\n);

$Interpolation ||= "t";
$Replace = "-" if (!defined $Replace);
my @Input=@ARGV;
my @file=glob("./*.txt");
my $n=0;
foreach (@file)
{
	$n++ if ($_ =~ /Input1\.txt/);
	$n++ if ($_ =~ /Input2\.txt/);
	$n++ if ($_ =~ /Parameter\.txt/);
}
if ($n==3)
{
	warn "there is Input1.txt, Input2.txt and Parameter.txt here\n";
	$Ref="Input1.txt";
	my @tmp=glob("Input*.txt");
	@Input=() if (@tmp>0);
	undef @Input;
	foreach (@tmp)
	{
		push @Input,$_ unless ($_=~/Input1\.txt/);
	}
	read_parameter("Parameter.txt",\$k1,\$k2,\$J);
}
die $usage if (((@ARGV==0)&&($n!=3))||($Help)||(!defined $k2)||(!defined $J));

$k1=~s/(\d+)/($1-1)/eg if (defined $k1);
$k2=~s/(\d+)/($1-1)/eg if (defined $k2);
$J=~s/(\d+)/($1-1)/eg if (defined $J);
my @A=check_dash(split /\,/,$k1) if (defined $k1);
my @B=split /\:/,$k2  if (defined $k2);
my @C=split /\:/,$J  if (defined $J);
my %jcol;

if ((@B==1)&&(@Input>1))
{
	my $b=$B[0];
	$b=~s/(\d+)/($1+1)/eg;
	warn ("for all the files, -ink parametr set as $b\n");
	for (my $i=0;$i<@Input;$i++)
	{
		@{$jcol{$i}}=check_dash(split /\,/,$B[0]);
	}
}
elsif (@B!=@Input)
{
	if (($B[0]=~/all/)&&(@B==2))
	{
		my @b=split /\,/,$B[1];
		@b=check_dash(@b);
		for (my $i=0;$i<@Input;$i++)
		{
			@{$jcol{$i}}=@b;
		}
	}
	else
	{
		die "-ink parameter: $k2 not equal to the number of files!";
	}
}
else
{
	my $i=0;
	map{@{$jcol{$i++}}=check_dash(split /\,/,$_)}@B;
}

my %join;
my $All=0;
if ((@C==1)&&(@Input>1))
{
	if ($C[0] ne "all")
	{
		my $c=$C[0];
		$c=~s/(\d+)/($1+1)/eg;
		warn ("for all the files, -j parameter set as $c\n");
	}
	else
	{
		warn "for all the files, -j parameter set as $C[0]\n";
	}
	for (my $i=0;$i<@Input;$i++)
	{
		@{$join{$i}}=check_dash(split /\,/,$C[0]);
	}
}
elsif (@C!=@Input)
{
	if (($C[0]=~/all/)&&(@C==2))
	{
		my @c=split /\,/,$C[1];
		@c=check_dash(@c);
		for (my $i=0;$i<@Input;$i++)
		{
			@{$join{$i}}=@c;
		}
	}
	else
	{
		die "-j parameter: $J not equal to the number of files!";
	}
}
else
{
	my $i=0;
	map{@{$join{$i++}}=check_dash(split /\,/,$_)}@C;
}

for (my $i=0;$i<@Input;$i++)
{
	open (CHK,$Input[$i]) || die "Can open $Input[$i] for reading!\n";
	my $line=<CHK>;
	chomp $line;
	my @t=();
	if ($Interpolation eq "t")
	{
		@t=split /\t/,$line;
	}
	elsif ($Interpolation eq "t+")
	{
		@t=split /\t+/,$line;
	}
	elsif ($Interpolation eq "s")
	{
		@t=split /\s/,$line;
	}
	elsif ($Interpolation eq "s+")
	{
		@t=split /\s+/,$line;
	}
	else
	{
		@t=split/$Interpolation/,$_;
	}
	for (my $m=0;$m<@{$join{$i}};$m++)
	{
		if (${$join{$i}}[$m] eq "all")
		{
			my @tmp=@t;
			delete @tmp[@{$jcol{$i}}];
			@{$join{$i}}=();
			undef @{$join{$i}};
			my $e=0;
			for (my $d=0;$d<@tmp;$d++)
			{
				${$join{$i}}[$e++]=$d if (defined $tmp[$d]);
			}
			last;
		}
	}
	close CHK;
}

my %table;

for (my $i=0;$i<@Input;$i++)
{
	open (IN,$Input[$i]) || die "Can open $Input[$i] for reading!\n";
	while (<IN>)
	{
		s/\s+$//;
		my @col=();
		if ($Interpolation eq "t")
		{
			@col=split /\t/,$_;
		}
		elsif ($Interpolation eq "t+")
		{
			@col=split /\t+/,$_;
		}
		elsif ($Interpolation eq "s")
		{
			@col=split /\s/,$_;
		}
		elsif ($Interpolation eq "s+")
		{
			@col=split /\s+/,$_;
		}
		else
		{
			@col=split/$Interpolation/,$_;
		}
		next if (($_ eq "")||(@col==0));
		my $key=join "\t",@col[@{$jcol{$i}}];
		next if (!defined $key);
		if (exists $table{$key}{$i})
		{
			$table{$key}{$i} .= " ";
			$table{$key}{$i} .= (join "\t",@col[@{$join{$i}}]);
		}
		else
		{
			$table{$key}{$i} = (join "\t",@col[@{$join{$i}}]);
		}
	}
	close IN;
}

if (defined $Ref)
{
	open (REF,$Ref) || die "Can open $Ref for reading!\n";
	die $usage if (!defined $k1);
	open (OUT,">Output.txt") || die "Can write to Output.txt\n" if ($n==3);
	while (<REF>)
	{
		s/\s+$//;
		my @col=();
		if ($Interpolation eq "t")
		{
			@col=split /\t/,$_;
		}
		elsif ($Interpolation eq "t+")
		{
			@col=split /\t+/,$_;
		}
		elsif ($Interpolation eq "s")
		{
			@col=split /\s/,$_;
		}
		elsif ($Interpolation eq "s+")
		{
			@col=split /\s+/,$_;
		}
		else
		{
			@col=split/$Interpolation/,$_;
		}
		my $key=join "\t",@col[@A];
		my $out="";
		if (defined $Ignore)
		{
			my $I=0;
			for (my $i=0;$i<@Input;$i++)
			{
				if (exists $table{$key}{$i})
				{
					$out .= $table{$key}{$i}."\t";
				}
				else
				{
					$I=1;
					last;
				}
			}
			next if ($I==1);
		}
		else
		{
			for (my $i=0;$i<@Input;$i++)
			{
				if (exists $table{$key}{$i})
				{
					#$out .= " ";
					$out .= $table{$key}{$i}."\t";
				}
				else
				{
					$out .= ("$Replace\t"x(@{$join{$i}}))
				}
			}
		}
		$out=~s/\s+$//;
		if ($n==3)
		{
			print OUT "$_\t$out\n";
		}
		else
		{
			print "$_\t$out\n";
		}
	}
	close REF;
	close OUT if ($n==3);
}
else
{
	foreach my $key(keys %table)
	{
		my $out;
		if (defined $Ignore)
		{
			my $I=0;
			for (my $i=0;$i<@Input;$i++)
			{
				if (exists $table{$key}{$i})
				{
					$out .= $table{$key}{$i}."\t";
				}
				else
				{
					$I=1;
					last;
				}
			}
			next if ($I==1);
		}
		else
		{
			for (my $i=0;$i<@Input;$i++)
			{
				if (exists $table{$key}{$i})
				{
					#$out .= " ";
					$out .= $table{$key}{$i}."\t";
				}
				else
				{
					$out .= ("$Replace\t"x(@{$join{$i}}))
				}
			}
		}
		$out=~s/\s+$//;
		print "$key\t$out\n";
	}
}

sub read_parameter
{
	my ($file,$F1,$F2,$D)=@_;
	open (P,$file) || die "Can open Parameter.txt for reading!\n";
	while(<P>)
	{
		s/\s+$//;
		$$F1=$1 if ($_=~/f1\=(\S+)/);
		$$F2=$1 if ($_=~/f2\=(\S+)/);
		$$D=$1 if ($_=~/d\=(\S+)/);
	}
	close P;
	die "Parameter file is wrong! For example:
f1=1,2		#column number of Input1.txt
f2=3,4		#column number of Input2.txt
d=5,6,7		#merge columns of Input2.txt into the front of Input1.txt
the number in right of the equal sign means the column number of the input file
" if ((!defined $F1)||(!defined $F2)||(!defined $D));
}

sub check_dash
{
	my @array=@_;
	for (my $i=0;$i<@array;$i++)
	{
		if ($array[$i]=~/(\d+)\-(\d+)/)
		{
			my ($s,$e)=($1<$2)?($1,$2):($2,$1);
			my @tmp=();
			for (my $j=$s;$j<=$e;$j++)
			{
				push @tmp,$j;
			}
			splice(@array,$i,1,@tmp);
		}
	}
	return @array;
}

__END__
