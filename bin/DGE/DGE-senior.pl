#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $program_name = $1 if ( $0 =~ /([^\/]+)$/);
my $usage = <<USAGE;
$program_name -- Gene Expression Clustering

Author: BENM <binxiaofeng\@gmail.com>
Version: 0.1 alpha
Date: 2010-07-08

$program_name [options] [functions]
	matrix <file>		input Matrix data file
	maxEd <int>		set limited Euclid distance
	help			help & usage

Example:

USAGE

my ($Matrix_file,$Normal,$MaxEd,$Help);
my %opts;
GetOptions
(
	\%opts,
	"matrix:s"=>\$Matrix_file,
	"maxEd:i"=>\$MaxEd,
	"help"=>\$Help
);
die $usage if (($Help)||(!defined $Matrix_file));
$MaxEd ||= 10;

my @Matrix=();
my @Tot=();
my %Euclid_Triangle=();
my %Euclid_Array=();
&read_matrix($Matrix_file);
&cal_Ed if (@Matrix>=2);
#print Dumper %Euclid_Array;
my $cutoff=train_cutoff(keys %Euclid_Array);

my %Cluster;
my %Clusterred;
my $c=0;
foreach my $ed_value(sort{$a<=>$b}keys %Euclid_Array)
{
	last if ($ed_value>$cutoff);
	foreach my $m (sort{$a->[0]<=>$b->[0]} @{$Euclid_Array{$ed_value}})
	{
		my ($i,$j)=@$m;
		next if ((exists $Clusterred{$i})||(exists $Clusterred{$j}));
		last if ($Euclid_Triangle{$i}{$j}>$cutoff);
		
		$Cluster{$c}{$i}=$Euclid_Triangle{$i}{$j};
		$Cluster{$c}{$j}=$Euclid_Triangle{$i}{$j};
		
		$Clusterred{$i}=1;
		$Clusterred{$j}=1;
		
		my $x;
		foreach my $k (sort{$Euclid_Triangle{$i}{$a}<=>$Euclid_Triangle{$i}{$b}}keys %{$Euclid_Triangle{$i}})
		{
			next if (exists $Cluster{$c}{$k});
			$x=$k if ($Euclid_Triangle{$i}{$k}<=$cutoff);
			last;
		}
		my $y;
		foreach my $k (sort{$Euclid_Triangle{$j}{$a}<=>$Euclid_Triangle{$j}{$b}}keys %{$Euclid_Triangle{$j}})
		{
			next if (exists $Cluster{$c}{$k});
			$y=$k if ($Euclid_Triangle{$j}{$k}<=$cutoff);
			last;
		}
		
		my ($I,$Y);
		if (defined $y)
		{
			my $score=0;
			foreach my $A(keys %{$Cluster{$c}})
			{
				($I,$Y)=($A<$y)?($A,$y):($y,$A);
				$score++ if ((exists $Euclid_Triangle{$I}{$Y})&&($Euclid_Triangle{$I}{$Y}<=$cutoff));
			}
			if ($score==(keys %{$Cluster{$c}}))
			{
				$Cluster{$c}{$y}=$Euclid_Triangle{$I}{$Y};
				$Clusterred{$y}=1;
			}
		}
		my ($J,$X);
		if (defined $x)
		{
			my $score=0;
			foreach my $B(keys %{$Cluster{$c}})
			{
				($J,$X)=($B<$x)?($B,$x):($x,$B);
				$score ++ if ((exists $Euclid_Triangle{$J}{$X})&&($Euclid_Triangle{$J}{$X}<=$cutoff));
			}
			if ($score==(keys %{$Cluster{$c}}))
			{
				$Cluster{$c}{$x}=$Euclid_Triangle{$J}{$X};
				$Clusterred{$x}=1;
			}
		}
		
		if ((!defined $x)&&(!defined $y))
		{
			$c++;
			next;
		}
	}
}

foreach my $C(sort{$a<=>$b} keys %Cluster)
{
	my $out=join "-",(sort{$a<=>$b}keys %{$Cluster{$C}});
	print "Cluster-$C:\t$out\n";
}

#######################################################################################
#---------------------------------- Main  Function ----------------------------------#
######################################################################################

sub read_matrix
{
	my $file=shift;
	open (IN,$file) || die "Can't loading this matrix file for reading\n";
	while(<IN>)
	{
		chomp;
		my @t=split /\t+/,$_;
		next if (($t[0] eq "ID")||($t[0] =~ /Tag/));
		my @copynumber=();
		for (my $m=1; $m<@t; $m++)
		{
			if (($t[$m]=~/^\d+[\.]\d+$/)||($t[$m]=~/^\d+$/))
			{
				$Tot[$m-1]+=$t[$m];
				$copynumber[$m-1]=$t[$m];
			}
			elsif ($t[$m]=~/(\d+)/)
			{
				$Tot[$m-1]+=$1;
				$copynumber[$m-1]=$1;
			}
			elsif ($t[$m]=~/\((\d+)[\,\;]/)
			{
				$Tot[$m-1]+=$1;
				$copynumber[$m-1]=$1;
			}
			else
			{
				$Tot[$m-1]+=0;
				$copynumber[$m-1]=0;
			}
		}
		push @Matrix,[@copynumber];
	}
	close IN;
}

sub cal_Ed
{
	for (my $i=0;$i<@Matrix-1;$i++)
	{
		my $n=@{$Matrix[$i]};
		for (my $j=$i+1;$j<@Matrix;$j++)
		{
			my $tot_square=0;
			for (my $k=0;$k<$n;$k++)
			{
				$tot_square+=($Matrix[$i][$k]-$Matrix[$j][$k])**2;
			}
			my $Ed_value=int(100000*sqrt($tot_square/(@{$Matrix[$j]}-1)))/100000;
			$Euclid_Triangle{$i}{$j}=$Ed_value if ($Ed_value<=$MaxEd);
			push @{$Euclid_Array{$Ed_value}},[$i,$j];
		}
	}
}

sub train_cutoff
{
	my @array=sort{$b<=>$a}@_;
	my $cut_off=0;
	my ($n,$total,$median,$mean,$square,$sd)=(0,0,0,0,0,0);
	$n=@array;
	for (my $i=0; $i<$n; $i++)
	{
		$total+=$array[$i];
		$square+=$array[$i]**2;
	}
	$median=$array[int($n/2)];
	$mean=$total/$n;
	$sd=sqrt($square-$total**2/$n)/($n-1);
	if (($mean<$MaxEd)&&($median<$MaxEd))
	{
		if (($mean>$median-$sd)&&($mean<$median+$sd))
		{
			$cut_off=$mean+$sd;
		}
		elsif ($mean<$median-$sd)
		{
			$cut_off=$median-$sd;
		}
		elsif ($mean>$median-$sd)
		{
			$cut_off=$median+$sd;
		}
	}
	$cut_off=(($cut_off>$MaxEd)||($cut_off==0))?$MaxEd:$cut_off;
	return $cut_off;
}

