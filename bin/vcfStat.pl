#!/usr/bin/perl -w
use strict;

my @vcf=glob("*.vcf");

my @sub=("G->A","A->G","T->C","C->T","C->A","T->A","C->G","T->G","A->C","G->C","A->T","G->T");
my %Snp;
my %Indel;
my %hash;
my %SNPtype=();
my %Chr=();
my @chr;
foreach my $f(@vcf)
{
	my $id=$1 if ($f=~/([^\.]+)\.gatk/);
	open (IN,$f) || die $!;
	my ($snp,$indel) = (0,0);
	while(<IN>)
	{
		next if (/^\#/);
		my @t=split;
		if ($_=~/INDEL/)
		{
			$Indel{"$t[0]-$t[1]-$t[3]"}{$id}=$t[4];
			$indel++;
		}
		else
		{
			$Snp{"$t[0]-$t[1]-$t[3]"}{$id}=$t[4];
			$snp++;
			$SNPtype{"$t[3]\-\>$t[4]"}{$id}++;
			push @chr,$t[0] if (!exists $Chr{$t[0]});
			$Chr{$t[0]}++;
		}
	}
	if ($snp>0)
	{
		$hash{$id}{'snp'}=$snp;
	}
	else
	{
		$hash{$id}{'snp'}=0;
	}
	if ($indel>0)
	{
		$hash{$id}{'indel'}=$indel;
	}
	else
	{
		$hash{$id}{'indel'}=0;
	}
	close IN;
}

my @SM=sort keys %hash;
print ("VAR\t",(join "\t",@SM),"\n");
print "snp";
foreach my $sm(@SM)
{
	print ("\t",$hash{$sm}{'snp'});
}
print "\n";
print "indel";
foreach my $sm(@SM)
{
	print ("\t",$hash{$sm}{'indel'});
}
print "\n\n";
print ("#substitutions\t",(join "\t",@SM),"\n");
foreach my $type(@sub)
{
	print $type;
	foreach my $sm(@SM)
	{
		if (exists $SNPtype{$type}{$sm})
		{
			print ("\t",$SNPtype{$type}{$sm});
		}
		else
		{
			print "\t0";
		}
	}
	print "\n";
}
print "\n";
print ("#CHR\tSNP\n");
foreach my $c(@chr)
{
	print "$c\t$Chr{$c}\n";
}
print "\n";
print "#CHR\tPOS\tREF\tALT\n";
my %linkeness;
foreach my $loc(sort keys %Snp)
{
	if (keys %{$Snp{$loc}}>=2)
	{
		print (join "\t",(split /\-/,$loc));
		my @SM=sort keys %{$Snp{$loc}};
		foreach my $sm(@SM)
		{
			print "\t$sm:$Snp{$loc}{$sm}";
		}
		print "\n";
		foreach my $A(@SM)
		{
			foreach my $B(@SM)
			{
				next if ($A eq $B);
				if ($Snp{$loc}{$A} eq $Snp{$loc}{$B})
				{
					$linkeness{$A}{$B}++;
					$linkeness{$B}{$A}++;
				}
			}
		}
	}
}
print "\n";
print "#SampleA\tSampleB\tSameSNPLocus\tNormalizationRate( SameSNPLoc/(SampleA+SampleB-SameSNPLoc) )\n";
my %count=();
foreach my $A(sort keys %linkeness)
{
	foreach my $B(sort keys %{$linkeness{$A}})
	{
		next if (exists $count{$A}{$B} || exists $count{$B}{$A});
		$count{$A}{$B}=1;
		$count{$B}{$A}=1;
		print ("$A\t$B\t$linkeness{$A}{$B}\t",($linkeness{$A}{$B}/($hash{$A}{'snp'}+$hash{$B}{'snp'}-$linkeness{$A}{$B})),"\n");
	}
}
print "\n";
foreach my $loc(sort keys %Indel)
{
	if (keys %{$Indel{$loc}}>=2)
	{
		print (join "\t",(split /\-/,$loc));
		foreach my $sm(@SM)
		{
			print "\t$sm:$Indel{$loc}{$sm}";
		}
		print "\n";
	}
}