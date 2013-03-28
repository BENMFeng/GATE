#!/usr/bin/perl -w
use strict;

die "perl $0 <IN1:gene.gff> <IN2:aln.bam> <IN:sam2bed.pl>\n" if (@ARGV!=3);

my ($Gff,$Aln,$sam2bed)=@ARGV;
my %Gene;
open (IN,$Gff);
my $geneName="";
while(<IN>)
{
	chomp;
	my @t=split /\t+/,$_;
	next if (/^\#/ || @t<8);
	#if ($t[2] =~ /gene/i || $t[2] =~ /mRNA/i)
	if ($t[2] =~ /gene/i)
	{
		if ($t[8]=~/Name\=([^\=\;\s]+)\.g/ || $t[8]=~/Name\=([^\=\;\s]+)/)
		{
			$geneName=$1;
			die "$_\n" if (exists $Gene{$t[0]}{$1}{'g'} && (${$Gene{$t[0]}{$1}{'g'}}[0]!=$t[3] || ${$Gene{$t[0]}{$1}{'g'}}[1]!=$t[4]) );
			push @{$Gene{$t[0]}{$geneName}{'g'}},($t[3],$t[4]);
		}
	}
	elsif ($t[2] =~ /CDS/i || $t[2] =~ /exon/i || $t[2]=~/UTR/i)
	{
		#if ($t[8]=~/Parent\=([^\=\;\s]+)/)
		#{
		next if (exists $Gene{$t[0]}{$geneName}{'c'} && (${$Gene{$t[0]}{$geneName}{'c'}}[-1][0]==$t[3] || ${$Gene{$t[0]}{$geneName}{'c'}}[-1][1]==$t[4]) );
		push @{$Gene{$t[0]}{$geneName}{'c'}},[$t[3],$t[4]];
		#}
	}
}
close IN;

foreach my $id(keys %Gene)
{
	foreach my $g(keys %{$Gene{$id}})
	{
		my ($alnlen,$cdslen)=(0,0);
		die "$g\n" if (!exists $Gene{$id}{$g}{'c'} || @{$Gene{$id}{$g}{'c'}}==0);
		for (my $i=0;$i<@{$Gene{$id}{$g}{'c'}};$i++)
		{
			my ($s,$e)=@{${$Gene{$id}{$g}{'c'}}[$i]};
			my $cmd=qq(samtools view -F 4 -q 10 $Aln $id:$s-$e | $sam2bed -u -n |awk '{a+=\$3-\$2+1}END{print a}');
			my $count=`$cmd`;
			chomp $count;
			if (defined $count && $count ne '' && $count>0)
			{
				$alnlen+=$count;
			}
			$cdslen+=$e+$s+1;
		}
		print ("$g\t",$alnlen/$cdslen,"\n");
	}
}