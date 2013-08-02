#Author: BENM <binxiaofeng@gmail.com>
#v1.0, 2013-06-03
#!/usr/bin/perl 
use strict;
use warnings;

die "usage:\nplotCoverageDepth.pl <coverage_depth_cov.tab> <output_prefix>\nnRequired R\nAuthor:BENM <binxiaofeng\@gmail.com> v1.0, 2013-06-03\n" if (@ARGV != 2);

my ($tabfile,$prefix)=@ARGV;
my (%cov_sigbase,%cov_allbase,%cov_sample,%cov_sample_sd,%cov_possion,%cov_possion_sd);
my ($sum_sig,$sum_all,$sum_sample_sd,$sum_possion_sd,$num);

open (IN,$tabfile)|| die "$!\n";
while(<IN>){
	chomp;
	my @d=split;
	$cov_sigbase{$d[0]}=$d[1];
	$cov_allbase{$d[0]}=$d[0]*$d[1];
	$sum_sig+=$d[1];
	$sum_all+=$d[0]*$d[1];
	$num++;
}
close IN;

my $aver_cov=$sum_all/$sum_sig;
my %cov_gt=("1"=>0,"2"=>0,"4"=>0,"5"=>0,"10"=>0,"20"=>0,"30"=>0,"40"=>0,"50"=>0,);
my $out_data1="$prefix.coverage_plot.out";
my $out_data2="$prefix.poisson_plot.txt";
open O,">$out_data1" || die "Can't open the file:$!\n";
open O2,">$out_data2" || die "Can't open the file:$!\n";
foreach (sort{$a<=>$b} keys %cov_sigbase){
	$cov_sample{$_}=$cov_sigbase{$_}/$sum_sig;
	$cov_sample_sd{$_}=($_-$aver_cov)*($_-$aver_cov)*$cov_sample{$_};
	$cov_possion{$_}=&poisson($_,$aver_cov);
	$cov_possion_sd{$_}=($_-$aver_cov)*($_-$aver_cov)*$cov_possion{$_};
	if($cov_possion{$_} >= 1E-10){
		$sum_sample_sd+=$cov_sample_sd{$_};
		$sum_possion_sd+=$cov_possion_sd{$_};
		print O2 "$_\t$cov_sample{$_}\t$cov_possion{$_}\n";
	}
	foreach my $cov (keys %cov_gt){
		if($_ >= $cov){
			$cov_gt{$cov}+=$cov_sigbase{$_};
		}
	}
	print O "$_\t$cov_sigbase{$_}\t$cov_allbase{$_}\t$sum_sig\t$cov_sample{$_}\t$cov_possion{$_}\t$cov_sample_sd{$_}\t$cov_possion_sd{$_}\n";
}
close O2;
close O;

my $out_data3="$prefix.coverage_plot.stat";
open O3,">$out_data3" || die "Can't open the file:$!\n ";
my @a=keys %cov_gt;
my @b=sort { $a <=> $b } @a;
foreach my $cov ( @b){
	print O3 "Depth >= $cov : ",$cov_gt{$cov}/$sum_sig*100;
	print O3 "\n";
}
print O3 "###################################\n";
print O3 "sum of single base : $sum_sig\nsum of all base : $sum_all\nmean of coverage : $aver_cov\n";
print O3 "sum of sample_sd >= 1E-10 : $sum_sample_sd\nsum of poisson_sd >= 1E-10 : $sum_possion_sd\n";
print O3 "square root of sum of sample_sd >= 1E-10 : ",sqrt($sum_sample_sd);
print O3 "\nsquare root of sum of poisson_sd >= 1E-10 : ",sqrt($sum_possion_sd);
print O3 "\nratio : ",(sqrt($sum_sample_sd))/(sqrt($sum_possion_sd));
print O3 "\n";
close O3;

sub poisson{
	my $x=shift;
	my $m=shift;
	my $Rscript=<<RS;
	ps<-dpois($x, $m, log = FALSE)
	as.character(ps)->ps
	write(ps,"$prefix.poisson_result.txt")
RS
	open OUT,">$prefix.Rscript.R" || die "Can't open the file:$!\n";
	print OUT "$Rscript";
	close OUT;
	system("Rscript $prefix.Rscript.R");
	system("rm $prefix.Rscript.R");
	open IN,"$prefix.poisson_result.txt" || die "Can't open the file:$!\n";
	my $p=<IN>;
	chomp $p;
	close IN;
	system("rm $prefix.poisson_result.txt");
	return $p;
}

my $Rscript=<<RS1;
	f<-read.table("./$prefix.poisson_plot.txt",head=F)
	pdf("./$prefix.coverage_plot.pdf",width=16,height=12);
	par(mar=c(5,5,5,5))
	plot(f[,1],f[,3],cex.lab=1.5,cex.axis=1.2,font.lab=2,font.axis=2,main="Depth Distribution",xlab="Depth",ylab="Percent",pch=".",cex.main=2.5,col="red")
	sp<-spline(f[,1],f[,3],n=1000)
	lines(sp,col="red",lwd=2)
	sp1<-spline(f[,1],f[,2],n=1000)
	lines(sp1,col="blue",lwd=2)
	legend("topright",pch=19,col=c("red","blue"),c("poisson","$prefix"))
RS1

open O1,">$prefix.Rscript.R" || die "Can't open the file:$!\n";
print O1 "$Rscript";
close O1;
system("Rscript $prefix.Rscript.R");
#system("rm $prefix.Rscript.R $prefix.poisson_plot.txt");

