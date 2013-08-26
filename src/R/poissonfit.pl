#Author: BENM <binxiaofeng@gmail.com>
#v1.0, 2013-06-17
#!/usr/bin/perl 
use strict;
use warnings;
use Data::Dumper;

die "usage:\npoissonfit.pl <file.tab> <output_prefix>\nRequired R\nAuthor:BENM <binxiaofeng\@gmail.com> v1.0, 2013-06-07\n" if (@ARGV != 2);

my ($tabfile,$prefix)=@ARGV;
my (%sig,%all,%sample,%sample_sd,%possion,%possion_sd);
#
my @sum_sig;
my @sum_all;
my @sum_sample_sd;
my @sum_possion_sd;

my $num=0;
open (IN,$tabfile)|| die "$!\n";
while(<IN>){
	chomp;
	my @t=split /\t+/,$_;
	for (my $i=1;$i<@t;$i++) {
		$sig{$t[0]}{$i-1}=$t[$i];
		$all{$t[0]}{$i-1}=$t[0]*$t[$i];
		$sum_sig[$i-1]+=$t[$i];
		$sum_all[$i-1]+=$t[0]*$t[$i];
	}
	$num++;
}
close IN;

my @avg=();
my $col_num=scalar(@sum_sig);
for (my $i=0;$i<$col_num;$i++) {
	$avg[$i]=$sum_all[$i]/$sum_sig[$i];
}
#my $aver_cov=$sum_all/$sum_sig;
my %gt=("1"=>0,"2"=>0,"4"=>0,"5"=>0,"10"=>0,"20"=>0,"30"=>0,"40"=>0,"50"=>0,);
my $out_data1="$prefix.plot.out";
my $out_data2="$prefix.poisson_plot.txt";
open O1,">$out_data1" || die "Can't open the file:$!\n";
open O2,">$out_data2" || die "Can't open the file:$!\n";
foreach my $n (sort{$a<=>$b} keys %sig){
	my $out=$n;
	for (my $i=0;$i<$col_num;$i++) {
		$sample{$n}=$sig{$n}{$i}/$sum_sig[$i];
		$sample_sd{$n}=($n-$avg[$i])*($n-$avg[$i])*$sample{$n};
		$possion{$n}=&poisson($n,$avg[$i]);
		$possion_sd{$n}=($n-$avg[$i])*($n-$avg[$i])*$possion{$n};
		if($possion{$n} >= 1E-10){
			$sum_sample_sd[$i]+=$sample_sd{$n};
			$sum_possion_sd[$i]+=$possion_sd{$n};
			$out.="\t$sample{$n}\t$possion{$n}";
		} else {
			$out.="\t0\t0";
		}
		foreach my $k (keys %gt){
			if($n >= $k){
				$gt{$k}+=$sig{$n};
			}
		}
		my $j=$i+1;
		print O1 "Col:$i\t$n\t$sig{$n}{$i}\t$all{$n}{$i}\t$sum_sig[$i]\t$sample{$n}\t$possion{$n}\t$sample_sd{$n}\t$possion_sd{$n}\n";
	}
	print O2 "$out\n";
}
close O1;
close O2;


my $out_data3="$prefix.plot.stat";
open O3,">$out_data3" || die "Can't open the file:$!\n ";
foreach my $n ( sort { $a <=> $b } keys %gt){
	my $out=" >= $n :";
	for (my $i=0;$i<@sum_sig;$i++) {
		$out.="\t".sprintf("%.4f",($gt{$n}/$sum_sig[$i]*100));
	}
	print O3 "$out\n";
}
my @raw_data=();
my @possion_data=();
for (my $i=0;$i<$col_num;$i++) {
	my $j=$i+1;
	print O3 "########## -- col$j -- ##########\n";
	print O3 "sum of single base : $sum_sig[$i]\nsum of all base : $sum_all[$i]\nmean of coverage : $avg[$i]\n";
	print O3 "sum of sample_sd >= 1E-10 : $sum_sample_sd[$i]\nsum of poisson_sd >= 1E-10 : $sum_possion_sd[$i]\n";
	print O3 "square root of sum of sample_sd >= 1E-10 : ",sqrt($sum_sample_sd[$i]);
	print O3 "\nsquare root of sum of poisson_sd >= 1E-10 : ",sqrt($sum_possion_sd[$i]);
	print O3 "\nratio : ",(sqrt($sum_sample_sd[$i]))/(sqrt($sum_possion_sd[$i]));
	print O3 "\n";
	
	my $m=2*$i+2;
	my $n=2*$i+3;
	push @raw_data,"data\[\,$m\]";
	push @possion_data,"data\[\,$n\]";
}
close O3;
my $DATA1=join ",",@raw_data;
my $DATA2=join ",",@possion_data;
my $tentimes=$num*10;
my $Rscript=qq(
library\("ggplot2"\)
data<-read.table\("./$prefix.poisson_plot.txt",head=F\)
bardata=data.frame\(x=rep\(1:$num,times=$col_num\),y=c\($DATA1\),c=rep\(1:$col_num,each=$num\)\)
);
my @sp_data;
for (my $i=0;$i<$col_num;$i++) {
	my $j=$i+1;
	my $n=2*$i+3;
	$Rscript.=qq(
sp$j<-spline\(data[,1],data[,$n],n=$tentimes\)
sp_data$j = data.frame\(sp$j\)
	);
	push @sp_data,"sp\_data$j\[\,2\]";
};
my $DATA3=join ",",@sp_data;

$Rscript.=qq(
poissonfit=data.frame\(x=rep\(1:$tentimes,times=$col_num\)/10,y=c\($DATA3\),c=rep\(1:$col_num,each=$tentimes\)\)
ggplot\(data=bardata,aes\(x,y,fill=factor\(c\)\)\)\+
  geom_bar\(stat="identity",position="dodge",\)\+
  scale_x_continuous\(breaks=c\(1:$num\)\)\+
  opts\(axis.text.x = theme_text\(angle = 90, hjust = 1,size = 10\)\)+
  geom_line\(data=poissonfit,aes\(color=factor\(c\)\),size=1\)\+
ggsave\(filename="./$prefix.poisson.png"\)
);

open O1,">$prefix.Rscript.R" || die "Can't open the file:$!\n";
print O1 "$Rscript";
close O1;
system("Rscript $prefix.Rscript.R");
#system("rm $prefix.Rscript.R $prefix.poisson_plot.txt");


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

