#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;


my %opts;
my ($kmer,$prefix,$help);
GetOptions(
	\%opts,
	"k:i"=>\$kmer,
	"prefix:s"=>\$prefix,
	"help!"=>\$help,
);
die "perl $0 <IN1:target.fa> <IN2:gaps.fa> [-prefix output_prefix] [-k 5]\n" if (@ARGV<2 || $help);


my ($target_file,$gap_file)=@ARGV;
$kmer ||= 5;
my %hash=();
my %GapSeq=();
my %exist=();
my $Max=$kmer;
fetchKmer($gap_file,$kmer);

my ($Name,$Seq)=("","");
open(FA, $target_file) or die $!;
open(OUT, ">$prefix.fa") if (defined $prefix);
while (<FA>) {
	s/\s*$//;
	if (/>(.*)/) {
		if ($Seq ne "" && length($Seq)>0) {
			my $fillgapseq=matchSeq($Name,$Seq);
			if (defined $prefix) {
				print OUT ">$Name\n";
				print OUT Display($fillgapseq);
			}
		}
		$Name=$1;
		$Seq="";
	} else {
		$Seq.=$_;
	}
}
close FA;

if ($Seq ne "" && length($Seq)>0) {
	my $fillgapseq=matchSeq($Name,$Seq);
	if (defined $prefix) {
		print OUT ">$Name\n";
		print OUT Display($fillgapseq);
	}
}

close OUT if (defined $prefix);


sub fetchKmer {
	my ($file,$k)=@_;
	my ($name,$seq)=("","");
	open(GAP,$file) or die $!;
	while (<GAP>) {
		s/\s*$//;
		if (/>(.*)/) {
			if ($seq ne "" && length($seq)>0) {
				subkmer($seq,$name,$k,0);
				my $seq_rc=rc($seq);
				subkmer($seq_rc,$name,$k,2);
			}
			$name=$1;
			$seq="";
		} else {
			$GapSeq{$name}.=uc($_);
			$seq.=uc($_);
		}
	}
	close GAP;
	if ($seq ne "" && length($seq)>0) {
		subkmer($seq,$name,$k,0);
		my $seq_rc=rc($seq);
		subkmer($seq_rc,$name,$k,2);
	}
}

sub rc {
	my $seq=shift;
	my $out;
	$out=reverse($seq);
	$out=~tr/ACGTacgt/TGCAtgca/;
	return $out;
}

sub subkmer {
	my ($seq,$name,$k,$n)=@_;
	my $max=($k*2<length($seq))?$k*2:(length($seq)/2>$k)?length($seq)/2:$k+1;
	for (my $i=$k;$i<=$max;$i++) {
		for (my $j=$i;$j<=$max;$j++){
			my $subseq=substr($seq,0,$j);
			my $subseq_rc=reverse($subseq);
			$subseq_rc=~tr/ACGTactg/TGCAtgca/;
			if (!exists $hash{$subseq}  && !exists $exist{$subseq} ) {
				#print STDERR "$subseq\n";
				$hash{$subseq}{$n+0}=$name;
				$hash{$subseq}{$n+1}="0-$j";
				$Max=($j+1>$Max)?$j+1:$Max;
				$exist{$subseq}=1;
				
			} else {
			#	#print $subseq,"\t",$n,"\t",$name,"\t",$i."-".($i+$kmer-1),"\n";
				delete $hash{$subseq};
				$exist{$subseq}=1;
			}
			if (!exists $hash{$subseq_rc} && !exists $exist{$subseq_rc}) {
				$hash{$subseq_rc}{$n+0}=$name;
				$hash{$subseq_rc}{$n+1}="$j-0";
				$exist{$subseq_rc}=1;
				$Max=($j+1>$Max)?$j+1:$Max;
			} else {
				delete $hash{$subseq_rc};
				$exist{$subseq}=1;
			}
		}
	}
}

sub matchSeq {
	my ($name,$seq)=@_;
	my $new_seq=$seq;
	my @gaps=();
	while ($seq=~/([Nn]+)/g) {
		push @gaps,[$-[0],$-[0]+length($1)-1];
	}
	my @cns=();
	if ($gaps[0][0]>0) {
		push @cns,[0,$gaps[0][0]-1];
	}
	my $pre_end=$gaps[0][1]+1;
	for (my $i=1;$i<$#gaps;$i++) {
		push @cns,[$pre_end,$gaps[$i][0]-1];
		$pre_end=$gaps[$i][1]+1;
	}
	if ($gaps[-1][1]<length($seq)-1) {
		push @cns,[$pre_end,length($seq)-1];
	}
	my %pass;
	for (my $i=0;$i<$#cns;$i++) {
		my $j=$i+1;
		#print "## $cns[$i][0]-$cns[$i][1]\t$cns[$j][0]-$cns[$j][1]\n";
		my $head=uc(substr($seq,$cns[$i][0],$cns[$i][1]-$cns[$i][0]+1));
		my $tail=uc(substr($seq,$cns[$j][0],$cns[$j][1]-$cns[$j][0]+1));
		my $max1=(length($head)>$Max)?$Max:length($head);
		my @head_seq=();
		my @tail_seq=();
		my @strand1=();
		my @strand2;
		my @pos1=();
		my @pos2=();
		for (my $m=$kmer;$m<$max1;$m++) {
			for (my $n=0;$n<=length($head)-$m;$n++) {
				my $subseq=substr($head,$n,$m);
				#print "$subseq\n";
				if (exists $hash{$subseq}) {
					#print "$subseq\t",length($head)-$m,"\t$m\n";
					#print Dumper $hash{$subseq};
					push @head_seq,$subseq;
					push @strand1,"+";
					push @pos1,[length($head)-$n,length($head)-($n+$m-1)];
				}
				my $subseq_rc=rc($subseq);
				if (exists $hash{$subseq_rc}) {
					#print "$subseq_rc\t",length($head)-$m,"\t$m\n";
					#print Dumper $hash{$subseq_rc};
					push @head_seq,$subseq_rc;
					push @strand1,"-";
					push @pos1,[length($head)-($n+$m-1),length($head)-$n];
				}
			}
		}
		my $max2=(length($tail)>$Max)?$Max:length($tail);
		for (my $m=$kmer;$m<$max2;$m++) {
			for (my $n=0;$n<=length($tail)-$m;$n++) {
				my $subseq=substr($tail,$n,$m);
				#print "$subseq\n";
				if (exists $hash{$subseq}) {
					#print "$subseq\t",length($head)-$m,"\t$m\n";
					#print Dumper $hash{$subseq};
					push @tail_seq,$subseq;
					push @strand2,"+";
					push @pos2,[$n,($n+$m-1)];
				}
				my $subseq_rc=rc($subseq);
				if (exists $hash{$subseq_rc}) {
					#print "$subseq_rc\t",length($head)-$m,"\t$m\n";
					#print Dumper $hash{$subseq_rc};
					push @tail_seq,$subseq_rc;
					push @strand2,"-";
					push @pos2,[$n+$m-1,$n];
				}
			}
		}
		for (my $m=0;$m<@head_seq;$m++) {
			my ($flag1,$flag2)=(0,0);
			for (my $n=0;$n<@tail_seq;$n++) {
				if (exists $hash{$head_seq[$m]}{2} && !exists $hash{$head_seq[$m]}{0}) {
					$flag1=2;
				} else {
					$flag1=0;
				}
				if (exists $hash{$tail_seq[$n]}{2} && !exists $hash{$tail_seq[$n]}{0}) {
					$flag2=2;
				} else {
					$flag2=0;
				}
				my $name1=$hash{$head_seq[$m]}{$flag1+0};
				my $loc1=$hash{$head_seq[$m]}{$flag1+1};
				my $name2=$hash{$tail_seq[$n]}{$flag2+0};
				my $loc2=$hash{$tail_seq[$n]}{$flag2+1};
				if ($name1 eq $name2 && $flag1 ne $flag2) {
					my $pm1=(sort{$b<=>$a}@{$pos1[$m]})[0];
					my $pm2=(sort{$b<=>$a}@{$pos2[$n]})[0];
					die "$name1\n" if (!exists $GapSeq{$name1});
					if ($pm1+$pm2+$cns[$j][0]-$cns[$i][1] == length($GapSeq{$name1})) {
						#print STDERR "$head_seq[$m]\t$strand1[$m]\t$pos1[$m][0]-$pos1[$m][1]\t$tail_seq[$n]\t$strand2[$n]\t$pos2[$n][0]-$pos2[$n][1]\n";
						#print STDERR "$name1\t$flag1\t$loc1\t$name2\t$flag2\t$loc2\n";
						my $fragment=substr($seq,$cns[$i][1]-$pm1+1,$pm1+$pm2+$cns[$j][0]-$cns[$i][1]);
						if (!exists $pass{$name1}) {
							print "## $cns[$i][0]-$cns[$i][1] - N(..) - $cns[$j][0]-$cns[$j][1]\n";
							print "$fragment\n";
							my ($x,$y)=(0,0);
							while ($fragment=~/([Nn]+)/g) {
								$x=$-[0];
								$y=length($1);
							}
							die "$fragment\n" if ($x==0 || $y==0);
							if ($flag1==0 && $flag2==2) {
								print (lc($GapSeq{$name1}),"\n\n");
								my $fillgaps=substr(lc($GapSeq{$name1}),$x,$y);
								substr($new_seq,$cns[$i][1]+1,$y,$fillgaps);
							} else {
								my $rc=reverse($GapSeq{$name1});
								$rc=~tr/ACGTtgca/TGCAtgca/;
								print (lc($rc),"\n\n");
								my $fillgaps=substr(lc($rc),$x,$y);
								substr($new_seq,$cns[$i][1]+1,$y,$fillgaps);
							}
							$pass{$name1}=1;
						}
					}
				}
			}
		}
	}
	return $new_seq;
}

sub Display 
{
	my $seq_p=shift;
	my $num ||= (@_) ? shift : 100;
	my $disp;
	$seq_p =~ s/\s$//;
	for (my $i=0;$i<length($seq_p);$i+=$num) {
		$disp .= substr($seq_p,$i,$num)."\n";
	}
	return $disp;
}