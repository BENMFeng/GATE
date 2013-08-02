#!/usr/bin/perl -w
#Author: BENM <binxiaofeng@gmail.com>
#v0.1.2, 2013-03-21
use strict;
use Getopt::Long;
my %opts;
my ($Uniq,$Num,$RG,$Help);
GetOptions(\%opts,"uniq"=>\$Uniq,"n"=>\$Num,"rg"=>\$RG,"help"=>\$Help);
die "perl $0 <IN:SAM|pipe_SAM_out> [-u required uniquely mapped] [-n use number instead reads name] [-rg add RG in readsId] [-help]\n" if ($Help);
my $k=0;
while(<>){
	my @t=split /\t+/,$_;
	next if (($t[1] & 0x4) == 4 || $t[5] eq "*");
	my $nbBestHit = 0; # intialize like unique hits if BWA wasn't used for mapping.
	my $xt = 'U';
	if($_ =~ m/XT:A:(\w+)\s/) { $xt = $1;} 
	if($_ =~ m/X1:i:(\d+)\s/) { $nbBestHit = $1;} 
	next if ((defined $Uniq) && ($xt ne "U" || $nbBestHit>0 || $_=~/XA\:Z\:/));
	my $rg=$1 if (/RG\:Z\:(\S+)/);
	my $strand=(($t[1] & 0x10)>0)?"-":"+";
	if (defined $Num) {
		$k++;
		$t[0]=sprintf("%09d",$k);
	}
	else
	{
		$t[0].=(($t[1] & 0x80) == 128)?"\/2":"\/1";
	}
	my $loc=$t[3];
	my $pos=$t[3];
	my $i=0;
	if (defined $RG)
	{
		$t[0]="$rg:$t[0]" if (defined $rg && $rg ne "");
	}
	if ($t[5] eq "=") {
		$pos+=length($t[9]);
		print ("$t[2]\t$loc\t",$pos-1,"\t$t[0]\t$i\t$strand\n");
	} else {
		while($t[5]=~/(\d+)(\D+)/g){
			my ($match,$cigar)=($1,$2);
			if (($cigar eq "M")||($cigar eq "S")){
				$pos+=$match;
			}elsif($cigar ne "I"){
				if ($pos>$loc){print ("$t[2]\t$loc\t",$pos-1,"\t$t[0]\t$i\t$strand\n");$i++;}
				$pos+=$match;
				$loc=$pos+1;
			}
			else{
				if ($pos>$loc){print ("$t[2]\t$loc\t",$pos-1,"\t$t[0]\t$i\t$strand\n");$i++;}
				if ($pos>$loc){$loc=$pos}
			}
		}
		if ($pos>$loc){print ("$t[2]\t$loc\t",$pos-1,"\t$t[0]\t$i\t$strand\n");}
	}
}
