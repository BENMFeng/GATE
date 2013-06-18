#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my ($Edge,$Detail,$Len_p,$Len_t,$SubStr,$TrimMode,$Homopolymer,$APmode,$APdb,$Mismatch,$Prefix,$Verbose,$Help);
my %opts;
GetOptions(\%opts,"trim_detail:s"=>\$Detail,"edge:i"=>\$Edge,"len_p:s"=>\$Len_p,"len_t:i"=>\$Len_t,"mis:i"=>\$Mismatch,
"sub:i"=>\$SubStr,"trim_mode:s"=>\$TrimMode,"trim_polymer:s"=>\$Homopolymer,"ap_mode"=>\$APmode,"ap_db:s"=>\$APdb,"prefix:s"=>\$Prefix,
"verbose"=>\$Verbose,"help"=>\$Help);
my $usage = qq(
trim_seq.pl -- trim fasta or fastq according to filter_adapter details
Usage: perl $0 [option] <Input files>
Option: --trim_detail <file>		filter_adapter detail file, include: "reads_id   reads_len   reads_start   reads_end   adapter_id"
        --edge <int>			adapter or primer align positon distance from the edge of reads, default: 2
        --len_p <float>			set the minimal bases percentage rate of whole reads after trimmed, default: 0.5
        --len_t <int>			set the minimal bases pair of the reads, default: 35
        --sub <int-int>			set reads substr: left_length-right_length
        --trim_mode <3|5|both>		trim 3\' or 5\' or both ends of reads, default: both
        --trim_polymer <Type:Len>	trim homerpolymers, default: NULL (example: AT:6)
        --ap_mode			trim all adapters & primers by merged mapping sequences
        --ap_db				input AP file for strictly trimming
        --mis <int>			allow mismatch and gaps bases, default: 2
        --prefix <str>			set output prefix
        --verbose			report runing process
        --help				help information
Example: perl trim_fastq.pl *.fastq -trim_detail PE.adapter-primer.detail --edge 3 --trim_mode 3 --ap_mode
Author : BENM <BinxiaoFeng\@gmail.com>
Version: 0.2.1 beta
Data: v1.0 Nov-10-2010
Update: v1.1 Dec-27-2010
Update: v1.2 Apr-28-2012
Update: v1.3 May-02-2012
Update: v1.4 May-03-2012
Update: v1.5 Oct-29-2012
Update: v1.6 Nov-02-2012
Update: v1.7 Nov-19-2012
Update: v1.8 Nov-21-2012
Update: v1.9 Dec-03-2012
Update: v2.0 Dec-08-2012
Update: v2.1 Apr-16-2013
);

die $usage if ((@ARGV==0)||($Help));
#reads_id reads_len reads_start reads_end AP_id AP_len AP_start AP_end align_len mismatch gap strand
#HWI-EAS59:1:1:384:318#0/1	40	6	17	PE_Adapter1-5'+	33	2	13	11	0	0	+
#HWI-EAS59:1:2:1645:544#0/1	40	30	40	PE_Adapter1-5'+	33	3	13	10	0	0	+
#HWI-EAS59:1:1:992:104#0/1	40	13	24	PE_Adapter1-5'+	33	0	11	11	0	0	+

$Edge ||= 2;
$Len_p ||= 0.5;
$Len_t ||= 35;
$Mismatch ||= 2;
$TrimMode ||= "both";
my ($Left,$Right) = (0,0);
($Left,$Right) = split /\-/,$SubStr if (defined $SubStr);

my ($PolymerBase,$PolymerLen)=split /\:/,$Homopolymer if (defined $Homopolymer);
die "--trim_polymer setting error: $Homopolymer, example: AT:10\n" if (defined $Homopolymer && ($PolymerBase!~/[ACGT]/ || $PolymerLen !~/\d+/));

my ($totalAP,$totalAln,$totalMis,$totalGap)=(0,0,0,0);

my %APSeq=();
my %APname=();
my @APtitle=();
my $apcount=0;
if (defined $APdb && -f $APdb){
	my ($name,$seq)=("","");
	open (AP,$APdb) || die $!;
	while(<AP>){
		chomp;
		next if ($_ eq "");
		if (/^\>(.*)/ || eof){
			$seq.=$_ if (eof);
			if ($seq ne "" && length($seq)>=5) {
				my $seed=substr($seq,0,5);
				push @{$APSeq{$seed}},[$name,$seq];
			}
			$name=$1;
			$APname{$name}=$apcount;
			push @APtitle,$apcount;
			$apcount++;
			$seq="";
		} else {
			$seq.=$_;
		}
	}
	close AP;
}

my %Trim=();
my $ExpetTrimReads=0;
if (defined $Detail)
{
	print STDERR (get_time()," Reading adapter align details file: $Detail...\n")  if (defined $Verbose);
	open (D,$Detail) || die "Can't open $Detail file for reading\n";
	while(<D>)
	{
		chomp;
		my @t=split /\t+/,$_;
		next if (/^\#/ || @t<10);
		chomp($t[0]);
		$totalAln+=$t[8];
		$totalMis+=$t[9];
		$totalGap+=$t[10];
		next if (($t[0]=~/reads\_id/)||($t[9]+$t[10]>$Mismatch));
		my $seqname=$1 if ($t[0]=~/\:(\d+\:\d+\:\d+\:\d+)/ || $t[0]=~/^(\S+)/);
		$t[4]=$APname{$t[4]} if (exists $APname{$t[4]});
		if (!defined $APmode)
		{
			if ((($t[2]<=$Edge)||($t[1]-$t[3]<=$Edge))&&(($t[6]<=$Edge)||($t[5]-$t[7]<=$Edge)))
			{
				if ((!exists $Trim{$seqname})||(${$Trim{$seqname}}[2]-${$Trim{$seqname}}[1]+1<$t[8]))
				{
					@{$Trim{$seqname}}=@t[1..4];
					$totalAP++;
				}
			}
		}
		else
		{
			if ((($t[2]<=$Edge)||($t[1]-$t[3]<=$Edge)||($t[8]==$t[5]))||(($t[6]<=$Edge)||($t[5]-$t[7]<=$Edge)))
			{
				my $indexLen=0;
				$indexLen=length($1) if ($t[0]=~/\:([ACGT]+)$/);
				my $overlapLen=overlap([@{$Trim{$seqname}}[1..2]],[@t[2..3]]) if (exists $Trim{$seqname});
				if ( (exists $Trim{$seqname}) && (($overlapLen>=-$indexLen)||($APmode)) )
				{
					${$Trim{$seqname}}[1]=(${$Trim{$seqname}}[1]<$t[2])?${$Trim{$seqname}}[1]:$t[2];
					${$Trim{$seqname}}[2]=(${$Trim{$seqname}}[2]>$t[3])?${$Trim{$seqname}}[2]:$t[3];
					${$Trim{$seqname}}[3]=($t[4] !~ /${$Trim{$seqname}}[3]/) ? ${$Trim{$seqname}}[3] : ${$Trim{$seqname}}[3]."-".$t[4];
					$totalAP++;
				}
				elsif ( (!exists $Trim{$seqname}) || (${$Trim{$seqname}}[2]-${$Trim{$seqname}}[1]+1<$t[8]) )
				{
					@{$Trim{$seqname}}=@t[1..4];
					$totalAP++;
				}
			}
		}
	}
	close D;
	$ExpetTrimReads=scalar(keys %Trim);
}
else
{
	print STDERR (get_time()," No details file input!\n")  if (defined $Verbose);
}

my $i=0;
foreach my $file(@ARGV)
{
	print STDERR (get_time()," Dealing with $file...\n")  if (defined $Verbose);
	open (IN,$file) || die "Can't open $file fastq file for reading\n" if ($file!~/\.gz$/);
	open (IN,"zcat $file|") || die "Can't open $file fastq file for reading\n" if ($file=~/\.gz$/);
	my $outprefix=(defined $Prefix)? $Prefix : (split /\//,$file)[-1];
	open (T,">$outprefix.trim.out") || die "Can't write to $file.trim.out\n";
	open (F,">$outprefix.filter.out") || die "Can't write to $file.filter.out\n";
	my ($trim_left,$trim_right,$filter,$unpolluted)=(0,0,0,0);
	my $deal=0;
	while(<IN>)
	{
		if (/^([\@\>])(\S+)/)
		{
			$i++;
			print STDERR (get_time()," Sequences: $i\n")  if ((defined $Verbose)&&($i%1000000==0));
			my $type = ($1 eq "@") ? 0 : 1;
			my $name=$1 if ($_=~/\:(\d+\:\d+\:\d+\:\d+)/ || $_=~/^([\@\>])(\S+)/);
			s/\s+$//;
			my $Name=$_;
			if (exists $Trim{$name})
			{
				$deal=1;
			}
			else
			{
				$deal=0;
			}
			$_=<IN>;
			s/\s+$//;
			my $seq=$_;
			my $qual;
			if ($type==0)
			{
				<IN>;
				$_=<IN>;
				s/\s+$//;
				$qual=$_;
			}
			if (defined $qual && length($qual) != length($seq))
			{
				if (length($qual)>length($seq))
				{
					substr($qual,0,length($qual)-length($seq),"");
				}
				else
				{
					$qual.="#"x(length($seq)-length($qual));
				}
			}
			if ((defined $SubStr) && (($Left>0) || ($Right>0)) )
			{
				my ($left_seq,$left_qual,$right_seq,$right_qual)=("","","","");
				if (defined $Left && $Left>0)
				{
					$left_seq=substr($seq,0,$Left);
					$left_qual=substr($qual,0,$Left);
				}
				if (defined $Right && $Right>0)
				{
					$right_seq=substr($seq,length($seq)-$Right,$Right);
					$right_qual=substr($qual,length($seq)-$Right,$Right);
				}
				$seq=$left_seq.$right_seq;
				$qual=$left_qual.$right_qual;
			}
			if ($deal==1)
			{
				my $db_title=(defined $APdb && ${$Trim{$name}}[3]=~/^\d+$/ && defined $APtitle[${$Trim{$name}}[3]])?$APtitle[${$Trim{$name}}[3]]:${$Trim{$name}}[3];
				if ( (${$Trim{$name}}[0]-${$Trim{$name}}[2]<=$Edge) || (defined $APmode && ${$Trim{$name}}[0]-${$Trim{$name}}[2]<${$Trim{$name}}[1]-1) )
				{
					if (($TrimMode =~ /both/i)||($TrimMode eq "3"))
					{
						$Name .= (" $db_title 3\'_trimmed: ".(length($seq)-${$Trim{$name}}[1]));
						if (length($seq)-${$Trim{$name}}[1]>0)
						{
							substr($seq,${$Trim{$name}}[1]-$Left,length($seq)-${$Trim{$name}}[1],"");
							if ($type==0)
							{
								substr($qual,${$Trim{$name}}[1]-$Left,length($qual)-${$Trim{$name}}[1],"");
							}
						}
						else
						{
							$seq="";
							$qual="" if ($type==0);
						}
						$trim_right++;
					}
					else
					{
						$Name .= " $db_title filtered";
						if ($type==0)
						{
							print F "$Name\n$seq\n+\n$qual\n";
						}
						else
						{
							print F "$Name\n$seq\n";
						}
						$filter++;
						next;
					}
				}
				elsif (${$Trim{$name}}[1]<=$Edge || (defined $APmode && ${$Trim{$name}}[0]-${$Trim{$name}}[2]>${$Trim{$name}}[1]-1))
				{
					if (($TrimMode =~ /both/i)||($TrimMode eq "5"))
					{
						my $index=$1 if ($Name=~/\:([ACGT]+)$/);
						my $seq_len=length($seq);
						if ($APmode)
						{
							$seq="";
						}
						elsif (${$Trim{$name}}[2]<length($seq)+$Left)
						{
							substr($seq,0,${$Trim{$name}}[2]-$Left,"");
						}
						if (length($seq)>0)
						{
							if (defined $index && $index ne "" && length($index)>3)
							{
								my @Index=split "",$index;
								for (my $i=0;$i<@Index-1;$i++)
								{
									my @tmp=@Index;
									$tmp[$i]="." if ($i>0);
									my $pattern=join "",@tmp;
									if ($seq=~/^$pattern/)
									{
										$seq="";
										${$Trim{$name}}[2]=$seq_len;
										$Name .= "match Index after trimmed: $pattern";
										last;
									}
								}
							}
							if ($type==0)
							{
								if (${$Trim{$name}}[2]<length($qual)+$Left)
								{
									substr($qual,0,${$Trim{$name}}[2]-$Left,"");
								}
								else{
									$qual="";
								}
							}
						}
						else
						{
							$qual="" if ($type==0);
						}
						$Name .= (" $db_title 5\'_trimmed: ".(${$Trim{$name}}[2]-${$Trim{$name}}[1]));
						$trim_left++;
					}
					else
					{
						$Name .= " $db_title filtered";
						if ($type==0)
						{
							print F "$Name\n$seq\n+\n$qual\n";
						}
						else
						{
							print F "$Name\n$seq\n";
						}
						$filter++;
						next;
					}
				}
				trim_polymer(\$Name,\$seq,\$qual) if (defined $Homopolymer);
				if ((length($seq)<${$Trim{$name}}[0]*$Len_p)||(length($seq)<$Len_t))
				{
					$Name .= " $db_title filtered";
					if ($type==0)
					{
						print F "$Name\n$seq\n+\n$qual\n";
					}
					else
					{
						print F "$Name\n$seq\n";
					}
					$filter++;
				}
				else
				{
					if (defined $APdb){
						strict_trim(\$Name,\$seq,\$qual);
					}
					if (length($seq)>=$Len_t)
					{
						if ($type==0)
						{
							print T "$Name\n$seq\n+\n$qual\n";
						}
						else
						{
							print T "$Name\n$seq\n";
						}
					}
					else
					{
						if ($type==0)
						{
							print F "$Name\n$seq\n+\n$qual\n";
						}
						else
						{
							print F "$Name\n$seq\n";
						}
					}
				}
				@{$Trim{$name}}=();
				delete $Trim{$name};
			}
			else
			{
				if (defined $APdb){
					strict_trim(\$Name,\$seq,\$qual);
				}
				if (defined $Homopolymer)
				{
					trim_polymer(\$Name,\$seq,\$qual);
				}
				if (length($seq)<$Len_t)
				{
					$Name .= " filtered";
					if ($type==0)
					{
						print F "$Name\n$seq\n+\n$qual\n";
					}
					else
					{
						print F "$Name\n$seq\n";
					}
					$filter++;
				}
				else
				{
					if ($type==0)
					{
						print T "$Name\n$seq\n+\n$qual\n";
					}
					else
					{
						print T "$Name\n$seq\n";
					}
					$unpolluted++;
				}
			}
		}
	}
	close IN;
	close T;
	close F;
	print STDERR (get_time()," Output report...\n")  if (defined $Verbose);
	if (defined $Detail)
	{
		open (OD,">$outprefix.filter-trim.stat") || die "Can't write to $file.filter-trim.stat.out\n";
		print OD "Total matched AP length: $totalAln\n";
		print OD "Total mismatch length between reads and AP: $totalMis\n";
		print OD "Total gap length between reads and AP: $totalGap\n";
		print OD ("Estimate error rate: ",int(1000000*($totalMis+$totalGap)/$totalAln+0.4999)/10000,"%\n");
		print OD ("Contain AP sequences number: $totalAP\n");
		print OD ("Expected trim reads number: $ExpetTrimReads\n");
		print OD "5\' trimmed: $trim_left\n";
		print OD "3\' trimmed: $trim_right\n";
		print OD ("All trimmed: ",($trim_left+$trim_right),"\n");
		print OD "Filtered: $filter\n";
		print OD "No polluted: $unpolluted\n";
		print STDERR (get_time()," Done!\n")  if (defined $Verbose);
	}
}

sub overlap
{
	my ($block1_p,$block2_p) = @_;
	my ($combine_start, $combine_end) = (sort {$a<=>$b} ($block1_p->[0],$block2_p->[0],$block1_p->[1],$block2_p->[1]))[0,-1];
	my $overlap_size = ($block1_p->[1]-$block1_p->[0]+1) + ($block2_p->[1]-$block2_p->[0]+1) - ($combine_end-$combine_start+1);
	return $overlap_size;
}

sub get_time
{
	my  ($sec,$min,$hour,$mday,$mon,$year) = (localtime)[0..5];
	($sec,$min,$hour,$mday,$mon,$year) = (
		sprintf("%02d", $sec),
		sprintf("%02d", $min),
		sprintf("%02d", $hour),
		sprintf("%02d", $mday),
		sprintf("%02d", $mon + 1),
		$year + 1900
	);
	return "[$year-$mon-$mday $hour:$min:$sec]";
}

sub strict_trim{
	my ($name,$seq,$qual)=@_;
	foreach my $seed(keys %APSeq) {
		while($$seq=~/($seed)/g){
			my $leftpart=$`;
			my $rightpart=qq($1$');
			my $len=length($rightpart);
			for (my $i=0;$i<@{$APSeq{$seed}};$i++){
				my ($spname,$apseq)=@{${$APSeq{$seed}}[$i]};
				if ($len<length($apseq)){
					if ($rightpart eq substr($apseq,0,$len)){
						$$name .= " $spname 3\'_trimmed: $len";
						$$seq=$leftpart;
					}
				} else {
					if ($rightpart =~ /($apseq)/) {
						my $trimlen=length(qq($1$'));
						$$name .= " $spname 3\'_trimmed: $trimlen";
						$$seq=$`;
					}
				}
			}
		}
	}
	if (defined $$qual && $$qual ne "")
	{
		if (length($$seq)<length($$qual))
		{
			substr($$qual,length($$seq),length($$qual)-length($$seq),"");
		}
	}
}

sub trim_polymer
{
	my ($name,$seq,$qual)=@_;
	my $polymerlen = $PolymerLen-1;
	while ($$seq=~/([^N\-]{1})(\1{$polymerlen,})/g){
		my ($element,$repeat,$rail,$pos)=($1,$2,$',$-[1]);
		my $primary_element=$element;
		while ($element=~/([ACGT]{1,}?)(\1{1,})/g)
		{
			$primary_element=$1;
			if ($rail=~/^(($primary_element)+)/)
			{
				$repeat.=$1;
			}
		}
		my $len=length("$element$repeat");
		if ($PolymerBase=~/$primary_element/i && $len>=$PolymerLen)
		{
			if ($pos+$len>=length($$seq)||($pos==0&&$TrimMode ne "3"))
			{
				if ($pos==0)
				{
					$$name.=" $primary_element:$len:$pos-".($pos+$len-1);
					substr($$seq,$pos,$len,"N"x($len));
				}
				else
				{
					$$name.=" $primary_element:$len:$pos-".($pos+$len-1);
					substr($$seq,$pos,$len,"");
					substr($$qual,$pos,$len,"") if (defined $$qual && $$qual ne "");
				}
			}
			else
			{
				if ($TrimMode ne "3")
				{
					$$name.=" $primary_element:$len:$pos-".($pos+$len-1);
					substr($$seq,$pos,$len,"N"x($len));
				}
			}
		}
	}
	if ($TrimMode ne "3" && $$seq=~/^(N+)/)
	{
		$$seq=$';
		if (defined $$qual && $$qual ne "")
		{
			substr($$qual,0,length($$qual)-length($$seq),"");
		}
	}
}
