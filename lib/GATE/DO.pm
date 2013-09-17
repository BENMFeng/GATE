##########################################################################
#  Copyright (c) 2012 - 2013 - BENM(Binxiao) Feng                        #
#  All Rights Reserved                                                   #
#  Send all comments to BENM - BinxiaoFeng\@gmail.com                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  (at your option) any later version.                                   #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#  GNU General Public License for more details.                          #
#  You should have received a copy of the GNU General Public License     #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. #
##########################################################################

package GATE::DO;
use strict;

use vars qw($VERSION);
$VERSION = "1.1d,09-15-2013";
package GATE::Element;
#package GATE::Extension;
use GATE::ALIAS;
use GATE::Error;
use FindBin qw($Bin $Script);    ## Find me? I (Binxiao) am here.
use lib "$FindBin::Bin/../lib";  ## Bin live in the lib
use File::Basename qw(basename dirname);
use Cwd;

#########################################################
#                                                       #
#                    Configure                          #
#                                                       #
#########################################################

sub parseConfig($) {
	my $self = shift;
	my $file=GATE::Error::checkPath($self->{'-config'});
	my ($name,$index,$barcode,$lib,$mergepe)=("","","","",""); 
	my %RG;
	my %countSample;
	open (IN,$file) || die "Can't open config file:$file for reading\n";
	while(<IN>) {
		chomp;
		next if ($_ eq "" || $_=~/^\#/);
		if (/[\[\<](.*)[\]\>]/) {
			$name=GATE::ALIAS::alias($1);
			if ($name eq "LIB")
			{
				($index,$barcode,$lib,$mergepe)=("","","","");
				delete @RG{keys %RG};
				%RG=();
			}
			next;
		}
		next if ($name eq "");
		if ($name eq "INPUT") {
			if (/^([A-Z]{2})\=([^\=]+)/) {
				my ($rg,$info)=($1,$2);
				$RG{$rg}=$info;
				$lib=$1 if (/LB\=([^\=]+)/);
			} elsif (/([^\=]+)\=([^\=]+)/) {
				my ($lb,$path)=($1,$2);
				$self->{'LIB'}{$lib}{$name}{$lb}++;
				$self->{$lib}{$lb}[$self->{'LIB'}{$lib}{$name}{$lb}-1]=GATE::Error::checkPath($path);
			}
		} elsif ($name eq "LIB") {
			if (/Label\=([^\=]+)/ || /LIB\=([^\=]+)/) {
				$lib=$1;
			} elsif(/Index\=([^\=]+)/) {
				$index=$1;
			} elsif (/Barcode\=([^\=]+)/) {
				$barcode=$1;
			} elsif (/MergePE\=([^\=]+)/) {
				$mergepe=$1;
			} elsif (/^fq/i || /^fa/i) {
				if (!defined $lib) {
					die get_time()."Error with configure file: $file, undefined LB";
				}
				if (/([^\=]+)\=([^\=]+)/){
					my ($type,$path)=($1,$2);
					$countSample{"$name $lib $type"}++;
					${$self->{$name}{$lib}{$type}}[$countSample{"$name $lib $type"}-1]=GATE::Error::checkPath($path);
					$self->{$lib}{$type}{$countSample{"$name $lib $type"}-1}{'Index'}=$index if (defined $index && $index ne "");
					$self->{'idx'}{$lib}=1 if (defined $index && $index ne "");
					$self->{$lib}{$type}{$countSample{"$name $lib $type"}-1}{'Barcode'}=$barcode if (defined $barcode && $barcode ne "");
					$self->{'bar'}{$lib}=1 if (defined $barcode && $barcode ne "");
					$self->{$lib}{$type}{$countSample{"$name $lib $type"}-1}{'MergePE'}=$mergepe if (defined $mergepe && $mergepe =~ /^[TY]/i);
					$self->{'pemerge'}{$lib}=1 if (defined $mergepe && $mergepe ne "");
					$RG{'LB'}=$lib if (defined $lib && !exists $RG{'LB'});
					foreach my $rg(keys %RG)
					{
						$self->{$lib}{$type}{$countSample{"$name $lib $type"}-1}{$rg} = $RG{$rg};
					}
				}
			} elsif (/^([A-Z]{2})\=([^\=]+)/) {
				my ($rg,$info)=($1,$2);
				$RG{$rg}=$info;
				$lib=$1 if (/LB\=([^\=]+)/);
			}
		} elsif ($name =~ /rule/i) {
			if (/([^\=]+)\=([^\=]+)/) {
				my ($rule,$condition)=($1,$2);
				if ($rule =~ /skip/i) {
					my @cond=split /\s+/,$condition;
					map{$self->{"$name:$rule"}{$_}=1}@cond;
				} else {
					$self->{"$name:$rule"}=$condition;
				}
			}
		} else {
			if (/([^\=]+)\=([^\=]+)/) {
				my ($lib,$path)=($1,$2);
				if ($name =~ /database/ || $name =~ /software/) {
					$path=GATE::Error::checkPath($path); 
					$self->{"$name:$lib"} = $path;
				}
			}
		}
	}
	close IN;
	return bless($self);
}

sub parseDir($) {
	my $self=shift;
	my $Workdir = $self->{'-workdir'};
	`mkdir -p $Workdir` if (defined $Workdir && !-d $Workdir);
	$Workdir=~s/[\/\s]+$//;
	my $workpath=getcwd();
	$workpath=~s/[\/\s]+$//;
	$Workdir="$workpath/$Workdir" if ($Workdir !~ /^\//);
	$Workdir=~s/[\/\s\.]+$//;
	$self->{'-workdir'}=$Workdir;
}

sub getlibInput ($) {
	my $lib=shift;
	my %Input=();
	for my $i(sort keys %$lib) {
		if (exists $lib->{$i}) {
			my ($A,$B)=($1,$2) if ($i=~/(\w+)([12])/);
			foreach my $reads(@{$lib->{$i}}) {
				if (defined $A && defined $B) {
					my $C=($B==1)?2:1;
					if (exists $lib->{"$A$C"}) {
						push @{$Input{$B}},$reads;
					} else {
						push @{$Input{0}},$reads;
					}
				} else {
					push @{$Input{0}},$reads;
				}
			}
		}
	}
	return %Input;
}

sub correctJavaCmd
{
	my ($tools,$heap,$Djavaio)=@_;
	my $cmd="java";
	my $init="40m";
	if (defined $tools && $tools !~ /^java/ && $tools !~ /\-jar/)
	{
		if (defined $heap)
		{
			if ($heap=~/(\d+)(\w+)/) {
				my ($mem,$unit)=($1,$2);
				$init=int($mem/10)."$unit";
			} elsif ($heap=~/(\d+)$/) {
				$heap.="m";
				my ($mem,$unit)=($1,$2);
				$init=int($mem/10).$unit;
			}
			$init="40m" if ($init =~ /^0/);
			$cmd.=" -Xms$init -Xmx$heap";
		}
		if (defined $Djavaio)
		{
			$cmd.=" -Djava.io.tmpdir=$Djavaio";
		}
		$cmd.=" -jar $tools";
	}
	else
	{
		$cmd=$tools;
		if (defined $heap && $cmd !~ /\-Xmx/)
		{
			if ($heap=~/(\d+)(\w+)/) {
				my ($mem,$unit)=($1,$2);
				$init=int($mem/10)."unit";
			} elsif ($heap=~/(\d+)$/) {
				$heap.="m";
				my ($mem,$unit)=($1,$2);
				$init=int($mem/10).$unit;
			}
			$init="40m" if ($init =~ /^0/);
			$cmd=~s/java/java\ -Xms$init \-Xmx$heap\ /;
		}
		if (defined $Djavaio && $cmd !~ /\-Djava/)
		{
			$cmd=~s/\-jar/\-Djava\.io\.tmpdir\=$Djavaio\ \-jar/;
		}
	}
	return $cmd;
}

sub checkFileFormat ($) {
	my $file=shift;
	if (-B $file) {
		if ($file=~/bam$/i)
		{
			my $index=GATE::Error::checkIndex('samtools',$file);
			if ($index==1)
			{
				return "bam";
			}
		} elsif ($file=~/fasta.gz$/i || $file=~/fa.gz$/i) {
			return "fasta.gz";
		} elsif ($file=~/fastq.gz$/i || $file=~/fq.gz$/i) {
			return "fastq.gz";
		} elsif ($file=~/gz$/i) {
			open (IN,"zcat $file|");
			my $line=<IN>;
			if ($line=~/^\>/) {
				return "fasta.gz";
			} elsif ($line=~/^\@/) {
				return "fastq.gz";
			}
			close IN;
		} elsif ($file=~/sff/i) {
			return "sff";
		} elsif ($file=~/bfa/i) {
			return "bfa";
		} elsif ($file=~/bfq/i) {
			return "bfq";
		}
	} elsif (-T $file) {
		if ($file=~/sam/i) {
			return "sam";
		} elsif ($file=~/fasta$/i || $file=~/fa$/i || $file=~/seq/i || $file=~/fna/i || $file=~/fas/i) {
			return "fasta";
		} elsif ($file=~/fastq$/i || $file=~/fq/i) {
			return "fastq";
		} else {
			open (IN,$file);
			my $line=<IN>;
			if ($line=~/^\>/) {
				return "fasta";
			} elsif ($line=~/^\@/) {
				return "fastq";
			}
			close IN;
		}
	}
}

sub print_check_process {
	my @ary=@_;
	my $out="";
	my $process = join " ",@ary;
	my $grep="";
	foreach my $p (@ary) {
		$grep .= " \| grep $p";
	}
	$out .= qq(user=`whoami`\np=\`ps -u \$user -f $grep | grep -v grep\`\nwhile [ "\$p" != "" ]\ndo\n\techo "$process is not finish yet! sleep 120s"\n\tsleep 120\n\tp=\`ps -u \$user -f $grep | grep -v grep\`\ndone\n);
	return $out;
}

sub get_time ($) {
	my $self=shift;
	my  ($sec,$min,$hour,$mday,$mon,$year) = (localtime)[0..5];
	($sec,$min,$hour,$mday,$mon,$year) = (
		sprintf("%02d", $sec),
		sprintf("%02d", $min),
		sprintf("%02d", $hour),
		sprintf("%02d", $mday),
		sprintf("%02d", $mon + 1),
		$year + 1900
	);
	$self->{year}=$year;
	$self->{month}=$mon;
	$self->{day}=$mday;
	$self->{hour}=$hour;
	$self->{minute}=$min;
	$self->{second}=$sec;
	return "## $year-$mon-$mday $hour:$min:$sec\n";
}

#########################################################
#                                                       #
#                    Data Processing                    #
#                                                       #
#########################################################

#phred_020425
sub runPhred ($) {
	my $self = shift;
	if (!exists $self->{"software:phred"} && `which phred` eq "") {
		return "";
	}
	my $phred=(exists $self->{"software:phred"})?checkPaht($self->{"software:phred"}):`which phred`;
	chomp $phred;
	my $phredpar=$self->{"database:phredpar"} if (exists $self->{"database:phredpar"});
	my $phred_cmd =  qq(echo `date`; echo "run Phred"\n);
	$phred_cmd .= qq(export PHRED_PARAMETER_FILE=$phredpar\n) if (defined $phredpar);
	$phred_cmd .= qq(export phred=$phred\n) if (defined $phred);
	my $para = (exists $self->{"setting:phred"}) ? $self->{"setting:phred"} : qq(-trim_cutoff 0.01 -trim_alt \\'\\' -trim_fasta);
	$phred_cmd .= qq(export phredpara="$para"\n) if (defined $para);
	$phred_cmd .= qq(cd $self->{"-workdir"}\n);
	$phred_cmd .= qq(export bc_outdir=$self->{"setting:bc_outdir"}\n);
	$phred_cmd .= qq([[ -d \${bc_outdir} || mkdir \${bc_outdir}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:bc_outdir"}));
	$phred_cmd .= qq(cd \${bc_outdir}\n);
	$phred_cmd .= qq([[ -d phred ]] || mkdir phred\n);
	$phred_cmd .= qq(cd phred\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	my $gotseq=0;
	foreach my $lib(@libraries) {
		my @dir;
		if (exists $self->{"LIB"}{$lib}{"sanger"}) {
			push @dir, ${$self->{"LIB"}{$lib}{"sanger"}}[0];
		} elsif (exists $self->{"LIB"}{$lib}{"3730"}) {
			push @dir, ${$self->{"LIB"}{$lib}{"3730"}}[0];
		} elsif (exists $self->{"LIB"}{$lib}{"ab1"}) {
			push @dir, ${$self->{"LIB"}{$lib}{"3730"}}[0];
		}
		if (@dir > 0) {
			for my $directory (@dir) {
				my $dirname = (split /\//,$directory)[-1]; 
				$phred_cmd .= qq([[ -d $lib ]] || mkdir $lib\n);
				$phred_cmd .= qq(cd $lib\n);
				$phred_cmd .= qq(ln -s $directory ; mkdir -p $dirname/qual_dir ; mkdir -p $dirname/seq_dir\n);
				$phred_cmd .= qq(\${phred} \${phredpara}-id $dirname/ -sd $dirname/seq_dir/ -qd $dirname/qual_dir/\n);
				$phred_cmd .= qq(cat $dirname/seq_dir/*.seq > $lib.$dirname.seq\n);
				$phred_cmd .= qq(cat $dirname/qual_dir/*.qual > $lib.$dirname.qual\n);
				push @{$self->{"LIB"}{$lib}{"seq"}},qq($self->{"-workdir"}/$self->{"setting:bc_outdir"}/$lib/$lib.$dirname.seq);
				push @{$self->{"LIB"}{$lib}{"qual"}},qq($self->{"-workdir"}/$self->{"setting:bc_outdir"}/$lib/$lib.$dirname.qual);
				$phred_cmd .= qq(cd ..\n);
				$gotseq++;
			}
		}
	}
	$phred_cmd .= qq(cd ..\ncd ..\n);
	return $phred_cmd if ($gotseq!=0);
}

#CASAVA_1.8.2
sub runCASAVA ($) {
	my $self = shift;
	my $CASAVA_PATH="";
	if (exists $self->{'setting:CASAVA_PATH'}) {
		$CASAVA_PATH=GATE::Error::checkPath($self->{'setting:CASAVA_PATH'});
	} elsif (exists $self->{'setting:CASAVAPATH'}){
		$CASAVA_PATH=GATE::Error::checkPath($self->{'setting:CASAVAPATH'});
	} elsif (exists $self->{'software:CASAVA'}){
		$CASAVA_PATH=GATE::Error::checkPath($1) if (/(\S+)\/bin/);
	}
	my $para = (exists  $self->{'setting:CASAVA'}) ? $self->{'setting:CASAVA'} : '--fastq-cluster-count 0 --mismatches 1'; 
	my $casava_cmd = qq(echo `date`; echo "run CASAVA"\n);
	$casava_cmd .= qq(export CASAVA_PATH=$CASAVA_PATH\n) if (defined $CASAVA_PATH);
	$casava_cmd .= qq(export CASAVA_para="$para"\n);
	$casava_cmd .= qq([[ -d \${bc_outdir} || mkdir \${bc_outdir}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:bc_outdir"}));
	$casava_cmd .= qq(cd \${bc_outdir}\n);
	$casava_cmd .= qq([[ -d casava ]] || mkdir casava\n);
	$casava_cmd .= qq(cd casava\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	my $gotseq=0;
	foreach my $lib(@libraries) {
		$casava_cmd .= qq([[ -d $lib ]] || mkdir $lib\n);
		$casava_cmd .= qq(cd $lib\n);
		my @input=@{$self->{'LIB'}{$lib}{"input-dir"}};
		for (my $i=0;$i<@input;$i++) {
			my $name;
			if ($input[$i] =~ /\/([^\/]+)\/Data/) {
				$name = $1;
			} else {
				$name = "$lib.".($i+1);
			}
			$casava_cmd .= qq(mkdir $name\n);
			$casava_cmd .= qq(\${CASAVA_PATH}/bin/configureBclToFastq.pl --input-dir $input[$i] \\\n);
			$casava_cmd .= qq(--outdir-dir ./$name/Demultiplexed --sample-sheet ./$name/sample_sheet/$name.csv \\\n);
			$casava_cmd .= qq(\${CASAVA_para}\n);
			$gotseq++;
		}
		$casava_cmd .= qq(cd ..\n);
	}
	$casava_cmd .= qq(cd ..\n);
	$casava_cmd .= qq(cd ..\n);
	return $casava_cmd if ($gotseq>0);
}

sub runLifeScope($) {
	
}

sub runSolidCaller ($) {
	
}

sub runDataAnalysis ($) {
	
}

sub splitSequence ($) {
	my $self=shift;
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		
	}
}

#selectIdxFastq.pl v1.1, 2013-05-03
sub selectIdxFastq ($) {
	my $self = shift;
	if (!exists $self->{"software:selectIdxFastq"} || !defined $self->{"software:selectIdxFastq"} || (!exists $self->{'idx'} && !exists $self->{'bar'})) {
		return "";
	}
	my $selectIdxFastq=GATE::Error::checkPath($self->{"software:selectIdxFastq"});
	my $para=(exists $self->{"setting:selectIdxFastq"})?$self->{"setting:selectIdxFastq"}:"-mis 1 -qual 30";
	$para=~s/\-.+prefix\s+\S+//;
	my $multirun=GATE::Error::checkPath($self->{"software:multithreads-run"}) if (exists $self->{"software:multithreads-run"});
	my $multirun_sh=(-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"})) ? qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/selectIdxFastq.$$.sh) : qq($self->{"-workdir"}/selectIdxFastq.$$.sh);
	my $Idx_cmd = qq(echo `date`; echo "run selectIdxFastq"\n);
	if (defined $multirun) {
		open (IDX,">$multirun_sh");
		$Idx_cmd .= qq(export multirun="$multirun"\n);
	}
	$Idx_cmd .= qq(cd $self->{"-workdir"}\n);
	$Idx_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$Idx_cmd .= qq(export selectIdxFastq=$selectIdxFastq\n);
	$Idx_cmd .= qq(export selectIdxFastq_para="$para"\n);
	$Idx_cmd .= qq(mkdir $self->{"setting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}));
	$Idx_cmd .= qq(cd $self->{"setting:qc_outdir"}\n);
	my $Idx_cmd_multi_head=$Idx_cmd;
	my $Idx_cmd_multi = "";
	my $withIdx=0;
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		next if (!exists $self->{'idx'}{$lib} && !exists $self->{'bar'}{$lib});
		$Idx_cmd_multi .= $Idx_cmd_multi_head;
		$Idx_cmd .= qq(echo `date` && echo "$lib"\n);
		$Idx_cmd .= qq([[ -d $lib ]] || mkdir $lib\n)  unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
		$Idx_cmd_multi .= qq(echo `date` && echo "$lib" && [[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
		if (defined $multirun) {
			print IDX qq(echo `date` && echo "$lib" && );
			print IDX qq([[ -d $lib ]] || mkdir $lib && cd $lib && ) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
		}
		$Idx_cmd .= "cd $lib\n";
		$Idx_cmd_multi .= "cd $lib\n";
		for my $lbmark(sort keys %{$self->{"LIB"}{$lib}}) {
			my ($lb,$i)=($1,$2) if ($lbmark=~/^(\w+)([12])/);
			$lb=$lbmark if (!defined $lb);
			my $j=1;
			if (defined $i && $i == 2) {
				next if (exists $self->{"LIB"}{$lib}{"$lb$i"});
			}
			my (@fq1,@fq2,@fq3)=();
			if (defined $i && $i == 1)
			{
				$j=2;
				if (exists $self->{"LIB"}{$lib}{"$lb$j"})
				{
					foreach my $reads(@{$self->{"LIB"}{$lib}{$lbmark}})
					{
						push @fq1,$reads;
					}
					foreach my $reads(@{$self->{"LIB"}{$lib}{"$lb$j"}})
					{
						push @fq2,$reads;
					}
				}
				else
				{
					foreach my $reads(@{$self->{"LIB"}{$lib}{$lbmark}})
					{
						push @fq3,$reads;
					}
				}
			}
			else
			{
				foreach my $reads(@{$self->{"LIB"}{$lib}{$lbmark}})
				{
					push @fq3,$reads;
				}
			}
			if (@fq2>0 && @fq1==@fq2){
				for (my $k=0;$k<@fq2;$k++){
					my $reads1=$fq1[$k];
					my $reads2=$fq2[$k];
					if ( (exists $self->{$lib}{$lbmark}{$k}{"Index"} && $self->{$lib}{$lbmark}{$k}{"Index"} ne "") ||
					    (exists $self->{$lib}{$lbmark}{$k}{"Barcode"} && $self->{$lib}{$lbmark}{$k}{"Barcode"} ne "") ){
						my $index=$self->{$lib}{$lbmark}{$k}{"Index"} if (exists $self->{$lib}{$lbmark}{$k}{"Index"});
						my $barcode=$self->{$lib}{$lbmark}{$k}{"Barcode"} if (exists $self->{$lib}{$lbmark}{$k}{"Barcode"});
						my $name1=(split /\//,$reads1)[-1];
						my $name2=(split /\//,$reads2)[-1];
						my $out1=$1 if ($name1=~/(\S+)\.fastq$/i || $name1=~/(\S+)\.fq$/i || $name1=~/(\S+)\.fastq\.gz$/i || $name1=~/(\S+)\.fq\.gz$/i);
						$out1.="\.$index" if (defined $index);
						$out1.="\.$1" if (defined $barcode && $barcode=~/^([^\:\s])+/);
						$out1.="\.fastq";
						my $out2=$1 if ($name2=~/(\S+)\.fastq$/i || $name2=~/(\S+)\.fq$/i || $name2=~/(\S+)\.fastq\.gz$/i || $name2=~/(\S+)\.fq\.gz$/i);
						$out2.="\.$index" if (defined $index);
						$out2.="\.$1" if (defined $barcode && $barcode=~/([^\:\s])+$/);
						$out2.="\.fastq";
						$Idx_cmd .= qq(\${selectIdxFastq} \${selectIdxFastq_para} -fastq1 $reads1 -fastq2 $reads2 );
						$Idx_cmd .= qq( -index $index ) if (defined $index);
						$Idx_cmd .= qq( -barcode $barcode ) if (defined $barcode);
						$Idx_cmd .= qq(\n);
						$Idx_cmd_multi .= qq(\${selectIdxFastq} \${selectIdxFastq_para} -fastq1 $reads1 -fastq2 $reads2);
						$Idx_cmd_multi .= qq( -index $index ) if (defined $index);
						$Idx_cmd_multi .= qq( -barcode $barcode ) if (defined $barcode);
						$Idx_cmd_multi .= qq( && );
						if (defined $multirun) {
							print IDX qq($selectIdxFastq $para -fastq1 $reads1 -fastq2 $reads2);
							print IDX qq( -index $index ) if (defined $index);
							print IDX qq( -barcode $barcode ) if (defined $barcode);
							print IDX qq( && );
						}
						die qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out2 is existent!\n) if (-f qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out1));
						die qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out2 is existent!\n) if (-f qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out2));
						${$self->{"LIB"}{$lib}{$lbmark}}[$k]=qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out1);
						${$self->{"LIB"}{$lib}{"$lb$j"}}[$k]=qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out2);
						$withIdx++;
					}
				}
			} else {
				if (@fq2>0 && @fq1!=@fq2)
				{
					push @fq3,@fq1;
					push @fq3,@fq2;
				}
				if (@fq3>0)
				{
					for (my $k=0;$k<@fq3;$k++){
						my $reads1=$fq3[$k];
						if ( (exists $self->{$lib}{$lbmark}{$k}{"Index"} && $self->{$lib}{$lbmark}{$k}{"Index"} ne "") ||
						    (exists $self->{$lib}{$lbmark}{$k}{"Barcode"} && $self->{$lib}{$lbmark}{$k}{"Barcode"} ne "") ){
							my $index=$self->{$lib}{$lbmark}{$k}{"Index"} if (exists $self->{$lib}{$lbmark}{$k}{"Index"});
							my $barcode=$self->{$lib}{$lbmark}{$k}{"Barcode"} if (exists $self->{$lib}{$lbmark}{$k}{"Barcode"});
							my $out1=(split /\//,$reads1)[-1];
							$out1="$1" if ($out1=~/(\S+)\.fastq$/i || $out1=~/(\S+)\.fq$/i || $out1=~/(\S+)\.fastq\.gz$/i || $out1=~/(\S+)\.fq\.gz$/i);
							$out1.="\.$index" if (defined $index);
							$out1.="\.$1" if (defined $barcode && $barcode=~/^([^\:\s])+/);
							$out1.="\.fastq";
							$Idx_cmd .= qq(\${selectIdxFastq} \${selectIdxFastq_para}  -fastq $reads1 );
							$Idx_cmd .= qq( -index $index ) if (defined $index);
							$Idx_cmd .= qq( -barcode $barcode ) if (defined $barcode);
							$Idx_cmd .= qq(\n);
							$Idx_cmd_multi .= qq(\${selectIdxFastq} \${selectIdxFastq_para} -fastq $reads1 );
							$Idx_cmd_multi .= qq( -index $index ) if (defined $index);
							$Idx_cmd_multi .= qq( -barcode $barcode ) if (defined $barcode);
							$Idx_cmd_multi .= qq( && );
							if (defined $multirun) {
								print IDX qq($selectIdxFastq $para -fastq );
								print IDX qq( -index $index ) if (defined $index);
								print IDX qq( -barcode $barcode ) if (defined $barcode);
								print IDX qq( && );
							}
							die qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out1 is existent!\n) if (-f qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out1));
							${$self->{"LIB"}{$lib}{$lbmark}}[$k]=qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out1);
							$withIdx++;
						}
					}
				}
			}
		}
		$Idx_cmd .= "cd ..\n";
		if (defined $Idx_cmd_multi){
			$Idx_cmd_multi =~ s/\s+$//;
			$Idx_cmd_multi =~ s/\&\&$//;
		}
		if ($withIdx % $self->{"setting:multithreads"}!=0)
		{
			$Idx_cmd_multi .= " &\n";
		} else {
			$Idx_cmd_multi .= "\n";
		}
		$Idx_cmd_multi .= "cd ..\ncd ..\n";
		if (defined $multirun) {
			print IDX "cd ..\n";
		}
	}
	if (defined $multirun){
		close IDX;
	}
	$Idx_cmd .= "cd ..\n";
	if ($withIdx>0){
		if (exists $self->{"rule:multimode"}) {
			chomp $Idx_cmd_multi;
			$Idx_cmd_multi=~s/\&+$//;
			$Idx_cmd_multi.="\n";
			$Idx_cmd_multi.=print_check_process('selectIdxFastq');
			if (defined $multirun) {
				my $multirun_cmd=$Idx_cmd_multi_head;
				$multirun_cmd .= qq(${multirun} $multirun_sh -nt $self->{"setting:multithreads"}\n);
				$multirun_cmd .= print_check_process('selectIdxFastq');
				return $multirun_cmd;
			} else {
				return $Idx_cmd_multi;
			}
		}else{
			return $Idx_cmd;
		}
	}
}

#mergeOverlapPE.pl v0.1.3 alpha 2013-03-20
sub mergeOverlapPE($) {
	my $self = shift;
	if ( (!exists $self->{"software:mergeOverlapPE"} && (!exists $self->{"software:bwa-pemerge"})) || (!exists $self->{'pemerge'}) ) {
		return "";
	}
	if (exists $self->{"software:mergeOverlapPE"}) {
		my $mergeOverlapPE=$self->{"software:mergeOverlapPE"};
		my $moppara = $self->{"setting:mergeOverlapPE"};
		my $mop_cmd = qq(echo `date`; echo "run mergeOverlapPE"\n);
		my $multirun=GATE::Error::checkPath($self->{"software:multithreads-run"}) if (exists $self->{"software:multithreads-run"});
		my $multirun_sh=(-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"})) ? qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/mergeOverlapPE.$$.sh) : qq($self->{"-workdir"}/mergeOverlapPE.$$.sh);
		if (defined $multirun) {
			open (MOP,">$multirun_sh");
			$mop_cmd .= qq(export multirun="$multirun"\n);
		}
		$mop_cmd .= qq(cd $self->{"-workdir"}\n);
		$mop_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
		$mop_cmd .= qq(mkdir $self->{"setting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}));
		$mop_cmd .= qq(cd $self->{"setting:qc_outdir"}\n);
		my $mop_cmd_multi_head=$mop_cmd;
		my $mop_cmd_multi="";
		my $withPE=0;
		my @libraries=sort keys %{$self->{'LIB'}};
		foreach my $lib(@libraries) {
			next if (!exists $self->{'pemerge'}{$lib});
			$mop_cmd_multi.= $mop_cmd_multi_head;
			$mop_cmd_multi .= qq(echo `date`; echo "$lib"\n);
			$mop_cmd .= qq(echo `date`; echo "$lib"\n);
			$mop_cmd .= qq([[ -d $lib ]] || mkdir $lib\n)  unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
			$mop_cmd_multi .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
			if (defined $multirun) {
				print MOP qq(echo `date` && echo "$lib" && );
				print MOP qq([[ -d $lib ]] || mkdir $lib && ) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
				print MOP qq(cd $lib && );
			}
			$mop_cmd .= "cd $lib\n";
			$mop_cmd_multi .= "cd $lib && ";
			my %fq=getlibInput($self->{"LIB"}{$lib});
			my @Reads1=();
			@Reads1=@{$fq{1}} if (exists $fq{1});
			my @Reads2=();
			@Reads2=@{$fq{2}} if (exists $fq{2});
			if (@Reads1==0 || @Reads2==0)
			{
				$mop_cmd .= "cd ../\n";
				$mop_cmd_multi .= "cd ../";
				if ($withPE % $self->{"setting:multithreads"}!=0)
				{
					$mop_cmd_multi .= " &\n";
				} else {
					$mop_cmd_multi .= "\n";
				}
				next;
			}
			for (my $i=0;$i<@Reads1;$i++)
			{
				if (exists $self->{$lib}{'fq1'}{$i}{"MergePE"} && GATE::Error::boolean($self->{$lib}{'fq1'}{$i}{"MergePE"})==1 ) {
					my $prefix=(@Reads1>1)?"$lib-".($i+1):$lib;
					$mop_cmd .= "$mergeOverlapPE $Reads1[$i] $Reads2[$i] -prefix $prefix $moppara\n";
					$mop_cmd_multi .= "$mergeOverlapPE $Reads1[$i] $Reads2[$i] -prefix $prefix $moppara && ";
					if (defined $multirun) {
						print MOP "$mergeOverlapPE $Reads1[$i] $Reads2[$i] -prefix $prefix $moppara && ";
					}
					${$self->{"LIB"}{$lib}{'fq1'}}[$i]=qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$prefix\_R1.fastq);
					${$self->{"LIB"}{$lib}{'fq2'}}[$i]=qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$prefix\_R2.fastq);
					push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$prefix\_merged.fastq);
					foreach my $k(keys %{$self->{$lib}{'fq1'}{$i}})
					{
						$self->{$lib}{'fq'}{@{$self->{"LIB"}{$lib}{'fq'}}-1}{$k}=$self->{$lib}{'fq1'}{$i}{$k};
					}
					$withPE++;
					if ($withPE % $self->{"setting:multithreads"} != 0)
					{
						$mop_cmd_multi .= " &\n";
					} else {
						$mop_cmd_multi .= "\n";
					}
				}
			}
			$mop_cmd .= "cd ..\n";
			$mop_cmd_multi .= "cd ..\ncd ..\n";
			if (defined $multirun) {
				print MOP "cd ..\n";
			}
		}
		$mop_cmd .= "cd ..\n";
		close MOP if (defined $multirun);
		if ($withPE>0){
			if (exists $self->{"rule:multimode"}) {
				$mop_cmd_multi =~ s/\s+$//;
				$mop_cmd_multi =~ s/\&+$//;
				$mop_cmd_multi .= "\n";
				$mop_cmd_multi .= print_check_process('mergeOverlapPE');
				if (defined $multirun) {
					my $multirun_cmd = $mop_cmd_multi_head;
					$multirun_cmd .= qq(${multirun} $multirun_sh -nt $self->{"setting:multithreads"}\n);
					$multirun_cmd .= print_check_process('mergeOverlapPE');
					return $multirun_cmd;
				} else {
					return $mop_cmd_multi;
				}
			}else{
				return $mop_cmd;
			}
		}
	}elsif (exists $self->{"software:bwa-pemerge"}) {
		my $bwapemerge=$self->{"software:bwa-pemerge"};
		my $moppara = (exists $self->{"rule:bwa-pemerge"}) ? $self->{"rule:bwa-pemerge"} : "";
		my $mop_cmd = qq(echo `date`; echo "run bwa-pemerge"\n);
		$mop_cmd .= qq(cd $self->{"-workdir"}\n);
		$mop_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
		$mop_cmd .= qq(mkdir $self->{"setting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}));
		$mop_cmd .= qq(cd $self->{"setting:qc_outdir"}\n);
		my $mop_cmd_multi_head=$mop_cmd;
		my $mop_cmd_multi="";
		my $withPE=0;
		my @libraries=sort keys %{$self->{'LIB'}};
		foreach my $lib(@libraries) {
			next if (!exists $self->{'pemerge'}{$lib});
			$mop_cmd_multi.= $mop_cmd_multi_head;
			$mop_cmd_multi .= qq(echo `date`; echo "$lib"\n);
			$mop_cmd .= qq(echo `date`; echo "$lib"\n);
			$mop_cmd .= qq([[ -d $lib ]] || mkdir $lib\n)  unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
			$mop_cmd_multi .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
			$mop_cmd .= "cd $lib\n";
			$mop_cmd_multi .= "cd $lib\n";
			my %fq=getlibInput($self->{"LIB"}{$lib});
			my @Reads1=();
			@Reads1=@{$fq{1}} if (exists $fq{1});
			my @Reads2=();
			@Reads2=@{$fq{2}} if (exists $fq{2});
			if (@Reads1==0 || @Reads2==0)
			{
				$mop_cmd .= "cd ../\n";
				$mop_cmd_multi .= "cd ../";
				if ($withPE % $self->{"setting:multithreads"}!=0)
				{
					$mop_cmd_multi .= " &\n";
				} else {
					$mop_cmd_multi .= "\n";
				}
				next;
			}
			for (my $i=0;$i<@Reads1;$i++)
			{
				if ( exists $self->{$lib}{'fq1'}{$i}{"MergePE"} && GATE::Error::boolean($self->{$lib}{'fq1'}{$i}{"MergePE"})==1 ) {
					my $prefix=(@Reads1>1)?"$lib-".($i+1):$lib;
					$mop_cmd .= "$bwapemerge $moppara -m $Reads1[$i] $Reads2[$i] > $prefix\_merged.fastq\n";
					$mop_cmd_multi .= qq($bwapemerge $moppara -m -t $self->{"setting:multithreads"} $Reads1[$i] $Reads2[$i] > $prefix\_merged.fastq && );
					$mop_cmd .= qq($bwapemerge $moppara -u $Reads1[$i] $Reads2[$i] | awk '\{if \(NR\%8<4\)\{print \$0 > "$prefix\_R1.fastq"\}else\{print \$0 > "$prefix\_R2.fastq"\}\}'\n);
					$mop_cmd_multi .= qq($bwapemerge $moppara -u -t $self->{"setting:multithreads"} $Reads1[$i] $Reads2[$i] | awk '\{if \(NR\%8<4\)\{print \$0 > "$prefix\_R1.fastq"\}else\{print \$0 > "$prefix\_R2.fastq"\}\}' && );
					${$self->{"LIB"}{$lib}{'fq1'}}[$i]=qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$prefix\_R1.fastq);
					${$self->{"LIB"}{$lib}{'fq2'}}[$i]=qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$prefix\_R2.fastq);
					push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$prefix\_merged.fastq);
					foreach my $k(keys %{$self->{$lib}{'fq1'}{$i}})
					{
						$self->{$lib}{'fq'}{@{$self->{"LIB"}{$lib}{'fq'}}-1}{$k}=$self->{$lib}{'fq1'}{$i}{$k};
					}
					$withPE++;
				}
			}
			$mop_cmd .= "cd ..\n";
			$mop_cmd_multi =~ s/\s+$//;
			$mop_cmd_multi =~ s/\&\&$//;
			if ($withPE % $self->{"setting:multithreads"} != 0)
			{
				$mop_cmd_multi .= " &\n";
			} else {
				$mop_cmd_multi .= "\n";
			}
			$mop_cmd_multi .= "cd ..\ncd ..\n";
		}
		$mop_cmd .= "cd ..\n";
		if ($withPE>0){
			if (exists $self->{"rule:multimode"}) {
				chomp $mop_cmd_multi;
				$mop_cmd_multi =~ s/\s+\&+$//;
				$mop_cmd_multi .= "\n";
				$mop_cmd_multi .= print_check_process('mergeOverlapPE');
				return $mop_cmd_multi;
			}else{
				return $mop_cmd;
			}
		}
	}
}

sub runQA($) {
	my $self = shift;
	if ( (!exists $self->{"software:SolexaQA"} || !defined $self->{"software:SolexaQA"}) &&
	    (!exists $self->{"software:check_fastq"} || !-e $self->{"software:check_fastq"}) ){
		return "";
	}
	if ( (exists $self->{"software:SolexaQA"} || -e $self->{"software:SolexaQA"}) ) {
		my $SolexaQA=GATE::Error::checkPath($self->{"software:SolexaQA"});
		my $para=$self->{"setting:SolexaQA"};
		#SolexaQA.pl ../BF3.1.fq -v -m -s 10000 -b -sanger -d ./ 
		my $qa_cmd = qq(echo `date`; echo "run SolexaQA"\n);
		$qa_cmd .= qq(cd $self->{"-workdir"}\n);
		$qa_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
		$qa_cmd .= qq(mkdir $self->{"setting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}));
		$qa_cmd .= qq(cd $self->{"setting:qc_outdir"}\n);
		my $qa_cmd_multi_head=$qa_cmd;
		my $qa_cmd_multi="";
		my $multi=0;
		my @libraries=sort keys %{$self->{'LIB'}};
		foreach my $lib(@libraries){
			$qa_cmd_multi.= $qa_cmd_multi_head;
			$qa_cmd_multi .= qq(echo `date`; echo "$lib"\n);
			$qa_cmd .= qq(echo `date`; echo "$lib"\n);
			$qa_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
			$qa_cmd_multi .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
			$qa_cmd .= "cd $lib\n";
			$qa_cmd .= "mkdir QA\n" if ($para !~ /\-d/);
			$qa_cmd_multi .= "cd $lib\n";
			$qa_cmd_multi .= "mkdir QA && " if ($para !~ /\-d/);
			for my $i(sort keys %{$self->{"LIB"}{$lib}}) {
				if (exists $self->{"LIB"}{$lib}{$i}) {
					my $j=0;
					foreach my $reads(@{$self->{"LIB"}{$lib}{$i}}) {
						if (($reads !~ /fastq$/i && $reads !~ /fq$/i) || ($reads =~ /gz$/)) {
							my $fq=(split /\//,$reads)[-1];
							$fq=$1 if ($fq=~/(\S+)\.gz/);
							if (!-f qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/QA/$fq))
							{
								$qa_cmd .= "gzip -cd $reads > $fq\n";
								$qa_cmd_multi .= "gzip -cd $fq && ";
							}
							$qa_cmd .= "$SolexaQA $fq ";
							$qa_cmd .= ($para=~/\-d/) ? " $para\n" : "$para -d QA\n";
							$qa_cmd_multi .= "$SolexaQA $fq ";
							$qa_cmd_multi .= ($para=~/\-d/) ? " $para && " : "$para -d QA && ";
						} else {
							$qa_cmd .= "ln -s $reads\n";
							$qa_cmd .= "$SolexaQA $reads ";
							$qa_cmd .= ($para=~/\-d/) ? "$para\n" : "$para -d QA\n";
							$qa_cmd_multi .= "ln -s $reads && ";
							$qa_cmd_multi .= "$SolexaQA $reads ";
							$qa_cmd_multi .= ($para=~/\-d/) ? " $para && " : "$para -d QA && ";
						}
					}
				}
			}
			$qa_cmd .= "cd ..\n";
			if ($multi % $self->{"setting:multithreads"}!=0){
				$qa_cmd_multi .= "&\n";
			}else{
				$qa_cmd_multi .= "\n";
			}
			$qa_cmd_multi .= "cd ..\ncd ..\n";
			$multi++;
		}
		$qa_cmd .= "cd ..\n";
		if (exists $self->{"rule:multimode"} && $self->{"rule:multimode"} =~ /y/i) {
			chomp $qa_cmd_multi;
			$qa_cmd_multi=~s/\&+$//;
			$qa_cmd_multi.="\n";
			return $qa_cmd_multi;
		} else {
			return ($qa_cmd);
		}
	}
	if ( (exists $self->{"software:check_fastq"} || -e $self->{"software:check_fastq"}) )
	{
		my $check_fastq=GATE::Error::checkPath($self->{"software:check_fastq"});
		my $para=$self->{"setting:check_fastq"};
		my $distribute_fqcheck=$self->{"setting:distribute_fqcheck"} if (exists $self->{"setting:distribute_fqcheck"});
		#SolexaQA.pl ../BF3.1.fq -v -m -s 10000 -b -sanger -d ./ 
		my $qa_cmd = qq(echo `date`; echo "run check_fastq"\n);
		$qa_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
		$qa_cmd .= qq(export check_fastq="$check_fastq"\n);
		$qa_cmd .= qq(cd $self->{"-workdir"}\n);
		$qa_cmd .= qq(mkdir $self->{"setting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}));
		$qa_cmd .= qq(cd $self->{"setting:qc_outdir"}\n);
		my $qa_cmd_multi_head=$qa_cmd;
		my $qa_cmd_multi="";
		my $multi=0;
		my @libraries=sort keys %{$self->{'LIB'}};
		foreach my $lib(@libraries){
			$qa_cmd_multi.= $qa_cmd_multi_head;
			$qa_cmd_multi .= qq(echo `date`; echo "$lib"\n);
			$qa_cmd .= qq(echo `date`; echo "$lib"\n);
			$qa_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
			$qa_cmd_multi .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
			$qa_cmd .= "cd $lib\n";
			$qa_cmd .= "mkdir QA\n" if ($para !~ /\-d/);
			$qa_cmd_multi .= "cd $lib && ";
			$qa_cmd_multi .= "mkdir QA && " if ($para !~ /\-d/);
			for my $i(sort keys %{$self->{"LIB"}{$lib}}) {
				if (exists $self->{"LIB"}{$lib}{$i}) {
					my $j=0;
					foreach my $reads(@{$self->{"LIB"}{$lib}{$i}}) {
						if (($reads !~ /fastq$/i && $reads !~ /fq$/i) || ($reads =~ /gz$/)) {
							my $fq=(split /\//,$reads)[-1];
							$fq=$1 if ($fq=~/(\S+)\.gz/);
							if (!-f qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/QA/$fq))
							{
								$qa_cmd .= "gzip -cd $reads > $fq\n";
								$qa_cmd_multi .= "gzip -cd $fq && ";
							}
							$qa_cmd .= "\${check_fastq} -i $fq ";
							$qa_cmd .= (defined $para) ? " $para > $fq.fqcheck\n" : " > $fq.fqcheck\n";
							$qa_cmd_multi .= "\${check_fastq} -i $fq ";
							$qa_cmd_multi .= (defined $para) ? " $para > $fq.fqcheck && " : " > $fq.fqcheck && ";
							if (defined $distribute_fqcheck)
							{
								$qa_cmd .= qq(\${distribute_fqcheck} $fq.fqcheck -o $fq.fqcheck\n);
								$qa_cmd_multi .= qq(\${distribute_fqcheck} $fq.fqcheck -o $fq.fqcheck && );
							}
						}
					}
				}
			}
			$qa_cmd .= "cd ..\n";
			$qa_cmd_multi .= "cd .. && ";
			$qa_cmd_multi .= "cd .. ";
			$multi++;
			if ($multi % $self->{"setting:multithreads"}!=0){
				$qa_cmd_multi .= "&\n";
			}else{
				$qa_cmd_multi .= "\n";
			}
		}
		$qa_cmd .= "cd ..\n";
		if (exists $self->{"rule:multimode"} && $self->{"rule:multimode"} =~ /y/i) {
			chomp $qa_cmd_multi;
			$qa_cmd_multi=~s/\s\&+$//;
			$qa_cmd_multi.="\n";
			return $qa_cmd_multi;
		} else {
			return ($qa_cmd);
		}
	}
}

sub runFltDup ($) {
	my $self = shift;
	if (!exists $self->{"software:filterPCRdup"} || !defined $self->{"software:filterPCRdup"}){
		return "";
	}
	my $filterPCRdup=GATE::Error::checkPath($self->{"software:filterPCRdup"});
	my $para=(exists $self->{"setting:filterPCRdup"})? $self->{"setting:filterPCRdup"} : "";
	my $fltpcr_cmd = qq(echo `date`; echo "run filterPCRdup"\n);
	$fltpcr_cmd .= qq(cd $self->{"-workdir"}\n);
	$fltpcr_cmd .= qq(mkdir $self->{"setting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}));
	$fltpcr_cmd .= qq(cd $self->{"setting:qc_outdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib (@libraries) {
		$fltpcr_cmd .= qq(echo `date`; echo "$lib"\n);
		$fltpcr_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
		$fltpcr_cmd .= "cd $lib\n";
		my %fq=getlibInput($self->{"LIB"}{$lib});
		my @Reads1=();
		@Reads1=@{$fq{1}} if (exists $fq{1});
		my @Reads2=();
		@Reads2=@{$fq{2}} if (exists $fq{2});
		my @Reads3=();
		push @Reads3,@{$fq{0}} if (exists $fq{0});
		if ((@Reads1==0 || @Reads2==0) && @Reads3==0)
		{
			$fltpcr_cmd .= "cd ..\n";
			next;
		}
		for (my $j=0;$j<@Reads1;$j++) {
			my $reads1=$Reads1[$j];
			my $reads2=$Reads2[$j] if (@Reads2>0);
			my $fqname1=(split /\//,$reads1)[-1];
			my $fq1=$reads1;
			if (($reads1 !~ /fastq$/i && $reads1 !~ /fq$/i) || ($reads1 =~ /gz$/)) {
				$fq1=$1 if ($fqname1=~/(\S+)\.gz/);
				if (!-f qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$fq1))
				{
					$fltpcr_cmd .= "gzip -cd $reads1 > $fq1\n" if ($reads1=~/gz$/);
				}
			}
			my $fqname2=(split /\//,$reads2)[-1] if (defined $reads2);
			my $fq2=$reads2 if (defined $reads2);
			if ( (defined $reads2) && ( ($reads2 !~ /fastq$/i && $reads2 !~ /fq$/i) || ($reads2 =~ /gz$/) ) ) {
				$fq2=$1 if ($fqname2=~/(\S+)\.gz/);
				if (!-f qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$fq2))
				{
					$fltpcr_cmd .= "gzip -cd $reads2 > $fq2\n" if ($reads2=~/gz$/);
				}
			}
			my $prefix="";
			my ($out1,$out2);
			if ($para eq "" || $para !~ /prefix/) {
				if (@Reads1>1) {
					$prefix="-prefix=$lib\_".($j+1);
					$out1="$lib\_".($j+1)."_uniq1.fastq";
				} elsif(@Reads1>0) {
					$prefix="-prefix=$lib";
					$out1="$lib\_uniq1.fastq";
				}
				if (@Reads2>1 || @Reads1>1) {
					$out2="$lib\_".($j+1)."_uniq2.fastq";
				} elsif(@Reads2>0) {
					$out2="$lib\_uniq2.fastq";
				}
			} elsif ($para =~ /prefix\=(\S+)/) {
				if (@Reads1>1) {
					$out1="$1\_".($j+1)."_uniq1.fastq";
				} elsif(@Reads1>0) {
					$out1="$1\_uniq1.fastq";
				}
				if (@Reads2>1 || @Reads1>1) {
					$out2="$1\_".($j+1)."_uniq2.fastq";
				} elsif(@Reads2>0)  {
					$out2="$1\_uniq2.fastq";
				}
			}
			$fltpcr_cmd .= "$filterPCRdup -fastq1=$fq1 ";
			$fltpcr_cmd .= " -fastq2=$fq2 " if (defined $reads2);
			$fltpcr_cmd .= " $prefix $para\n";
			if (@Reads1>0 && @Reads2>0) 
			{
				${$self->{"LIB"}{$lib}{'fq1'}}[$j]=qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out1);
				${$self->{"LIB"}{$lib}{'fq2'}}[$j]=qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out2);
			}
			elsif (@Reads1>0)
			{
				${$self->{"LIB"}{$lib}{'fq1'}}[$j]=qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out1);
			}
		}
		for (my $j=0;$j<@Reads3;$j++) {
			my $reads3=$Reads3[$j];
			my $fqname3=(split /\//,$reads3)[-1];
			my $fq3=$reads3;
			if (($reads3 !~ /fastq$/i && $reads3 !~ /fq$/i) || ($reads3 =~ /gz$/)) {
				$fq3=$1 if ($fqname3=~/(\S+)\.gz/);
				if (!-f qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$fq3))
				{
					$fltpcr_cmd .= "gzip -cd $reads3 > $fq3\n" if ($reads3=~/gz$/);
				}
			}
			my $prefix="";
			my $out3;
			if ($para eq "" || $para !~ /prefix/) {
				if (@Reads3>1) {
					$prefix="-prefix=$lib\_".($j+1);
					$out3="$lib\_".($j+1)."_uniq.fastq";
				} elsif(@Reads3>0) {
					$prefix="-prefix=$lib";
					$out3="$lib\_uniq.fastq";
				}
			} elsif ($para =~ /prefix\=(\S+)/) {
				if (@Reads3>1) {
					$out3="$1\_".($j+1)."_uniq.fastq";
				} elsif(@Reads3>0) {
					$out3="$1\_uniq.fastq";
				}
			}
			$fltpcr_cmd .= "$filterPCRdup -fastq1=$fq3 $prefix $para\n";
			${$self->{"LIB"}{$lib}{'fq'}}[$j]=qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$out3);
		}
		$fltpcr_cmd .= "cd ..\n";
	}
	$fltpcr_cmd .= "cd ..\n";
	return $fltpcr_cmd;
}

sub runFltAP ($) {
	my $self = shift;
	if (!exists $self->{"software:scanAP"} || !defined $self->{"software:scanAP"}){
		return "";
	}
	my $scanAP=GATE::Error::checkPath($self->{"software:scanAP"});
	my $fastqcut=$self->{"software:fastqcut"} if (exists $self->{"software:fastqcut"});
	my $fastqcut_para = $self->{"setting:fastqcut"} if (exists $self->{"setting:fastqcut"});
	my $align_matrix=$self->{"setting:align.mat"} if (exists $self->{"setting:align.mat"});
	my $para=(exists $self->{"setting:scanAP"})? $self->{"setting:scanAP"} : "";
	my $AP=GATE::Error::checkPath($self->{"database:AP"});
	my $trim_seq=GATE::Error::checkPath($self->{"software:trim_seq"}) if (exists $self->{"software:trim_seq"});
	my $trim_seq_para=$self->{"setting:trim_seq"};
	my $multirun=GATE::Error::checkPath($self->{"software:multithreads-run"}) if (exists $self->{"software:multithreads-run"});
	my $multirun_sh=(-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"})) ? qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/runFltAP.$$.sh) : qq($self->{"-workdir"}/runFltAP.$$.sh);
	my $fltap_cmd = qq(echo `date`; echo "run scanAP"\n);
	if (defined $multirun) {
		open (FLA,">$multirun_sh");
		$fltap_cmd .= qq(export multirun="$multirun"\n);
	}
	$fltap_cmd .=qq(cd $self->{"-workdir"}\n);
	$fltap_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$fltap_cmd .= qq(export AP=$AP\n);
	$fltap_cmd .= qq(mkdir $self->{"setting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}));
	$fltap_cmd .= qq(cd $self->{"setting:qc_outdir"}\n);
	my $fltap_cmd_multi_head=$fltap_cmd;
	my $fltap_cmd_multi="";
	my $multi=0;
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib (@libraries) {
		$fltap_cmd_multi .= $fltap_cmd_multi_head;
		$fltap_cmd_multi .= qq(echo `date` && echo "$lib"\n);
		$fltap_cmd .= qq(echo `date` && echo "$lib"\n);
		$fltap_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
		$fltap_cmd_multi .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
		$fltap_cmd .= "cd $lib\n";
		$fltap_cmd_multi .= "cd $lib && ";
		$fltap_cmd .= qq(cp $align_matrix ./\n) unless (!defined $align_matrix || -f qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/align.mat));
		$fltap_cmd_multi .=  qq(cp $align_matrix ./ && ) unless (!defined $align_matrix || -f qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/align.mat));
		if (defined $multirun) {
			print FLA qq(echo `date` && echo "$lib" && ); 
			print FLA qq([[ -d $lib ]] || mkdir $lib && )unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
			print FLA qq(cd $lib && );
		}
		my %Reads;
		foreach my $i(sort keys %{$self->{"LIB"}{$lib}}) {
			my ($reads1,$reads2,$filter1,$filter2)=("","","","");
			foreach my $reads(@{$self->{"LIB"}{$lib}{$i}}) {
				my $fqname=(split /\//,$reads)[-1];
				my $fq=$reads;
				if (($reads !~ /fastq$/i && $reads !~ /fq$/i) || ($reads =~ /gz$/)) {
					$fq=$1 if ($fqname=~/(\S+)\.gz/);
					$fqname = $fq;
					if (!-f qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$fq))
					{
						$fltap_cmd .= "gzip -cd $reads > $fq\n" if ($reads=~/gz$/);
						if ($reads=~/gz$/){
							$fltap_cmd_multi .= "gzip -cd $reads > $fq && ";
							if (defined $multirun) {
								print FLA "gzip -cd $reads > $fq && ";
							}
						}
					}
				}
				
				my ($stat,$detail)=("$fqname.stat","$fqname.detail");
				if ($fqname=~/(\S+)\.[^\.\s]+$/) {
					$stat="$1.stat";
					$detail="$1.detail";
				}
				if ( ( (!exists $self->{"setting:reuse"}) || (exists $self->{"setting:reuse"} && $self->{"setting:reuse"} =~/N/i) ) && !-f qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib/$detail))
				{
					$fltap_cmd .= "$scanAP -i $fq -a \$AP -s $stat -d $detail $para\n";
					$fltap_cmd_multi .= "$scanAP -i $fq -a \$AP -s $stat -d $detail $para && ";
					if (defined $multirun) {
						print FLA "$scanAP -i $fq -a \$AP -s $stat -d $detail $para && ";
					}
				}
				$fltap_cmd .= "$trim_seq $trim_seq_para -trim_detail $detail $fq\n" if (exists $self->{"software:trim_seq"});
				$fltap_cmd_multi .= "$trim_seq $trim_seq_para -trim_detail $detail $fq  && " if (exists $self->{"software:trim_seq"});
				if (defined $multirun) {
					print FLA "$trim_seq $trim_seq_para -trim_detail $detail $fq  && " if (exists $self->{"software:trim_seq"});
				}
				my $prefix=$1 if ($fqname=~/(\S+)\.[^\.\s]+$/);
				if (defined $fastqcut)
				{
					$fltap_cmd .= "$fastqcut $fqname.trim.out $fastqcut_para -prefix $fqname > $prefix.clean.fastq\n";
					$fltap_cmd_multi .= "$fastqcut $fqname.trim.out $fastqcut_para -prefix $fqname > $prefix.clean.fastq && ";
					if (defined $multirun) {
						print FLA "$fastqcut $fqname.trim.out $fastqcut_para -prefix $fqname > $prefix.clean.fastq && ";
					}
					push @{$Reads{$i}},["$prefix.clean.fastq","$fqname.filter.out $fqname.qcut.out"];
				}
				else
				{
					$fltap_cmd .= "mv $fqname.trim.out $prefix.clean.fastq\n";
					$fltap_cmd_multi .= "mv $fqname.trim.out $prefix.clean.fastq && ";
					if (defined $multirun) {
						print FLA "mv $fqname.trim.out $prefix.clean.fastq && ";
					}
					push @{$Reads{$i}},["$prefix.clean.fastq","$fqname.filter.out"];
				}
				#push @{$self->{"LIB"}{"$sm-clean"}{$i}},$self->{"-workdir"}."/QC/$sm/$prefix.clean.fastq";
			}
		}
		
		my ($p,$s1,$s2)=(0,0,0);
		foreach my $lbmark(sort keys %{$self->{"LIB"}{$lib}}) {
			my ($lb,$i)=($1,$2) if ($lbmark=~/^(\w+)([12])/);
			my $j;
			if (defined $i && $i==2) {
				$j=1;
				next if (exists $Reads{"$lb$j"});
			} elsif (defined $i && $i==1) {
				$j=2;
			}
			for (my $m=0;$m<@{$Reads{$lbmark}};$m++) {
				my ($reads1,$filter1)=@{${$Reads{$lbmark}}[$m]};
				my ($reads2,$filter2)=@{${$Reads{"$lb$j"}}[$m]} if (defined $i && defined $j && exists $self->{"LIB"}{$lib}{"$lb$j"});
				if (defined $reads1 && defined $reads2) {
					if (exists $self->{"software:fltfastq2pe"}){
						my $fltfastq2pe=GATE::Error::checkPath($self->{"software:fltfastq2pe"});
						my ($pair1,$pair2)=($reads1,$reads2);
						$pair1=~s/fastq$/pair.fastq/i;
						$pair2=~s/fastq$/pair.fastq/i;
						my ($single1,$single2)=($reads1,$reads2);
						$single1=~s/fastq$/single.fastq/i;
						$single2=~s/fastq$/single.fastq/i;
						$fltap_cmd .= "$fltfastq2pe -fastq1 $reads1 -fastq2 $reads2 $filter1 $filter2\n";
						$fltap_cmd .= "rm $reads1 $reads2\n" if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1 );
						$fltap_cmd_multi .= "$fltfastq2pe -fastq1 $reads1 -fastq2 $reads2 $filter1 $filter2 && ";
						$fltap_cmd_multi .= "rm $reads1 $reads2 && " if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1 );
						if (defined $multirun) {
							print FLA "$fltfastq2pe -fastq1 $reads1 -fastq2 $reads2 $filter1 $filter2 && ";
							print FLA "rm $reads1 $reads2 && " if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1 );
						}
						${$self->{"LIB"}{$lib}{"$lb$i"}}[$p]=$self->{"-workdir"}."/".$self->{"setting:qc_outdir"}."/$lib/$pair1";
						${$self->{"LIB"}{$lib}{"$lb$j"}}[$p]=$self->{"-workdir"}."/".$self->{"setting:qc_outdir"}."/$lib/$pair2";
						push @{$self->{"LIB"}{$lib}{"fq"}},$self->{"-workdir"}."/".$self->{"setting:qc_outdir"}."/$lib/$single1";
						foreach my $k(keys %{$self->{$lib}{"$lb$i"}{$p}})
						{
							$self->{$lib}{'fq'}{@{$self->{"LIB"}{$lib}{"fq"}}-1}{$k}=$self->{$lib}{"$lb$i"}{$p}{$k};
						}
						push @{$self->{"LIB"}{$lib}{"fq"}},$self->{"-workdir"}."/".$self->{"setting:qc_outdir"}."/$lib/$single2";
						foreach my $k(keys %{$self->{$lib}{"$lb$i"}{$p}})
						{
							$self->{$lib}{'fq'}{@{$self->{"LIB"}{$lib}{"fq"}}-1}{$k}=$self->{$lib}{"$lb$i"}{$p}{$k};
						}
						$p++;
					} else {
						${$self->{"LIB"}{$lib}{"$lb$i"}}[$s1]=$self->{"-workdir"}."/".$self->{"setting:qc_outdir"}."/$lib/$reads1";
						$s1++;
						${$self->{"LIB"}{$lib}{"$lb$j"}}[$s2]=$self->{"-workdir"}."/".$self->{"setting:qc_outdir"}."/$lib/$reads2";
						$s2++;
					}
				} elsif (defined $reads1 && $reads1 ne "") {
					${$self->{"LIB"}{$lib}{$lbmark}}[$s1]=$self->{"-workdir"}."/".$self->{"setting:qc_outdir"}."/$lib/$reads1";
					foreach my $k(keys %{$self->{$lib}{$lbmark}{$p}})
					{
						if (exists $self->{$lib}{$lbmark}{$s1}{$k})
						{
							$self->{$lib}{'fq'}{$s1}{$k}=$self->{$lib}{$lbmark}{$s1}{$k};
						}
						elsif (exists $self->{$lib}{"fq1"}{$s1}{$k})
						{
							$self->{$lib}{'fq'}{$s1}{$k}=$self->{$lib}{"fq1"}{$s1}{$k};
						}
						else
						{
							$self->{$lib}{'fq'}{$s1}{$k}=$self->{$lib}{"fq1"}{0}{$k};
						}
					}
					$s1++;
				} elsif (defined $reads2 && $reads2 ne "") {
					${$self->{"LIB"}{$lib}{$lbmark}}[$s2]=$self->{"-workdir"}."/".$self->{"setting:qc_outdir"}."/$lib/$reads2";
					foreach my $k(keys %{$self->{$lib}{$lbmark}{$p}})
					{
						if (exists $self->{$lib}{$lbmark}{$s2}{$k})
						{
							$self->{$lib}{'fq'}{$s2}{$k}=$self->{$lib}{"$lb$i"}{$s2}{$k};
						}
						elsif (exists $self->{$lib}{"fq1"}{$s2}{$k})
						{
							$self->{$lib}{'fq'}{$s2}{$k}=$self->{$lib}{"fq1"}{$s2}{$k};
						}
						else
						{
							$self->{$lib}{'fq'}{$s2}{$k}=$self->{$lib}{"fq1"}{0}{$k};
						}
					}
					$s2++;
				}
			}
		}
		$fltap_cmd .= "rm *.out\n" if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1 );
		$fltap_cmd_multi .= "rm *.out && " if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1 );
		if (defined $multirun) {
			print FLA "rm *.out && " if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1 );
			print FLA "cd .. && ";
		}
		$fltap_cmd .= qq(rm align.mat\n) unless (!defined $align_matrix || -f "align.mat");
		$fltap_cmd_multi .= qq(rm align.mat && ) unless (!defined $align_matrix || -f "align.mat");
		$fltap_cmd .= "cd ..\n";
		$multi++;
		$fltap_cmd_multi =~ s/\s+$//;
		$fltap_cmd_multi =~ s/\&\&$//;
		if ($multi % $self->{"setting:multithreads"}!=0){
			$fltap_cmd_multi .= " &\n";
		} else {
			$fltap_cmd_multi .= "\n";
		}
		$fltap_cmd_multi .= "cd ..\ncd ..\n";
	}
	if (defined $multirun) {
		print FLA "cd ..\n";
	}
	$fltap_cmd .= "cd ..\n";
	if (exists $self->{"rule:multimode"}) {
		chomp $fltap_cmd_multi;
		$fltap_cmd_multi=~s/\s\&+$//;
		$fltap_cmd_multi.="\n";
		$fltap_cmd_multi.= print_check_process('scanAP','trim_seq','fastqcut','fltfastq2pe');
		if (defined $multirun) {
			my $multirun_cmd = $fltap_cmd_multi_head;
			$multirun_cmd .= qq(${multirun} $multirun_sh -nt $self->{"setting:multithreads"}\n);
			$multirun_cmd .= print_check_process('scanAP','trim_seq','fastqcut','fltfastq2pe');
			return $multirun_cmd;
		} else {
			return $fltap_cmd_multi;
		}
	} else {
		return $fltap_cmd;
	}
}

sub runRSeQC ($) {
	my $self = shift;
	if (!exists $self->{"setting:RSeQCPATH"} || !defined $self->{"setting:RSeQCPATH"}){
		return "";
	}
	my $rseqc_path = $self->{"setting:RSeQCPATH"};
	my $refseqbed = $self->{"database:refseq-bed"};
	
	my $rseqc_cmd = qq(echo `date`; echo "run RSeQC"\n);
	$rseqc_cmd .= qq(cd $self->{"-workdir"}\n);
	$rseqc_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$rseqc_cmd .= qq(export PYTHONPATH=$self->{"setting:PYTHONPATH"}:\$PYTHONPATH\n) if (exists $self->{"setting:PYTHONPATH"});
	$rseqc_cmd .= qq(export PATH=$rseqc_path:\$PATH\n);
	unless (exists $self->{'cmd'}{'aln'} && $self->{'cmd'}{'aln'}==0) {
		$rseqc_cmd .= $self->runBWA("ref");
		$self->{'cmd'}{'aln'} = 0;
	}
	$self->{'cmd'}{'aln'}=0;
	$rseqc_cmd .= qq(export qc_outdir=$self->{"setting:qc_outdir"}\n);
	$rseqc_cmd .= qq([[ -d \${qc_outdir} ]] || mkdir \${qc_outdir}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}));
	$rseqc_cmd .= qq(cd $self->{"setting:qc_outdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib (@libraries) {
		$rseqc_cmd .= qq(echo `date`; echo "$lib"\n);
		$rseqc_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
		$rseqc_cmd .= "cd $lib\n";
		if (exists $self->{$lib}{"ref-bwabam"}) {
			foreach my $bam (@{$self->{$lib}{"ref-bwabam"}}) {
				my $bamName=$1 if ($bam=~/([^\/]+)\.bam$/);
				$rseqc_cmd .= "bam_stat.py  -i $bam > $bam.stat\n";
				$rseqc_cmd .= "clipping_profile.py -i $bam -o $bamName\n";
				if (defined $refseqbed) {
					$rseqc_cmd .= "geneBody_coverage.py -r $refseqbed  -i $bam -o $bamName\n";
					$rseqc_cmd .= "infer_experiment.py -r $refseqbed -i $bam > $bamName\n";
					$rseqc_cmd .= "inner_distance.py -i $bam -r $refseqbed -o $bamName\n";
					$rseqc_cmd .= "junction_annotation.py -i $bam -r $refseqbed -o $bamName\n";
					$rseqc_cmd .= "junction_saturation.py -i $bam -r $refseqbed -o $bamName\n";
				}
				#$rseqc_cmd .= "overlay_bigwig.py\n";
				#$rseqc_cmd .= "normalize_bigwig.py";
				$rseqc_cmd .= "read_distribution.py  -i $bam -r $refseqbed > $bamName\_read_distribution.txt";
				$rseqc_cmd .= "read_duplication.py -i $bam -o $bamName\_duplicate.png\n";
				$rseqc_cmd .= "read_GC.py -i $bam -o $bamName\n";
				$rseqc_cmd .= "read_NVC.py -i $bam -o $bamName\n";
				$rseqc_cmd .= "read_quality.py -i $bam -o $bamName\n";
				#$rseqc_cmd .= "read_hexamer.py";
				$rseqc_cmd .= "RPKM_count.py -d '1++,1--,2+-,2-+'   -i $bam -o $bamName\n";
				$rseqc_cmd .= "RPKM_saturation.py -d '1++,1--,2+-,2-+' -r $refseqbed -i $bam -o $bamName\n";
				$rseqc_cmd .= "split_bam.py -i $bam  -r $refseqbed -o $bamName\n";
			}
		}
	}
	return ($rseqc_cmd);
}

#fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] [-c contaminant file] seqfile1 .. seqfileN
sub runFastQC ($) {
	my $self = shift;
	if (!exists $self->{"software:FastQC"} || !defined $self->{"software:FastQC"}){
		return "";
	}
	my $fastqc=checkPath($self->{"software:FastQC"});
	my $fastqc_cmd = qq(echo `date`; echo "run FastQC"\n);
	$fastqc_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$fastqc_cmd .= qq(export fastqc="$fastqc"\n);
	my $para = (exists $self->{"setting:FastQC"}) ? $self->{"setting:FastQC"} : "" ;
	my $multirun=GATE::Error::checkPath($self->{"software:multithreads-run"}) if (exists $self->{"software:multithreads-run"});
	$fastqc_cmd .= qq(export multirun="$multirun"\n);
	my $multirun_sh=(-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"})) ? qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/runFastQC.$$.sh) : qq($self->{"-workdir"}/runFastQC.$$.sh);
	if (defined $multirun) {
		open (QC,">$multirun_sh");
		$fastqc_cmd .= qq(export multirun="$multirun"\n);
	}
	my $contaminat_file=(exists $self->{"database:contaminat_file"}) ? GATE::Error::checkPath($self->{"database:contaminat_file"}) :
						(exists $self->{"database:AP"}) ? GATE::Error::checkPath($self->{"database:AP"}) : "";
	$fastqc_cmd .= qq(export contamination="$contaminat_file"\n);
	$fastqc_cmd .= qq(cd $self->{"-workdir"}\n);
	$fastqc_cmd .= qq(export qc_outdir=$self->{"setting:qc_outdir"}\n);
	$fastqc_cmd .= qq([[ -d \${qc_outdir} ]] || mkdir \${qc_outdir}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}));
	$fastqc_cmd .= qq(cd $self->{"setting:qc_outdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	my $fastqc_cmd_multi="";
	my $fastqc_cmd_head=$fastqc_cmd;
	foreach my $lib (@libraries) {
		$fastqc_cmd_multi .= $fastqc_cmd_head;
		$fastqc_cmd .= qq(echo `date`; echo "$lib"\n);
		$fastqc_cmd_multi .= qq(echo `date`; echo "$lib"\n);
		$fastqc_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
		$fastqc_cmd_multi .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
		$fastqc_cmd .= "cd $lib\n";
		$fastqc_cmd_multi .= "cd $lib\n";
		if (defined $multirun) {
			print QC qq(echo `date`; echo "$lib" && );
			print QC qq([[ -d $lib ]] || mkdir $lib && ) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
			print QC qq(cd $lib && );
		}
		my %input=getlibInput($self->{"LIB"}{$lib});
		if (exists $input{1} && exists $input{2}){
			for (my $i=0;$i<@{$input{2}};$i++){
				my ($A,$B)=(${$input{1}}[$i],${$input{2}}[$i]);
				my $format1=checkFileFormat($A);
				my $format2=checkFileFormat{$B};
				if ($format1 eq $format2){
					$fastqc_cmd .= qq(\${fastqc} $para -f $format1 );
					$fastqc_cmd .= qq( -c \${contamination} ) if ($contaminat_file ne "");
					$fastqc_cmd .= qq( $A $B\n);
					$fastqc_cmd_multi .= qq(\${fastqc} $para -f $format1 );
					$fastqc_cmd_multi .= qq( -c \${contamination} ) if ($contaminat_file ne "");
					$fastqc_cmd_multi .= qq( $A $B && );
					if (defined $multirun) {
						print QC qq(\${fastqc} $para -f $format1 );
						print QC qq( -c \${contamination} ) if ($contaminat_file ne "");
						print QC qq( $A $B && );
					}
				} else {
					$fastqc_cmd .= qq(\${fastqc} $para -f $format1 );
					$fastqc_cmd .= qq( -c \${contamination} ) if ($contaminat_file ne "");
					$fastqc_cmd .= qq( $A\n);
					$fastqc_cmd .= qq(\${fastqc} $para -f $format2 );
					$fastqc_cmd .= qq( -c \${contamination} ) if ($contaminat_file ne "");
					$fastqc_cmd .= qq( $B\n);
					$fastqc_cmd_multi .= qq(\${fastqc} $para -f $format1 );
					$fastqc_cmd_multi .= qq( -c \${contamination} ) if ($contaminat_file ne "");
					$fastqc_cmd_multi .= qq( $A && );
					$fastqc_cmd_multi .= qq(\${fastqc} $para -f $format2 );
					$fastqc_cmd_multi .= qq( -c \${contamination} ) if ($contaminat_file ne "");
					$fastqc_cmd_multi .= qq( $B && );
					if (defined $multirun) {
						print QC qq(\${fastqc} $para -f $format1 );
						print QC qq( -c \${contamination} ) if ($contaminat_file ne "");
						print QC qq( $A && );
						print QC qq(\${fastqc} $para -f $format2 );
						print QC qq( -c \${contamination} ) if ($contaminat_file ne "");
						print QC qq( $B && );
					}
				}
			}
		}
		if (exists $input{0}) {
			for (my $i=0;$i<@{$input{0}};$i++){
				my $format=checkFileFormat{${$input{0}}[$i]};
				my $C = ${$input{0}}[$i];
				$fastqc_cmd .= qq(\${fastqc} $para -f $format );
				$fastqc_cmd .= qq( -c  $contaminat_file ) if ($contaminat_file ne "");
				$fastqc_cmd .= qq( $C\n);
				$fastqc_cmd_multi .= qq(\${fastqc} $para -f $format );
				$fastqc_cmd_multi .= qq( -c \${contamination} ) if ($contaminat_file ne "");
				$fastqc_cmd_multi .= qq( $C && );
				if (defined $multirun) {
					print QC qq(\${fastqc} $para -f $format );
					print QC qq( -c \${contamination} ) if ($contaminat_file ne "");
					print QC qq( $C && );
				}
			}
		}
		$fastqc_cmd_multi =~ s/[\&\s]+$//;
		$fastqc_cmd_multi =~ s/[\&\s]+$//;
		if ($multi % $self->{"setting:multithreads"}!=0){
			$fastqc_cmd_multi .= " &\n";
		} else {
			$fastqc_cmd_multi .= "\n";
		}
		$fastqc_cmd .= "cd ..\n";
		$fastqc_cmd_multi .= "cd ..\ncd ..\n";
		if (defined $multirun) {
			print QC "cd .. && ";
		}
		$fastqc_cmd .= "cd ..\n";
	}
	if (defined $multirun) {
		print QC qq(cd ..\n);
	}
	if (exists $self->{"rule:multimode"}) {
		chomp $fastqc_cmd_multi;
		$fastqc_cmd_multi=~s/[\s\&]+$//;
		$fastqc_cmd_multi=~s/[\s\&]+$//;
		$fastqc_cmd_multi.="\n";
		$fastqc_cmd_multi.= print_check_process('fastqc');
		if (defined $multirun) {
			my $multirun_cmd = $fastqc_cmd_multi_head;
			$multirun_cmd .= qq(${multirun} $multirun_sh -nt $self->{"setting:multithreads"}\n);
			$multirun_cmd .= print_check_process('fastqc','$multirun_sh');
			return $multirun_cmd;
		} else {
			return $fastqc_cmd_multi;
		}
	} else {
		return $fastqc_cmd;
	}
}

sub runFastX ($) {
	
}

sub runTagDust ($) {
	
}

sub ESTclean ($) {
	
}

sub runRC454 ($) {
	
}

sub runECHO ($) {
	
}

sub runLucy ($) {
	
}

sub runRNASeqQC ($) {
##RNA-SeQC can be run with or without a BWA-based rRNA level estimation mode.
##To run without (less accurate, but faster) use the command:
# java -jar RNASeQC.jar0 -s "TestId|ThousandReads.bam|TestDesc" \
# -t gencode.v7.annotation_goodContig.gtf -r Homo_sapiens_assembly19.fasta \
# -o ./testReport/ -strat gc -gc gencode.v7.gc.txt

##To run the more accurate but slower, BWA-based method :
# java -jar RNASeQC.jar -n 1000 -s "TestId|ThousandReads.bam|TestDesc" \
# -t gencode.v7.annotation_goodContig.gtf -r Homo_sapiens_assembly19.fasta
# -o ./testReport/ -strat gc -gc gencode.v7.gc.txt -BWArRNA human_all_rRNA.fasta
##Note: this assumes BWA is in your PATH. If this is not the case, use the -bwa flag to specify the path to BWA
	my $self = shift;
	if (!exists $self->{"software:RNA-SeQC"} || !defined $self->{"software:RNA-SeQC"}){
		return "";
	}
	my $RNASeQC = $self->{"software:RNA-SeQC"};
	my $reference=GATE::Error::checkPath($self->{"database:ref"});
	my $genecode=GATE::Error::checkPath($self->{"database:gencode-gtf"});
	my $gcfile=GATE::Error::checkPath($self->{"database:gencode-gc"});
	my $para=$self->{'setting:RNA-SeQC'};
	
	my $rnaseqqc_cmd = qq(echo `date`; echo "run RNA-SeQC"\n);
	$rnaseqqc_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$rnaseqqc_cmd .= qq(export qc_outdir=$self->{"setting:qc_outdir"}\n);
	$rnaseqqc_cmd .= qq(cd $self->{"-workdir"}\n);
	$rnaseqqc_cmd .= qq([[ -d \${qc_outdir} ]] || mkdir \${qc_outdir}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}));
	unless (exists $self->{'cmd'}{'aln'} && $self->{'cmd'}{'aln'}==0) {
		$rnaseqqc_cmd .= $self->runBWA("ref");
		$self->{'cmd'}{'aln'}=0;
	}
	$rnaseqqc_cmd .= qq(cd $self->{"setting:qc_outdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		$rnaseqqc_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
		$rnaseqqc_cmd .= "cd $lib\n";
		my $refbam=(@{$self->{$lib}{"ref-bwabam"}}>0)?join ",",@{$self->{$lib}{"ref-bwabam"}}:${$self->{$lib}{"ref-bwabam"}}[0];
		if (exists $self->{"database:rRNA"}) {
			#$rnaseqqc_cmd .= $self->runBWAaln("rRNA");
			$rnaseqqc_cmd .=  qq($RNASeQC -n 1000 -s "$lib|$refbam|$lib" -t $genecode -r $reference -o ./$lib\_RNA-SeQCReport/ -strat gc -gc $gcfile -BWArRNA $self->{"database:rRNA"});
		} else {
			$rnaseqqc_cmd .=  qq($RNASeQC -n 1000 -s "$lib|$refbam|$lib" -t $genecode -r $reference -o ./$lib\_RNA-SeQCReport/ -strat gc -gc $gcfile);
		}
		$rnaseqqc_cmd .= (defined $para) ? "$para\n" : "\n";
	}
	$rnaseqqc_cmd .= "cd ..\n";
	return ($rnaseqqc_cmd);
}

sub runSEECER($) {
#Invalid option: --
#   # This script runs the SEECER pipeline of 4 steps:
#   #
#   # 1. Replace Ns and strip off read IDs (to save memory).
#   # 2. Run JELLYFISH to count kmers.
#   # 3. Correct errors with SEECER.
#   # 4. Clean up and put back original read IDs.
#   
#   run_seecer.sh [options] read1 read2
#
#   read1, read2: are Fasta/Fastaq files.
#	  If only read1 is provided, the reads are considered singles.
#          Otherwise, read1 and read2 are paired-end reads.
#
#   Options:
#      -t <v> : *required* specify a temporary working directory.
#      -k <v> : specify a different K value (default = 17).
#      -j <v> : specify the location of JELLYFISH binary (default = ../jellyfish-1.1.4/bin/jellyfish).
#      -p <v> : specify extra SEECER parameters (default = '').
#      -s <v> : specify the starting step ( default = 1). Values = 1,2,3,4.
#      -h : help message
	my $self = shift;
	if (!exists $self->{"software:SEECER"} || !defined $self->{"software:SEECER"} || !exists $self->{"software:JELLYFISH"}){
		return "";
	}
	my $seecer = $self->{"software:SEECER"};
	my $jellyfish = $self->{"software:JELLYFISH"};
	my $seecerpara = $self->{"setting:SEECER"};
	my $seecer_cmd = qq(echo `date`; echo "run SEECER"\n);;
	$seecer_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$seecer_cmd .= qq(export seecer=$seecer\n);
	$seecer_cmd .= qq(export jellyfish=$jellyfish\n);
	$seecer_cmd .= qq(export seecerpara="$seecerpara"\n);
	$seecer_cmd .= qq(export qc_outdir=$self->{"setting:qc_outdir"}\n);
	$seecer_cmd .= qq(cd $self->{"-workdir"}\n);
	$seecer_cmd .= qq([[ -d \${qc_outdir} ]] || mkdir \${qc_outdir}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}));
	$seecer_cmd .= qq(cd $self->{"setting:qc_outdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		$seecer_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
		$seecer_cmd .= "cd $lib\n";
		my %fq=getlibInput($self->{"LIB"}{$lib});
		if (exists $fq{1} && exists $fq{2}){
			for (my $i=0;$i<@{$fq{2}};$i++) {
				my $fq1=${$fq{1}}[$i];
				my $fq2=${$fq{2}}[$i];
				$seecer_cmd .= qq(\${seecer} \${seecerpara} $fq1 $fq2\n);
			}
		}
		if (exists $fq{0} && @{$fq{0}}>0){
			for (my $i=0;$i<@{$fq{0}};$i++) {
				my $fq=${$fq{0}}[$i];
				$seecer_cmd .= qq(\${seecer} \${seecerpara} $fq\n);
			}
		}
		$seecer_cmd .= "cd ..\n";
	}
	return $seecer_cmd;
}

sub runBEDtools ($) {
	
}

#########################################################
#                                                       #
#                       Alignment                       #
#                                                       #
#########################################################

## MAQ v0.7.1
sub runMAQ ($$) {
#maq fasta2bfa in.ref.fasta out.ref.bfa
#maq sol2sanger in.solexa.fastq out.sanger.fastq
#maq fastq2bfq in.read.fastq out.read.bfq
#maq map [-a  maxins] aln.map in.ref.bfa in.read1.bfq in.read2.bfq 2>out.map.log
#maq mapmerge out.aln.map in.aln1.map in.aln2.map [...]
#maq mapview in.aln.map > out.aln.mapview
#maq assemble [-sp] [-m maxmis] [-Q maxerr] [-r hetrate] [-t coef] [-q minQ] [-M avgmaf] out.cns in.ref.bfa in.aln.map 2> out.cns.log
#maq indelpe in.ref.bfq in.aln.map > out.indelpe
#maq cns2snp in.cns > out.snp
#maq.pl SNPfilter [-d minDep] [-D maxDep] [-Q maxMapQ] [-q minCnsQ] [-w indelWinSize] [-n minNeiQ] [-F in.indelpe] [-f in.indelsoa] [-s minScore] [-m maxAcross] [-a] [-N maxWinSNP] [-W densWinSize] in.cns2snp.snp > out.filtered.snp 

	my $self=shift;
	my $ref=shift;
	if (!exists $self->{"software:maq"} || !defined $self->{"software:maq"} || !-e $self->{"software:maq"}){
		return "";
	}
	my $maq_cmd = qq(echo `date`; echo "run MAQ"\n);
	$maq_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	my $maq =GATE::Error::checkPath($self->{"software:maq"});
	$maq_cmd .= qq(export maq="$maq"\n);
	my $maqpl=GATE::Error::checkPath($self->{"software:maq.pl"});
	$maq_cmd .= qq(export maqpl="$maqpl"\n);
	my $samtools=GATE::Error::checkPath($self->{"software:samtools"});
	my $samtools_path=$1 if ($samtools =~ /(\S+)\/samtools/);
	my $maq2sam;
	$maq2sam=GATE::Error::checkPath("$samtools_path/misc/maq2sam-short") if (-f "$samtools_path/misc/maq2sam-short" && -e "$samtools_path/misc/maq2sam-short");
	$maq2sam=GATE::Error::checkPath($self->{"software:maq2sam"}) if (exists $self->{"software:maq2sam"});
	$maq_cmd .= qq(export samtools="$samtools"\n);
	$maq_cmd .= qq(export maq2sam="$maq2sam"\n) if (defined $maq2sam);
	my $reference=GATE::Error::checkPath($self->{"database:$ref"});
	
	my $mappara="";
	$mappara=$self->{'setting:maqmap'} if (exists $self->{'setting:maqmap'});
	$maq_cmd .= qq(export mappara="$mappara"\n) if (defined $mappara);
	my $assemblepara="";
	$assemblepara=$self->{'setting:maq_assemble'} if (exists $self->{'setting:maq_assemble'});
	$maq_cmd .= qq(export assemblepara="$assemblepara"\n) if (defined $assemblepara);

	my $workdir=GATE::Error::checkPath($self->{"-workdir"});
	$maq_cmd .= qq(export workdir="$workdir"\n);
	my $alndir= GATE::Error::checkPath($self->{"setting:aln_outdir"});
	$maq_cmd .= qq(export alndir="$alndir"\n);
	my $bfa   = "$1.bfa" if ($reference=~/(\S+).fa/);
	$maq_cmd .= "\${maq} fasta2bfa \$REFERENCE \n" unless (GATE::Error::checkIndex('maq',$reference)==1);
	$maq_cmd .= qq(export REFBFA="$bfa"\n);
	$maq_cmd .= qq(cd \${workdir}\n);
	$maq_cmd .= qq(cd \${alndir});
	
	my @libraries=sort keys %{$self->{'LIB'}};
	my $multi=0;
	foreach my $lib(@libraries) {
		$maq_cmd .= qq(echo `date`; echo "$lib"\n);
		$maq_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) if (!-e qq($self->{"-workdir"}/$self->{"setting:aln_outdir"}/$lib));
		$maq_cmd .= "cd $lib\n";
		my %fq=getlibInput($self->{"LIB"}{$lib});
		my @map;
		my $rg="";
		if (exists $fq{1} && exists $fq{2}){
			for (my $j=0;$j<@{$fq{2}};$j++) {
				if (@{$fq{2}}>1) {
					my $k=$j+1;
					my $ID=(exists $self->{$lib}{'fq1'}{$j}{'ID'})?$self->{$lib}{'fq1'}{$j}{'ID'}:"$lib-$k";
					my $SM=(exists $self->{$lib}{'fq1'}{$j}{'SM'})?$self->{$lib}{'fq1'}{$j}{'SM'}:$lib;
					my $LB=(exists $self->{$lib}{'fq1'}{$j}{'LB'})?$self->{$lib}{'fq1'}{$j}{'LB'}:$lib;
					my $PL=(exists $self->{$lib}{'fq1'}{$j}{'PL'})?$self->{$lib}{'fq1'}{$j}{'PL'}:"ILLUMINA";
					my $PI=(exists $self->{$lib}{'fq1'}{$j}{'PI'})?$self->{$lib}{'fq1'}{$j}{'PI'}:500;
					$rg.=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
					foreach my $rb(keys %{$self->{$lib}{'fq1'}{$j}})
					{
						if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
						{
							$rg=~s/\'$//;
							$rg.=qq(\\t$rb:$self->{$lib}{'fq1'}{$j}{$rb}');
						}
					}
					$rg=~s/\'$//;
					$rg.="\\n";
					my $fq1=${$fq{1}}[$j];
					my $fq2=${$fq{2}}[$j];
					my $bfq1="$1.bfq" if ($fq1=~/(\S+).fastq/ || $fq1=~/(\S+)\.fq/);
					my $bfq2="$1.bfq" if ($fq2=~/(\S+).fastq/ || $fq2=~/(\S+)\.fq/);
					$maq_cmd.=qq(\${maq} fastq2bfq $fq1 $bfq1) unless (GATE::Error::checkIndex('maq',$fq1)==1);
					$maq_cmd.=qq( &\n);
					$maq_cmd.=qq(\${maq} fastq2bfq $fq2 $bfq2\n) unless (GATE::Error::checkIndex('maq',$fq1)==1);
					$maq_cmd.=print_check_process('fastq2bfq',"$bfq1","$bfq2");
					if ($PI>600) {
						$maq_cmd.=qq(\${maq} map -A $PI $mappara $lib.pair.$k.aln.map \$REFBFA $bfq1 $bfq2 && );
					} else {
						$maq_cmd.=qq(\${maq} map -a $PI $mappara $lib.pair.$k.aln.map \$REFBFA $bfq1 $bfq2 && );
					}
					#$maq_cmd.=qq(\${maq2sam} $lib.$k.aln.map $lib.$k.aln.sam $rg && ) if (defined $maq2sam);
					#$maq_cmd.=qq(\${maq} assemble \${assemblepara} $lib.$k.cns \$REFBFA $lib.$k.aln.map && );
					$maq_cmd.=($multi % $self->{"setting:multithreads"} != 0) ? " &\n" : "\n";
					$multi++;
					push @map, "$lib.pair.$k.aln.map";
					
				}else{
					my $ID=(exists $self->{$lib}{'fq1'}{$j}{'ID'})?$self->{$lib}{'fq1'}{$j}{'ID'}:$lib;
					my $SM=(exists $self->{$lib}{'fq1'}{$j}{'SM'})?$self->{$lib}{'fq1'}{$j}{'SM'}:$lib;
					my $LB=(exists $self->{$lib}{'fq1'}{$j}{'LB'})?$self->{$lib}{'fq1'}{$j}{'LB'}:$lib;
					my $PL=(exists $self->{$lib}{'fq1'}{$j}{'PL'})?$self->{$lib}{'fq1'}{$j}{'PL'}:"ILLUMINA";
					my $PI=(exists $self->{$lib}{'fq1'}{$j}{'PI'})?$self->{$lib}{'fq1'}{$j}{'PI'}:500;
					$rg.=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
					foreach my $rb(keys %{$self->{$lib}{'fq'}{$j}})
					{
						if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
						{
							$rg=~s/\'$//;
							$rg.=qq(\\t$rb:$self->{$lib}{'fq'}{$j}{$rb}');
						}
					}
					$rg=~s/\'$//;
					$rg.="\\n";
					my $fq1=${$fq{1}}[$j];
					my $fq2=${$fq{2}}[$j];
					my $bfq1="$1.bfq" if ($fq1=~/(\S+).fastq/ || $fq1=~/(\S+)\.fq/);
					my $bfq2="$1.bfq" if ($fq2=~/(\S+).fastq/ || $fq2=~/(\S+)\.fq/);
					$maq_cmd.=qq(\${maq} fastq2bfq $fq1 $bfq1) unless (GATE::Error::checkIndex('maq',$fq1)==1);
					$maq_cmd.=qq( &\n);
					$maq_cmd.=qq(\${maq} fastq2bfq $fq2 $bfq2\n) unless (GATE::Error::checkIndex('maq',$fq1)==1);
					$maq_cmd.=print_check_process('fastq2bfq',"$bfq1","$bfq2");
					if ($PI>600) {
						$maq_cmd.=qq(\${maq} map -A $PI $mappara $lib.pair.aln.map \$REFBFA $bfq1 $bfq2 && );
					} else {
						$maq_cmd.=qq(\${maq} map -a $PI $mappara $lib.pair.aln.map \$REFBFA $bfq1 $bfq2 && );
					}
					#$maq_cmd.=qq(\${maq2sam} $lib.aln.map $lib.aln.sam $rg && ) if (defined $maq2sam);
					#$maq_cmd.=qq(\${maq} assemble \${assemblepara} $lib.cns \$REFBFA $lib.aln.map && );
					$maq_cmd.=($multi % $self->{"setting:multithreads"} != 0) ? " &\n" : "\n";
					$multi++;
					push @map, "$lib.pair.aln.map";
				}
			}
		}
		if (exists $fq{0}) {
			for (my $j=0;$j<@{$fq{0}};$j++) {
				if (@{$fq{0}}>1) {
					my $k=$j+1;
					my $ID=(exists $self->{$lib}{'fq'}{$j}{'ID'})?$self->{$lib}{'fq'}{$j}{'ID'}:"$lib-$k";
					my $SM=(exists $self->{$lib}{'fq'}{$j}{'SM'})?$self->{$lib}{'fq'}{$j}{'SM'}:$lib;
					my $LB=(exists $self->{$lib}{'fq'}{$j}{'LB'})?$self->{$lib}{'fq'}{$j}{'LB'}:$lib;
					my $PL=(exists $self->{$lib}{'fq'}{$j}{'PL'})?$self->{$lib}{'fq'}{$j}{'PL'}:"ILLUMINA";
					my $PI=(exists $self->{$lib}{'fq'}{$j}{'PI'})?$self->{$lib}{'fq'}{$j}{'PI'}:500;
					$rg.=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
					foreach my $rb(keys %{$self->{$lib}{'fq'}{$j}})
					{
						if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
						{
							$rg=~s/\'$//;
							$rg.=qq(\\t$rb:$self->{$lib}{'fq'}{$j}{$rb}');
						}
					}
					$rg=~s/\'$//;
					$rg.="\\n";
					my $fq=${$fq{0}}[$j];
					my $bfq="$1.bfq" if ($fq=~/(\S+).fastq/ || $fq=~/(\S+)\.fq/);
					$maq_cmd.=qq(\${maq} fastq2bfq $fq $bfq) unless (GATE::Error::checkIndex('maq',$fq)==1);
					$maq_cmd.=print_check_process('fastq2bfq','$bfq');
					$maq_cmd.=qq(\${maq} map $mappara $lib.single.$k.aln.map \$REFBFA $bfq && );
					$maq_cmd.=($multi % $self->{"setting:multithreads"} != 0) ? " &\n" : "\n";
					$multi++;
					push @map, "$lib.single.$k.aln.map";
					
				}else{
					my $ID=(exists $self->{$lib}{'fq'}{$j}{'ID'})?$self->{$lib}{'fq'}{$j}{'ID'}:"$lib";
					my $SM=(exists $self->{$lib}{'fq'}{$j}{'SM'})?$self->{$lib}{'fq'}{$j}{'SM'}:$lib;
					my $LB=(exists $self->{$lib}{'fq'}{$j}{'LB'})?$self->{$lib}{'fq'}{$j}{'LB'}:$lib;
					my $PL=(exists $self->{$lib}{'fq'}{$j}{'PL'})?$self->{$lib}{'fq'}{$j}{'PL'}:"ILLUMINA";
					my $PI=(exists $self->{$lib}{'fq'}{$j}{'PI'})?$self->{$lib}{'fq'}{$j}{'PI'}:500;
					$rg.=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
					foreach my $rb(keys %{$self->{$lib}{'fq'}{$j}})
					{
						if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
						{
							$rg=~s/\'$//;
							$rg.=qq(\\t$rb:$self->{$lib}{'fq'}{$j}{$rb}');
						}
					}
					foreach my $rb(keys %{$self->{$lib}{'fq'}{$j}})
					{
						if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
						{
							$rg=~s/\'$//;
							$rg.=qq(\\t$rb:$self->{$lib}{'fq'}{$j}{$rb}');
						}
					}
					$rg=~s/\'$//;
					$rg.="\\n";
					my $fq=${$fq{0}}[$j];
					my $bfq="$1.bfq" if ($fq=~/(\S+).fastq/ || $fq=~/(\S+)\.fq/);
					$maq_cmd.=qq(\${maq} fastq2bfq $fq $bfq) unless (GATE::Error::checkIndex('maq',$fq)==1);
					$maq_cmd.=print_check_process('fastq2bfq','$bfq');
					$maq_cmd.=qq(\${maq} map $mappara $lib.single.aln.map \$REFBFA $bfq && );
					$maq_cmd.=($multi % $self->{"setting:multithreads"} != 0) ? " &\n" : "\n";
					$multi++;
					push @map, "$lib.single.aln.map";
				}
			}
		}
		$maq_cmd.=print_check_process('maq','map',"$lib");
		my $libmap;
		if (@map>1) {
			my $allmap=join " ",@map;
			$maq_cmd.=qq(\${maq} mapmerge $lib.merge.aln.map $allmap && );
			$libmap = "$lib.merge.aln.map";
		} elsif (@map>0) {
			$libmap = $map[0];			
		}
		if (defined $libmap) {
			$maq_cmd.=qq(\${maq} assemble \${assemblepara} $lib.cns \$REFBFA $libmap && );
			$maq_cmd.=qq(\${maq}cns2snp $lib.cns > $lib.snp && );
			$maq_cmd.=qq(\${maq} indelpe \$REFBFA $libmap > out.indelpe && );
			$maq_cmd.=qq(\${maqpl} SNPfilter $lib.snp > $lib.filtered.snp && );
			$maq_cmd.=qq(\${maq2sam} $libmap $rg) if (defined $maq2sam); 
			$maq_cmd.=($multi % $self->{"setting:multithreads"} != 0) ? " &\n" : "\n";
			$multi++;
			push @{$self->{$lib}{"$ref-maqmap"}},$libmap;
		}
		$maq_cmd.="cd ..\n";
	}
	$maq_cmd.="cd ..\n";
	return $maq_cmd;
}

## BWA 0.7.4-r385
## SAMtools 0.1.19-44428cd
sub runBWA($$) {
	my $self=shift;
	my $ref=shift;
	if (!exists $self->{"software:bwa"} || !defined $self->{"software:bwa"} || !-e $self->{"software:bwa"}){
		return "";
	}
	$ref ||= 'ref';
	my $bwa_cmd = qq(echo `date`; echo "run BWA"\n);
	$bwa_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	my $bwa=GATE::Error::checkPath($self->{"software:bwa"});
	$bwa_cmd .= qq(export bwa="$bwa"\n);
	my $samtools=GATE::Error::checkPath($self->{"software:samtools"});
	$bwa_cmd .= qq(export samtools="$samtools"\n);
	my $reference=GATE::Error::checkPath($self->{"database:$ref"});
	$bwa_cmd .= qq(export REFERENCE="$reference"\n);
	my $alnpara=$self->{'setting:bwaaln'} if (exists $self->{'setting:bwaaln'});
	my $mempara=$self->{'setting:bwamem'} if (exists $self->{'setting:bwamem'});
	if ($alnpara!~/\-t\s+\d+/ && exists $self->{'setting:multithreads'})
	{
		$alnpara.=" -t ".$self->{'setting:multithreads'};
	}
	$bwa_cmd .= qq(export alnpara="$alnpara"\n) if (defined $alnpara);
	$bwa_cmd .= qq(export mempara="$mempara"\n) if (defined $mempara);
	my $sampepara="";
	$sampepara=$self->{'setting:sampe'} if (exists $self->{'setting:sampe'});
	my $samsepara="";
	$samsepara=$self->{'setting:samse'} if (exists $self->{'setting:samse'});
	$bwa_cmd .= qq(export heap="$self->{'setting:heap'}"\n);
	my ($picard,$FixMateInformation)=("","");
	if (exists $self->{"software:picard"})
	{
		$picard=GATE::Error::checkPath($self->{"software:picard"});
		my $picardpath=$1 if ($self->{"software:picard"}=~/(.*)\/[^\/\s]+$/);
		$bwa_cmd .= qq(export picard="$picard"\n);
		$picard = correctJavaCmd($picard,"\${heap}");
		$FixMateInformation=(exists $self->{"software:FixMateInformation"})?$self->{"software:FixMateInformation"}:qq($picardpath/FixMateInformation.jar);
	}
	
	my $workdir=GATE::Error::checkPath($self->{"-workdir"});
	$bwa_cmd .= qq(export workdir="$workdir"\n);
	my $alndir=GATE::Error::checkPath($self->{"setting:aln_outdir"});
	$bwa_cmd .= qq(export alndir="$alndir"\n);
	$bwa_cmd .= qq(cd \${workdir}\n);
	$bwa_cmd .= qq([[ -d \${alndir} ]] || mkdir \${alndir}\n) if (!-d qq($self->{"-workdir"}/$alndir));;
	$bwa_cmd .= qq(cd \${alndir}\n);
	$bwa_cmd .= "\${bwa} index -a bwtsw \$REFERENCE\n" unless (GATE::Error::checkIndex('bwa',$reference)==1);
	my $db=$1 if ($reference =~ /([^\/\.]+)\.fa/);
	my @libraries=sort keys %{$self->{'LIB'}};
	my $skip = 0;
	foreach my $lib(@libraries) {
		$bwa_cmd .= qq(echo `date`; echo "$lib"\n);
		$bwa_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) if (!-e qq($self->{"-workdir"}/$self->{"setting:aln_outdir"}/$lib));
		$bwa_cmd .= "cd $lib\n";
		if (exists $self->{$lib}{"$ref-bwabam"}) {
			$skip++;
			$bwa_cmd .= "cd ../\n";
			next;
		}
		my %fq=getlibInput($self->{"LIB"}{$lib});
		if (exists $fq{1} && exists $fq{2}){
			for (my $j=0;$j<@{$fq{2}};$j++) {
				my $bam="";
				if (defined $mempara) {
					my $rg="";
					if (@{$fq{2}}>1) {
						my $k=$j+1;
						my $ID=(exists $self->{$lib}{'fq1'}{$j}{'ID'})?$self->{$lib}{'fq1'}{$j}{'ID'}:"$lib-$k";
						my $SM=(exists $self->{$lib}{'fq1'}{$j}{'SM'})?$self->{$lib}{'fq1'}{$j}{'SM'}:$lib;
						my $LB=(exists $self->{$lib}{'fq1'}{$j}{'LB'})?$self->{$lib}{'fq1'}{$j}{'LB'}:$lib;
						my $PL=(exists $self->{$lib}{'fq1'}{$j}{'PL'})?$self->{$lib}{'fq1'}{$j}{'PL'}:"ILLUMINA";
						$rg=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
						foreach my $rb(keys %{$self->{$lib}{'fq1'}{$j}})
						{
							if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
							{
								$rg=~s/\'$//;
								$rg.=qq(\\t$rb:$self->{$lib}{'fq1'}{$j}{$rb}');
							}
						}
						$bwa_cmd .= "\${bwa} mem \${mempara} -R $rg \$REFERENCE ${$fq{1}}[$j] ${$fq{2}}[$j] | \${samtools} view -Sbh - -o $lib.pair.$k.bam\n";
						#$bwa_cmd .= "\${samtools} rmdup $lib.pair.$k.bam $lib.pair.$k.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.pair.$k.bam $lib.pair.$k.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.pair.$k.sort.bam\n";
						$bam="$lib.pair.$k.sort.bam";
					} else {
						my $ID=(exists $self->{$lib}{'fq1'}{$j}{'ID'})?$self->{$lib}{'fq1'}{$j}{'ID'}:$lib;
						my $SM=(exists $self->{$lib}{'fq1'}{$j}{'SM'})?$self->{$lib}{'fq1'}{$j}{'SM'}:$lib;
						my $LB=(exists $self->{$lib}{'fq1'}{$j}{'LB'})?$self->{$lib}{'fq1'}{$j}{'LB'}:$lib;
						my $PL=(exists $self->{$lib}{'fq1'}{$j}{'PL'})?$self->{$lib}{'fq1'}{$j}{'PL'}:"ILLUMINA";
						$rg=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
						foreach my $rb(keys %{$self->{$lib}{'fq1'}{$j}})
						{
							if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
							{
								$rg=~s/\'$//;
								$rg.=qq(\\t$rb:$self->{$lib}{'fq1'}{$j}{$rb}');
							}
						}
						$bwa_cmd .= "\${bwa} mem \${mempara} -R $rg \$REFERENCE ${$fq{1}}[$j] ${$fq{2}}[$j] | \${samtools} view -Sbh - -o $lib.pair.bam\n";
						#$bwa_cmd .= "\${samtools} rmdup $lib.pair.bam $lib.pair.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.pair.bam $lib.pair.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.pair.sort.bam\n";
						$bam="$lib.pair.sort.bam";
					}
				} elsif (defined $alnpara) {
					my $sai1=(split /\//,${$fq{1}}[$j])[-1];
					my $sai2=(split /\//,${$fq{2}}[$j])[-1];
					$sai1=~s/\.gz//i;
					$sai1=~s/\.fq//i;
					$sai1=~s/\.fastq//i;
					$sai2=~s/\.gz//i;
					$sai2=~s/\.fq//i;
					$sai2=~s/\.fastq//i;
					$sai1.=".sai";
					$sai2.=".sai";
					$bwa_cmd .= ("\${bwa} aln \${alnpara} \$REFERENCE ${$fq{1}}[$j] > $sai1\n");
					$bwa_cmd .= ("\${bwa} aln \${alnpara} \$REFERENCE ${$fq{2}}[$j] > $sai2\n");
					if (@{$fq{2}}>1) {
						my $k=$j+1;
						my $ID=(exists $self->{$lib}{'fq1'}{$j}{'ID'})?$self->{$lib}{'fq1'}{$j}{'ID'}:"$lib-$k";
						my $SM=(exists $self->{$lib}{'fq1'}{$j}{'SM'})?$self->{$lib}{'fq1'}{$j}{'SM'}:$lib;
						my $LB=(exists $self->{$lib}{'fq1'}{$j}{'LB'})?$self->{$lib}{'fq1'}{$j}{'LB'}:$lib;
						my $PI=(exists $self->{$lib}{'fq1'}{$j}{'PI'})?$self->{$lib}{'fq1'}{$j}{'PI'}:500;
						my $PL=(exists $self->{$lib}{'fq1'}{$j}{'PL'})?$self->{$lib}{'fq1'}{$j}{'PL'}:"ILLUMINA";
						my $rg=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
						foreach my $rb(keys %{$self->{$lib}{'fq1'}{$j}})
						{
							if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
							{
								$rg=~s/\'$//;
								$rg.=qq(\\t$rb:$self->{$lib}{'fq1'}{$j}{$rb}');
							}
						}
						$bwa_cmd .= qq(\${bwa} sampe $sampepara -a $PI -r $rg \$REFERENCE $sai1 $sai2 ${$fq{1}}[$j] ${$fq{2}}[$j] | \${samtools} view -Sbh - -o $lib.pair.$k.bam\n);
						#$bwa_cmd .= "\${samtools} rmdup $lib.pair.$k.bam $lib.pair.$k.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.pair.$k.bam $lib.pair.$k.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.pair.$k.sort.bam\n";
						$bam="$lib.pair.$k.sort.bam";
					} else {
						my $ID=(exists $self->{$lib}{'fq1'}{$j}{'ID'})?$self->{$lib}{'fq1'}{$j}{'ID'}:$lib;
						my $SM=(exists $self->{$lib}{'fq1'}{$j}{'SM'})?$self->{$lib}{'fq1'}{$j}{'SM'}:$lib;
						my $LB=(exists $self->{$lib}{'fq1'}{$j}{'LB'})?$self->{$lib}{'fq1'}{$j}{'LB'}:$lib;
						my $PI=(exists $self->{$lib}{'fq1'}{$j}{'PI'})?$self->{$lib}{'fq1'}{$j}{'PI'}:500;
						my $PL=(exists $self->{$lib}{'fq1'}{$j}{'PL'})?$self->{$lib}{'fq1'}{$j}{'PL'}:"ILLUMINA";
						my $rg=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
						foreach my $rb(keys %{$self->{$lib}{'fq1'}{$j}})
						{
							if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
							{
								$rg=~s/\'$//;
								$rg.=qq(\\t$rb:$self->{$lib}{'fq1'}{$j}{$rb}');
							}
						}
						$bwa_cmd .= "\${bwa} sampe $sampepara -a $PI -r $rg \$REFERENCE $sai1 $sai2 ${$fq{1}}[$j] ${$fq{2}}[$j] | \${samtools} view -Sbh - -o $lib.pair.bam\n";
						#$bwa_cmd .= "\${samtools} rmdup $lib.pair.bam $lib.pair.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.pair.bam $lib.pair.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.pair.sort.bam\n";
						$bam="$lib.pair.sort.bam";
					}
				}
## Fix mate
#java -Xmx${heap}m -Djava.io.tmpdir\=${tmp_folder}_fixmate \
# -jar ${picard}/FixMateInformation.jar \
# INPUT\=$PWDS/${subjectID}.srt.bam \
# OUTPUT\=$PWDS/${subjectID}.fxmt.bam \
# SO\=coordinate \
# CREATE_INDEX\=true  \
# VALIDATION_STRINGENCY\=SILENT
				if (exists $self->{"software:picard"}) {
					$bwa_cmd .= qq(export FixMateInformation="$FixMateInformation"\n);
					#$bwa_cmd .= qq(mkdir -p ./tmp_fixmate\n);
					#$bwa_cmd .= qq(export tmp_fixmate="./tmp_fixmate"\n);
					$FixMateInformation=correctJavaCmd($FixMateInformation,"\${heap}","./tmp_fixmate");
					my $FixMateInformationPara=(exists $self->{"setting:FixMateInformation"})? $self->{"setting:FixMateInformation"}:"SO=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT";
					my $fxmtbam=$bam;
					$fxmtbam=~s/bam$/fxmat\.bam/;
					$bwa_cmd .= qq($FixMateInformation INPUT\=$bam OUTPUT=$fxmtbam $FixMateInformationPara\n);
					if (exists $self->{"software:bamtools"})
					{
## Generate stat
#
#${bamtools} stats \
#  -insert \
#  -in $PWDS/${subjectID}.fxmt.bam \
#  > $PWDS/${subjectID}.fxmt.stats
						my $bamtools=GATE::Error::checkPath($self->{"software:bamtools"});
						$bwa_cmd .= qq(export bamtools\="$bamtools"\n);
						my $stats="$1.stats" if ($fxmtbam=~/(\S+)\.bam/);
						$bwa_cmd .= qq(\${bamtools} stats -insert -in $fxmtbam > $stats\n);
					}
					$bwa_cmd .= qq(rm -rf tmp_fixmate\n);
					$bam=$fxmtbam;
				}
				push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$bam";
				$bwa_cmd .= "rm *.sai *.pair.bam\n" if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1 );
			}
		}
		if (exists $fq{0} && @{$fq{0}}>0){
			for (my $j=0;$j<@{$fq{0}};$j++) {
				if (defined $mempara) {
					my $rg="";
					if (@{$fq{2}}>1) {
						my $k=$j+1;
						my $ID=(exists $self->{$lib}{'fq'}{$j}{'ID'})?$self->{$lib}{'fq'}{$j}{'ID'}:"$lib-$k";
						my $SM=(exists $self->{$lib}{'fq'}{$j}{'SM'})?$self->{$lib}{'fq'}{$j}{'SM'}:$lib;
						my $LB=(exists $self->{$lib}{'fq'}{$j}{'LB'})?$self->{$lib}{'fq'}{$j}{'LB'}:$lib;
						my $PL=(exists $self->{$lib}{'fq'}{$j}{'PL'})?$self->{$lib}{'fq'}{$j}{'PL'}:"ILLUMINA";
						$rg=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
						foreach my $rb(keys %{$self->{$lib}{'fq'}{$j}})
						{
							if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
							{
								$rg=~s/\'$//;
								$rg.=qq(\\t$rb:$self->{$lib}{'fq'}{$j}{$rb}');
							}
						}
						$bwa_cmd .= "\${bwa} mem \${mempara} -R $rg \$REFERENCE ${$fq{0}}[$j] | \${samtools} view -Sbh - -o $lib.single.$k.bam\n";
						#$bwa_cmd .= "\${samtools} rmdup $lib.single.$k.bam $lib.single.$k.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.single.$k.bam $lib.single.$k.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.single.$k.sort.bam\n";
						push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.single.$k.sort.bam";
					} else {
						my $ID=(exists $self->{$lib}{'fq'}{$j}{'ID'})?$self->{$lib}{'fq'}{$j}{'ID'}:$lib;
						my $SM=(exists $self->{$lib}{'fq'}{$j}{'SM'})?$self->{$lib}{'fq'}{$j}{'SM'}:$lib;
						my $LB=(exists $self->{$lib}{'fq'}{$j}{'LB'})?$self->{$lib}{'fq'}{$j}{'LB'}:$lib;
						my $PL=(exists $self->{$lib}{'fq'}{$j}{'PL'})?$self->{$lib}{'fq'}{$j}{'PL'}:"ILLUMINA";
						$rg=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
						foreach my $rb(keys %{$self->{$lib}{'fq'}{$j}})
						{
							if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
							{
								$rg=~s/\'$//;
								$rg.=qq(\\t$rb:$self->{$lib}{'fq'}{$j}{$rb}');
							}
						}
						$bwa_cmd .= "\${bwa} mem \${mempara} -R $rg \$REFERENCE ${$fq{0}}[$j] | \${samtools} view -Sbh - -o $lib.single.bam\n";
						#$bwa_cmd .= "\${samtools} rmdup $lib.single.bam $lib.single.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.single.bam $lib.single.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.single.sort.bam\n";
						push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.single.sort.bam";
					}
					$bwa_cmd .= "rm *.single.bam *.single.?.bam\n" if (exists $self->{"rule:Clean"});
				}elsif (defined $alnpara) {
					my $sai1=(split /\//,${$fq{0}}[$j])[-1];
					$sai1=~s/\.gz//i;
					$sai1=~s/\.fq//i;
					$sai1=~s/\.fastq//i;
					$sai1.=".sai";
					$bwa_cmd .= ("\$bwa aln \${alnpara} \$REFERENCE ${$fq{0}}[$j] > $sai1\n");
					if (@{$fq{0}}>1) {
						my $k=$j+1;
						my $ID=(exists $self->{$lib}{'fq'}{$j}{'ID'})?$self->{$lib}{'fq'}{$j}{'ID'}:"$lib-$k";
						my $SM=(exists $self->{$lib}{'fq'}{$j}{'SM'})?$self->{$lib}{'fq'}{$j}{'SM'}:$lib;
						my $LB=(exists $self->{$lib}{'fq'}{$j}{'LB'})?$self->{$lib}{'fq'}{$j}{'LB'}:$lib;
						my $PL=(exists $self->{$lib}{'fq'}{$j}{'PL'})?$self->{$lib}{'fq'}{$j}{'PL'}:"ILLUMINA";
						my $rg=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
						foreach my $rb(keys %{$self->{$lib}{'fq'}{$j}})
						{
							if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
							{
								$rg=~s/\'$//;
								$rg.=qq(\\t$rb:$self->{$lib}{'fq'}{$j}{$rb}');
							}
						}
						$bwa_cmd .= "\$bwa samse $samsepara -r $rg \$REFERENCE $sai1 ${$fq{0}}[$j] | \${samtools} view -Sbh - -o $lib.single.$k.bam\n";
						#$bwa_cmd .= "\${samtools} rmdup $lib.single.$k.bam $lib.single.$k.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -\@ 8 -m 3G  $lib.single.$k.bam $lib.single.$k.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.single.$k.sort.bam\n";
						push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.single.$k.sort.bam";
					} else {
						my $ID=(exists $self->{$lib}{'fq'}{$j}{'ID'})?$self->{$lib}{'fq'}{$j}{'ID'}:$lib;
						my $SM=(exists $self->{$lib}{'fq'}{$j}{'SM'})?$self->{$lib}{'fq'}{$j}{'SM'}:$lib;
						my $PL=(exists $self->{$lib}{'fq'}{$j}{'PL'})?$self->{$lib}{'fq'}{$j}{'PL'}:"ILLUMINA";
						my $LB=(exists $self->{$lib}{'fq'}{$j}{'LB'})?$self->{$lib}{'fq'}{$j}{'LB'}:$lib;
						my $rg=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
						foreach my $rb(keys %{$self->{$lib}{'fq'}{$j}})
						{
							if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
							{
								$rg=~s/\'$//;
								$rg.=qq(\\t$rb:$self->{$lib}{'fq'}{$j}{$rb}');
							}
						}
						$bwa_cmd .= "\$bwa samse $samsepara -r $rg \$REFERENCE $sai1 ${$fq{0}}[$j] | \${samtools} view -Sbh - -o $lib.single.bam\n";
						#$bwa_cmd .= "\${samtools} rmdup $lib.single.bam $lib.single.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.single.bam $lib.single.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.single.sort.bam\n";
						push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.single.sort.bam";
					}
					$bwa_cmd .= "rm *.sai *.single.bam *.single.?.bam\n" if (exists $self->{"rule:Clean"});
				}
			}
		}
		if (exists $self->{'rule:bwamerge'} && exists $self->{"software:picard"}) {
			my $MergeSamFiles="$1/MergeSamFiles.jar" if ($self->{"software:picard"}=~/(.*)\/[^\/\s]+$/);
			$MergeSamFiles=qq(java -Xmx\${heap} -Djava.io.tmpdir=./tmp_merge -jar $MergeSamFiles) if ($MergeSamFiles!~/^java/ && $MergeSamFiles !~ /\-jar/);
			my $MergeSamFilesPara=$self->{"setting:MergeSamFiles"};
			my $merge_bam=join " INPUT=",@{$self->{$lib}{"ref-bwabam"}};
			$bwa_cmd .= "$MergeSamFiles INPUT=$merge_bam $MergeSamFilesPara OUTPUT=$lib.merge.bam\n";
			my $bam="$lib.merge.bam";
			$bwa_cmd .= "\${samtools} index $lib.merge.bam\n";
			if (!exists $self->{rmdup} || GATE::Error::boolean($self->{rmdup})==1) {
				$bwa_cmd .= "\${samtools} rmdup $lib.merge.bam - | \${samtools} rmdup -S - - | \${samtools} sort - $lib.merge.rmdup.sort\n";
				$bam="$lib.merge.rmdup.sort.bam";
				$bwa_cmd .= "\${samtools} index $bam\n";
			}
			$bwa_cmd .= "rm -rf ./tmp_merge\n" if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1);
			@{$self->{$lib}{"$ref-bwabam"}}=();
			$bam=$self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$bam";
			push @{$self->{$lib}{"$ref-bwabam"}},$bam;
			if (exists $self->{"overlap"} && exists $self->{"sam2bed"} && exists $self->{"msort"}) {
				$bwa_cmd .= $self->stat_mappedreads("bam",$bam,"lib",$lib);
			}
			if (exists $self->{"sam2bed"}){
				my $bed="$1.bed" if ($bam=~/([^\/\s]+)\.bam$/);
				$bwa_cmd .= qq( [[ -f $bed ]] || \${samtools} view -F4 $bam | $self->{"sam2bed"} -n > $bed\n);
				push @{$self->{'LIB'}{$lib}{'INPUT'}{'BED'}},$self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$bed";
			}
		}
		if (exists $self->{'rule:UnmappedRealign'})
		{
			if (exists $self->{$lib}{"$ref-bwabam"} && exists $self->{$lib}{"software:sam2reads"} && exists $self->{$lib}{"software:fltfastq2pe"})
			{
				my $sam2reads=GATE::Error::checkPath($self->{$lib}{"software:sam2reads"});
				$bwa_cmd .= qq(export sam2reads="$sam2reads"\n);
				my $fltfastq2pe=GATE::Error::checkPath($self->{$lib}{"software:fltfastq2pe"});
				$bwa_cmd .= qq(export fltfastq2pe="$fltfastq2pe"\n);
				@{$self->{"LIB"}{$lib}{'fq1'}}=();
				@{$self->{"LIB"}{$lib}{'fq2'}}=();
				@{$self->{"LIB"}{$lib}{'fq'}}=();
				my $i=1;
				foreach my $bwabam(@{$self->{$lib}{"$ref-bwabam"}})
				{
					if (@{$self->{$lib}{"$ref-bwabam"}}>1)
					{
						$bwa_cmd .= qq(\${samtools} -f 4 $bwabam | \${sam2reads} -R1 $lib.$i-unmpped.R1.fastq -R2 $lib.$i-unmpped.R2.fastq -R3 $lib.$i-unmpped.R3.fastq -f fq\n);
						$bwa_cmd .= qq(\${fltfastq2pe} -fastq1 $lib.$i-unmpped.R1.fastq -fastq2 $lib.$i-unmpped.R2.fastq\n);
						push @{$self->{"LIB"}{$lib}{'fq1'}},qq($self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.$i-unmpped.R1.pair.fastq);
						push @{$self->{"LIB"}{$lib}{'fq2'}},qq($self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.$i-unmpped.R2.pair.fastq);
						push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.$i-unmpped.R1.single.fastq);
						push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.$i-unmpped.R2.single.fastq);
						push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.$i-unmpped.R3.fastq);
					}
					else
					{
						$bwa_cmd .= qq(\${samtools} -f 4 $bwabam | \${sam2reads} -R1 $lib.unmpped.R1.fastq -R2 $lib.unmpped.R2.fastq -R3 $lib.unmpped.R3.fastq -f fq\n);
						$bwa_cmd .= qq(\${fltfastq2pe} -fastq1 $lib.unmpped.R1.fastq -fastq2 $lib.unmpped.R2.fastq\n);
						push @{$self->{"LIB"}{$lib}{'fq1'}},qq($self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.unmpped.R1.pair.fastq);
						push @{$self->{"LIB"}{$lib}{'fq2'}},qq($self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.unmpped.R2.pair.fastq);
						push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.unmpped.R1.single.fastq);
						push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.unmpped.R2.single.fastq);
						push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.unmpped.R3.fastq);
					}
					$i++;
				}
			}
		}
		$bwa_cmd.="cd ..\n";
	}
	$bwa_cmd.="cd ..\n";
	if ($skip<@libraries) {
		return $bwa_cmd unless (exists $self->{'rule:skip'}{'aln'} && GATE::Error::boolean($self->{'rule:skip'}{'aln'})==1);
	} else {
		return "";
	}
}

sub runBowtie($$) {
	my ($self,$ref) = @_;
	if (!exists $self->{"software:bowtie"} || !-e $self->{"software:bowtie"} || !defined $self->{"software:bowtie"} || defined $self->{"software:tophat"}) {
		return "";
	}
	$ref ||= 'ref';
	my $bowtie_cmd = qq(echo `date`; echo "run TopHat"\n);
	$bowtie_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	my $bowtie=GATE::Error::checkPath($self->{'software:bowtie'});
	$bowtie_cmd.="export bowtie=$bowtie\n";
	my $samtools = GATE::Error::checkPath($self->{"software:samtools"});
	$bowtie_cmd .= qq(export samtools="$samtools"\n);
	my $reference=GATE::Error::checkPath($self->{"database:$ref"});
	$bowtie_cmd.="export reference=$reference\n";
	my $bowtiepara=$self->{'setting:bowtie'};
	$bowtie_cmd.=qq(export bowtiepara="$bowtiepara"\n);
	my $bowtie_version=$1 if ($bowtie=~/([^\/]+)$/);
	$bowtie_version=$1 if ($bowtiepara=~/(bowtie[2]?)/);
	if (defined $bowtie_version && $bowtie !~ /$bowtie_version$/)
	{
		$bowtie =~ s/[^\/]+$/$bowtie_version/;
		$reference=$1 if ($reference=~/(\S+)\.fa/i);
	}
	my $bowtie_build="$bowtie-build";
	$bowtie_cmd.="export bowtie_build=$bowtie_build\n";
	$bowtie_cmd .= qq(export REFERENCE=$reference\n);
	my $workdir=GATE::Error::checkPath($self->{"-workdir"});
	$bowtie_cmd.="export workdir=$workdir\n";
	my $alndir=GATE::Error::checkPath($self->{"setting:aln_outdir"});
	$bowtie_cmd .= qq(export alndir="$alndir"\n);
	$bowtie_cmd.=qq(cd \$workdir\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	$bowtie_cmd .= qq([[ -d \${alndir} ]] || mkdir \${alndir}\n) if (!-d qq($self->{"-workdir"}/$alndir));;
	$bowtie_cmd .= "cd \${alndir}\n";
#bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]
#bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S [<hit>]
	foreach my $lib(@libraries) {
		$bowtie_cmd .= qq(echo `date`; echo "$lib"\n);
		$bowtie_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$alndir/$lib));
		$bowtie_cmd .= "cd $lib\n";
		my %fq=getlibInput($self->{"LIB"}{$lib});
		my @bam=();
		if (exists $fq{1} && exists $fq{2}){
			my @fq1=@{$fq{1}};
			my @fq2=@{$fq{2}};
			my $fq1=join ",",@fq1;
			my $fq2=join ",",@fq2;
			my @rg=();
			my $ID=(exists $self->{$lib}{'fq1'}{0}{'ID'})?$self->{$lib}{'fq1'}{0}{'ID'}:$lib;
			my $SM=(exists $self->{$lib}{'fq1'}{0}{'SM'})?$self->{$lib}{'fq1'}{0}{'SM'}:$lib;
			push @rg,"SM:$SM";
			my $LB=(exists $self->{$lib}{'fq1'}{0}{'LB'})?$self->{$lib}{'fq1'}{0}{'LB'}:$lib;
			push @rg,"LB:$LB";
			my $PL=(exists $self->{$lib}{'fq1'}{0}{'PL'})?$self->{$lib}{'fq1'}{0}{'PL'}:"ILLUMINA";
			push @rg,"PL:$LB";
			foreach my $rb(keys %{$self->{$lib}{'fq'}{0}})
			{
				if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
				{
					push @rg,qq($rb:$self->{$lib}{'fq'}{0}{$rb}');
				}
			}
			if ($bowtie_version=~/2$/) {
				#--rg-id HWI-SN957Lane7 --rg SM:ZEN456A1 --rg LB:ZEN456A1LI5 --rg PI:400 --rg PL:ILLUMINA
				my $samrg="--rg-id $ID ".(join "--rg ",@rg);
				$bowtie_cmd .= qq(\${bowties} \${bowtiepara} $samrg -x \$REFERENCE -1 $fq1 -2 $fq2 -S | \${samtools} -Sbh - -o $lib.pair.bam\n);
				push @bam,qq($workdir/$alndir/$lib/$lib.pair.bam);
			} else {
				#--sam-RG ID:HWI-SN957Lane7 --sam-RG SM:ZEN456A1 --sam-RG LB:ZEN456A1LI5 --sam-RG PI:400 --sam-RG PL:ILLUMINA
				my $samrg="--sam-RG $ID ".(join "--sam-RG ",@rg);
				$bowtie_cmd .= qq(\${bowties} \${bowtiepara} $samrg \$REFERENCE -1 $fq1 -2 $fq2 | \${samtools} -Sbh - -o $lib.pair.bam\n);
				push @bam,qq($workdir/$alndir/$lib/$lib.pair.bam);
			}
		}
		if (exists $fq{0} && @{$fq{0}}>0){
			my $fq=join ",",@{$fq{0}};
			my @rg="";
			my $ID=(exists $self->{$lib}{'fq'}{0}{'ID'})?$self->{$lib}{'fq'}{0}{'ID'}:$lib;
			my $SM=(exists $self->{$lib}{'fq'}{0}{'SM'})?$self->{$lib}{'fq'}{0}{'SM'}:$lib;
			push @rg,"SM:$SM";
			my $LB=(exists $self->{$lib}{'fq'}{0}{'LB'})?$self->{$lib}{'fq'}{0}{'LB'}:$lib;
			push @rg,"LB:$LB";
			my $PL=(exists $self->{$lib}{'fq'}{0}{'PL'})?$self->{$lib}{'fq'}{0}{'PL'}:"ILLUMINA";
			push @rg,"PL:$LB";
			foreach my $rb(keys %{$self->{$lib}{'fq'}{0}})
			{
				if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
				{
					push @rg,qq(\\t$rb:$self->{$lib}{'fq'}{0}{$rb}');
				}
			}
			if ($bowtie_version=~/2$/) {
				my $samrg="--rg-id $ID ".(join "--rg ",@rg);
				$bowtie_cmd .= qq(\${bowties} \${bowtiepara} $samrg -x \$REFERENCE -U $fq | \${samtools} -Sbh - -o $lib.single.bam\n);
				push @bam,qq($workdir/$alndir/$lib/$lib.single.bam);
			} else {
				my $samrg="--sam-RG $ID ".(join "--sam-RG ",@rg);
				$bowtie_cmd .= qq(\${bowties} \${bowtiepara} $samrg \$REFERENCE $fq | \${samtools} -Sbh - -o $lib.single.bam\n);
				push @bam,qq($workdir/$alndir/$lib/$lib.single.bam);
			}
		}
		if (@bam>1)
		{
			if (exists $self->{"software:picard"}) {
				my $picardpath=$1 if ($self->{"software:picard"}=~/(.*)\/[^\/\s]+$/);
				my $MergeSamFiles="$picardpath/MergeSamFiles.jar";
				$MergeSamFiles=qq(java -Xmx\${heap} -Djava.io.tmpdir=./tmp_merge -jar $MergeSamFiles) if ($MergeSamFiles!~/^java/ && $MergeSamFiles !~ /\-jar/);
				my $MergeSamFilesPara=$self->{"setting:MergeSamFiles"};
				my $MarkDuplicates=(exists $self->{"software:MarkDuplicates"})?$self->{"software:MarkDuplicates"}:"$picardpath/MarkDuplicates.jar";
				$MarkDuplicates=correctJavaCmd($MarkDuplicates,"\${heap}","./tmp_rmdup");
				my $MarkDuplicatesPara=(exists $self->{"setting:MarkDuplicates"})?$self->{"setting:MarkDuplicates"}:'VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true ASSUME_SORTED=true';
				my $merge_bam=join " INPUT=",@bam;
				$bowtie_cmd .= "$MergeSamFiles INPUT=$merge_bam $MergeSamFilesPara OUTPUT=$lib.merge.bam\n";
				$bowtie_cmd .= "\${samtools} index $lib.merge.bam\n";
				my $bowtiebam="$lib.merge.bam";
				if (!exists $self->{'rule:rmdup'} || GATE::Error::boolean($self->{'rule:rmdup'}==1)) {
					$bowtie_cmd .=  qq($MarkDuplicates INPUT=$bowtiebam OUTPUT=$lib.merge.rmdup.bam M=$lib.duplicate_report.txt $MarkDuplicatesPara);
					$bowtie_cmd .= "\${samtools} sort $lib.merge.rmdup.bam $lib.merge.rmdup.sort\n";
					$bowtie_cmd .= "\${samtools} index $lib.merge.rmdup.sort.bam\n";
					$bowtiebam = "$lib.merge.rmdup.sort.bam";
				}
				$bowtie_cmd .= "rm -rf ./tmp_merge\n" if (exists $self->{"rule:Clean"});
				@{$self->{$lib}{"$ref-bowtiebam"}}=();
				push @{$self->{$lib}{"$ref-bowtiebam"}},"$workdir/$alndir/$lib/$bowtiebam";
			} else {
				my $merge_bam=join " ",@bam;
				$bowtie_cmd .= qq(\${samtools} view -H $bam[0] |grep -v "^\@RG" | grep -v "^\@PG" >> $lib.inh.sam\n);
				foreach my $tophatbam(@bam)
				{
					$bowtie_cmd .=  qq(\${samtools} view -H $tophatbam |grep RG >> $lib.inh.sam\n);
				}
				$bowtie_cmd .=  qq(\${samtools} view -H $bam[0] |grep PG >> $lib.inh.sam\n);
				$bowtie_cmd .=  "\${samtools} merge -f -nr -h $lib.inh.sam $lib.merge.bam $merge_bam\n";
				my $bowtiebam="$lib.merge.bam";
				if (!exists $self->{'rule:rmdup'} || GATE::Error::boolean($self->{'rule:rmdup'}==1)) {
					$bowtie_cmd .= "\${samtools} rmdup $lib.merge.bam - |\${samtools} rmdup -S - - | \${samtools} sort -@ 8 -m 1G - $lib.merge.rmdup.sort\n";
					$bowtie_cmd .= "\${samtools} index $lib.merge.rmdup.sort.bam\n";
					$bowtiebam = "$lib.merge.rmdup.sort.bam";
				}
				$bowtie_cmd .= "rm -rf ./tmp_merge\n" if (exists $self->{"rule:Clean"});
				@{$self->{$lib}{"$ref-bowtiebam"}}=();
				push @{$self->{$lib}{"$ref-bowtiebam"}},"$workdir/$alndir/$lib/$bowtiebam";
			}
		}
		else
		{
			@{$self->{$lib}{"$ref-bowtiebam"}}=@bam;
		}
	}
	$bowtie_cmd .= "cd ..\n";
	return ($bowtie_cmd) unless (exists $self->{'rule:skip'}{'aln'} && GATE::Error::boolean($self->{'rule:skip'}{'aln'})==1);
}

## SOAPaligner v2.21
sub runSOAP ($$) {
#<ExecutablePath>/2bwt-builder <FastaPath/YourFasta>
#./soap –a <reads_a> -D <index.files> -o <output></output>
#./soap –a <reads_a> -b <reads_b> -D <index.files> -o <PE_output> -2 <SE_output> -m <min_insert_size> -x <max_insert_size>
	my $self=shift;
	my $ref=shift;
	if (!exists $self->{"software:soap"} || !-e $self->{"software:soap"} || !defined $self->{"software:soap"}) {
		return "";
	}
	$ref ||= 'ref';
	my $soap_cmd = qq(echo `date`; echo "run SOAPaligner"\n);
	$soap_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	my $soap = GATE::Error::checkPath($self->{"software:soap"});
	my $soap_path = $1 if (defined $soap && $soap=~/(\S+)\/soap/);
	my $builder;
	if (defined $self->{"software:2bwt-builder"}) {
		$builder=GATE::Error::checkPath($self->{"software:2bwt-builder"});
	} else {
		$builder=GATE::Error::checkPath("$soap_path/2bwt-builder");
	}
	$soap_cmd .= qq(export soap="$soap"\n);
	$soap_cmd .= qq(export soapbuilder="$builder"\n);
	my $reference=GATE::Error::checkPath($self->{"database:$ref"});
	$soap_cmd .= "\${soapbuilder} \$reference\n" unless (GATE::Error::checkIndex('soap',$reference)==1);
	$soap_cmd .= "export REFERENCE=$reference\.index\.\n";
	my $soappara=$self->{'setting:soap'};
	$soap_cmd .= qq(export soappara="$soappara"\n);
	my $soap2sam = GATE::Error::checkPath($self->{"software:soap2sam.pl"}) if (exists $self->{"software:soap2sam.pl"});
	$soap_cmd .= qq(export soap2sam="$soap2sam") if (defined $soap2sam);
	my $samtools = GATE::Error::checkPath($self->{"software:samtools"});
	$soap_cmd .= qq(export samtools="$samtools"\n);
	$soap_cmd .= "\${samtools} faidx \$reference\n" unless (GATE::Error::checkIndex('samtools',$reference)==1);
	$soap_cmd .= qq(perl -e 'open \(FAI,"$reference.fai"\);while\(<FAI>\){chomp;my \@t=split;print "\\\@SQ\\tSN:\$t[0]\\tLN:\$t[1]\\n"}close FAI;" > $reference.inh.sam) if (defined $samtools);
	
	my $workdir=GATE::Error::checkPath($self->{"-workdir"});
	$soap_cmd.="export workdir=$workdir\n";
	$soap_cmd.=qq(cd \$workdir\n);
	my $soap_outdir = (exists $self->{"setting:soap_outdir"})? $self->{"setting:soap_outdir"} : (exists $self->{"setting:aln_outdir"}) ? $self->{"setting:aln_outdir"} : "soap";
	$soap_cmd .= "[[ -d $soap_outdir ]] || mkdir $soap_outdir\n" if (!-d qq($self->{"-workdir"}/$soap_outdir));
	$soap_cmd .= "cd $soap_outdir\n";
	my @libraries=sort keys %{$self->{'LIB'}};
	my @bam=();
	foreach my $lib(@libraries) {
		$soap_cmd .= "[[ -d $lib ]] || mkdir $lib\n" if (!-d qq($self->{"-workdir"}/$soap_outdir/$lib));
		$soap_cmd .= "cd $lib\n";
		my %fq=getlibInput($self->{"LIB"}{$lib});
		if (exists $fq{1} && exists $fq{2}){
			my @fq1=@{$fq{1}};
			my @fq2=@{$fq{2}};
			for (my $j=0;$j<@fq2;$j++) {
				my $reads_a = $fq1[$j];
				my $reads_b = $fq2[$j];

				my $PI=(exists $self->{$lib}{'fq1'}{$j}{'PI'})?$self->{$lib}{'fq1'}{$j}{'PI'}:500;
				my $m = $PI*0.9;
				my $x = $PI*1.1;
				if (@fq2>1) {
					my $k=$j+1;
					my $ID=(exists $self->{$lib}{'fq1'}{$j}{'ID'})?$self->{$lib}{'fq1'}{$j}{'ID'}:"$lib-$k";
					my $SM=(exists $self->{$lib}{'fq1'}{$j}{'SM'})?$self->{$lib}{'fq1'}{$j}{'SM'}:$lib;
					my $LB=(exists $self->{$lib}{'fq1'}{$j}{'LB'})?$self->{$lib}{'fq1'}{$j}{'LB'}:$lib;
					my $PL=(exists $self->{$lib}{'fq1'}{$j}{'PL'})?$self->{$lib}{'fq1'}{$j}{'PL'}:"ILLUMINA";
					my $rg=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
					foreach my $rb(keys %{$self->{$lib}{'fq'}{$j}})
					{
						if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
						{
							$rg=~s/\'$//;
							$rg.=qq(\\t$rb:$self->{$lib}{'fq'}{$j}{$rb}');
						}
					}
					$soap_cmd .= qq(\${soap} -a $reads_a -b $reads_b -D \${REFERENCE} -o $lib.$k.paired.soap -2 $lib.$k.unpaired.soap -m $m -x $x \{$soappara}\n);
					if (defined $soap2sam) {
						$soap_cmd .= qq(awk 'BEGIN{print $rg}{print}' $reference.ihg.sam > $lib.inh.sam\n);
						$soap_cmd .= qq(\${soap2sam} -p $lib.$k.paired.soap > $lib.$k.paired.sam\n);
						$soap_cmd .= qq(\${soap2sam} -p $lib.$k.unpaired.soap > $lib.$k.unpaired.sam\n);
						$soap_cmd .= qq(cat $lib.$k.paired.sam $lib.$k.unpaired.sam > $lib.$k.pair.sam\n);
						$soap_cmd .= qq(cat $lib.inh.sam $lib.$k.pair.sam | \${samtools} view -Sbh - -o $lib.$k.paired.bam && \${samtools} index $lib.$k.pair.bam && );
						$soap_cmd .= qq(\${samtools} sort -m 3000000000 $lib.$k.pair.bam $lib.$k.pair.sort && \${samtools} index $lib.$k.pair.sort.bam \n);
						push @bam,"$lib.$k.pair.sort.bam";
					}
				} else {
					my $ID=(exists $self->{$lib}{'fq1'}{$j}{'ID'})?$self->{$lib}{'fq1'}{$j}{'ID'}:$lib;
					my $SM=(exists $self->{$lib}{'fq1'}{$j}{'SM'})?$self->{$lib}{'fq1'}{$j}{'SM'}:$lib;
					my $LB=(exists $self->{$lib}{'fq1'}{$j}{'LB'})?$self->{$lib}{'fq1'}{$j}{'LB'}:$lib;
					my $PL=(exists $self->{$lib}{'fq1'}{$j}{'PL'})?$self->{$lib}{'fq1'}{$j}{'PL'}:"ILLUMINA";
					my $rg=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
					foreach my $rb(keys %{$self->{$lib}{'fq'}{$j}})
					{
						if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
						{
							$rg=~s/\'$//;
							$rg.=qq(\\t$rb:$self->{$lib}{'fq'}{$j}{$rb}');
						}
					}
					$soap_cmd .= qq(\${soap} -a $reads_a -b $reads_b -D \${REFERENCE} -o $lib.pair.soap -2 $lib.single.soap -m $m -x $x \{$soappara}\n);
					if (defined $soap2sam) {
						$soap_cmd .= qq(awk 'BEGIN{print $rg}{print}' $reference.ihg.sam > $lib.inh.sam\n);
						$soap_cmd .= qq(\${soap2sam} -p $lib.paired.soap > $lib.paired.sam\n);
						$soap_cmd .= qq(cat $lib.paired.sam $lib.unpaired.sam > $lib.pair.sam\n);
						$soap_cmd .= qq(cat $lib.inh.sam $lib.sam |samtools view -Sbh - -o $lib.paired.bam && samtools index $lib.paired.bam && );
						$soap_cmd .= qq(\${samtools} sort -m 3000000000 $lib.pair.bam $lib.pair.sort && \${samtools} index $lib.pair.sort.bam\n);
						push @bam,"$lib.pair.sort.bam";
					}

				}
			}
		}
		if (exists $fq{0} && @{$fq{0}}>0){
			my @fq=@{$fq{0}};
			for (my $j=0;$j<@fq;$j++) {
				my $reads_a = $fq{$j};
				if (@fq>1) {
					my $k=$j+1;
					my $ID=(exists $self->{$lib}{'fq'}{$j}{'ID'})?$self->{$lib}{'fq'}{$j}{'ID'}:"$lib-$k";
					my $SM=(exists $self->{$lib}{'fq'}{$j}{'SM'})?$self->{$lib}{'fq'}{$j}{'SM'}:$lib;
					my $LB=(exists $self->{$lib}{'fq'}{$j}{'LB'})?$self->{$lib}{'fq'}{$j}{'LB'}:$lib;
					my $PL=(exists $self->{$lib}{'fq'}{$j}{'PL'})?$self->{$lib}{'fq'}{$j}{'PL'}:"ILLUMINA";
					my $rg=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
					foreach my $rb(keys %{$self->{$lib}{'fq'}{$j}})
					{
						if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
						{
							$rg=~s/\'$//;
							$rg.=qq(\\t$rb:$self->{$lib}{'fq'}{$j}{$rb}');
						}
					}
					$soap_cmd .= qq(\${soap} -a $reads_a -D \${REFERENCE} -o $lib.$k.single.soap \{$soappara}\n);
					if (defined $soap2sam) {
						$soap_cmd .= qq(awk 'BEGIN{print $rg}{print}' $reference.ihg.sam > $lib.inh.sam\n);
						$soap_cmd .= qq(\${soap2sam} -p $lib.$k.single.soap > $lib.$k.single.sam\n);
						$soap_cmd .= qq(cat $lib.inh.sam $lib.$k.single.sam |samtools view -Sbh - -o $lib.$k.single.bam && samtools index $lib.$k.single.bam && ); 
						$soap_cmd .= qq(\${samtools} sort -m 3000000000 $lib.$k.single.bam $lib.$k.single.sort && \${samtools} index $lib.$k.single.sort.bam\n);
						push @bam,"$lib.$k.single.sort.bam";
					}
				} else {
					my $ID=(exists $self->{$lib}{'fq'}{$j}{'ID'})?$self->{$lib}{'fq'}{$j}{'ID'}:$lib;
					my $SM=(exists $self->{$lib}{'fq'}{$j}{'SM'})?$self->{$lib}{'fq'}{$j}{'SM'}:$lib;
					my $LB=(exists $self->{$lib}{'fq'}{$j}{'LB'})?$self->{$lib}{'fq'}{$j}{'LB'}:$lib;
					my $PL=(exists $self->{$lib}{'fq'}{$j}{'PL'})?$self->{$lib}{'fq'}{$j}{'PL'}:"ILLUMINA";
					my $rg=qq('\@RG\\tID:$ID\\tPL:$PL\\tLB:$LB\\tSM:$SM');
					foreach my $rb(keys %{$self->{$lib}{'fq'}{$j}})
					{
						if ($rb=~/^([A-Z]{2})$/ && $rb ne 'ID' && $rb ne 'SM' && $rb ne 'LB' && $rb ne 'PL')
						{
							$rg=~s/\'$//;
							$rg.=qq(\\t$rb:$self->{$lib}{'fq'}{$j}{$rb}');
						}
					}
					$soap_cmd .= qq(\${soap} -a $reads_a -D \${REFERENCE} -o $lib.single.soap \{$soappara}\n);
					if (defined $soap2sam) {
						$soap_cmd .= qq(awk 'BEGIN{print $rg}{print}' $reference.ihg.sam > $lib.inh.sam\n);
						$soap_cmd .= qq(\${soap2sam} -p $lib.single.soap > $lib.single.sam\n);
						$soap_cmd .= qq(cat $lib.inh.sam $lib.single.sam |samtools view -Sbh - -o $lib.single.bam && samtools index $lib.single.bam && );
						$soap_cmd .= qq(\${samtools} sort -m 3000000000 $lib.single.bam $lib.single.sort && \${samtools} index $lib.single.sort.bam\n);
						push @bam,"$lib.single.sort.bam";
					}
				}
			}
			if (@bam>1) {
				my $merge_bam="";
				my $MergeSamFiles;
				if (exists $self->{"software:picard"}) {
					$MergeSamFiles="$1/MergeSamFiles.jar" if ($self->{"software:picard"}=~/(.*)\/[^\/\s]+$/);
					$MergeSamFiles=qq(java -Xmx\${heap} -Djava.io.tmpdir=./tmp_merge -jar $MergeSamFiles) if ($MergeSamFiles!~/^java/ && $MergeSamFiles !~ /\-jar/);
					if (defined $MergeSamFiles && $MergeSamFiles ne "") {
						my $MergeSamFilesPara=(exists $self->{"setting:MergeSamFiles"})?$self->{"setting:MergeSamFiles"}:'USE_THREADING=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT';
						$merge_bam=join " INPUT\=",@bam;
						$soap_cmd .= "$MergeSamFiles INPUT\=$merge_bam $MergeSamFilesPara OUTPUT\=$lib.merge.bam\n";
						if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1 )
						{
							$soap_cmd .= qq(rm -rf ./tmp_merge);
							$soap_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
						}
						$self->{"$ref-bam"}="$workdir/$soap_outdir/$lib/$lib.merge.bam";
					}
				} 
				else {
					$merge_bam=join " ",@bam;
					$soap_cmd .=  qq(\${samtools} view -H $bam[0] |grep -v "^\@RG" | grep -v "^\@PG" > $lib.inh.sam);

					foreach my $soapbam(@bam)
					{
						$soap_cmd .=  qq(\${samtools} view -H $soapbam |grep "^\@RG" >> $lib.inh.sam\n);
					}
					$soap_cmd .=  qq(\${samtools} view -H $bam[0] |grep "^\@PG" >> $lib.inh.sam);

					$soap_cmd .=  "\${samtools} merge -f -nr -h $lib.inh.sam $lib.merge.bam $merge_bam";
					$self->{"$ref-bam"}="$workdir/$soap_outdir/$lib/$lib.merge.bam";
				}
			} else {
				$self->{"$ref-bam"}="$workdir/$soap_outdir/$lib/$bam[0]";
			}
			
		}
		$soap_cmd .= "cd ../\n";
	}
	$soap_cmd .= "cd ..\n";
	return ($soap_cmd);	
}

sub runTopHat($$) {
	my $self=shift;
	my $ref=shift;
	if (!exists $self->{"software:tophat"} || !-e $self->{"software:tophat"} || !defined $self->{"software:bowtie"}) {
		return "";
	}
	$ref ||= 'ref';
	my $tophat_cmd = qq(echo `date`; echo "run TopHat"\n);
	$tophat_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	my $tophat=GATE::Error::checkPath($self->{'software:tophat'});
	$tophat_cmd.="export tophat=$tophat\n";
	my $bowtie=GATE::Error::checkPath($self->{'software:bowtie'});
	$tophat_cmd.="export bowtie=$bowtie\n";
	my $samtools = GATE::Error::checkPath($self->{"software:samtools"});
	$tophat_cmd .= qq(export samtools="$samtools"\n);
	my $reference=GATE::Error::checkPath($self->{"database:$ref"});
	$tophat_cmd.="export reference=$reference\n";
	#my $refGene=GATE::Error::checkPath($self->{'database:refGene'});
	#$tophat_cmd.="export refGene=$refGene\n";
	my $tophatpara=$self->{'setting:tophat'};
	$tophat_cmd.=qq(export tophatpara="$tophatpara"\n);
	my $bowtie_version=$1 if ($bowtie=~/([^\/]+)$/);
	$bowtie_version=$1 if ($tophatpara=~/(bowtie[2]?)/);
	if (defined $bowtie_version && $bowtie !~ /$bowtie_version$/)
	{
		$bowtie =~ s/[^\/]+$/$bowtie_version/;
	}
	my $bowtie_build="$bowtie-build";
	$tophat_cmd.="export bowtie_build=$bowtie_build\n";
	if ($bowtie=~/bowtie2/)
	{
		my $db_index=$1 if ($reference=~/(\S+)\.fa/i);
		my $dbpath=$1 if ($reference=~/(\S+)\/[^\/\s]+$/);
		$tophat_cmd .= "cd $dbpath\n";
		#if ($db_index=~/\S+\/([^\/\s]+)$/)
		#{
		#	$tophat_cmd .= "ln -s \$reference $1\n" unless (-f $db_index);
		#}
		$tophat_cmd .= "\${bowtie_build} \$reference $db_index\n" unless (GATE::Error::checkIndex('bowtie2',$db_index)==1);
		$reference=$db_index;
	}
	$tophat_cmd .= qq(export REFERENCE=$reference\n);
	my $workdir=GATE::Error::checkPath($self->{"-workdir"});
	$tophat_cmd.="export workdir=$workdir\n";
	$tophat_cmd.=qq(cd \$workdir\n);
	
	my @libraries=sort keys %{$self->{'LIB'}};
	my $tophat_outdir = (exists $self->{"setting:tophat_outdir"})? $self->{"setting:tophat_outdir"} : "tophat";
	$tophat_cmd .= "[[ -d $tophat_outdir ]] || mkdir $tophat_outdir\n" if (!-d qq($self->{"-workdir"}/$tophat_outdir));
	$tophat_cmd .= "cd $tophat_outdir\n";
	foreach my $lib(@libraries) {
		$tophat_cmd .= qq(echo `date`; echo "$lib"\n);
		#$tophat_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/tophat/$lib));
		#$tophat_cmd .= "cd $lib\n";
		my %fq=getlibInput($self->{"LIB"}{$lib});
		my @bam=();
		if (exists $fq{1} && exists $fq{2}){
			my @fq1=@{$fq{1}};
			my @fq2=@{$fq{2}};
			my $fq1=join ",",@fq1;
			my $fq2=join ",",@fq2;
			my $ID=(exists $self->{$lib}{'fq1'}{0}{'ID'})?$self->{$lib}{'fq1'}{0}{'ID'}:$lib;
			my $SM=(exists $self->{$lib}{'fq1'}{0}{'SM'})?$self->{$lib}{'fq1'}{0}{'SM'}:$lib;
			my $LB=(exists $self->{$lib}{'fq1'}{0}{'LB'})?$self->{$lib}{'fq1'}{0}{'LB'}:$lib;
			my $PL=(exists $self->{$lib}{'fq1'}{0}{'PL'})?$self->{$lib}{'fq1'}{0}{'PL'}:"ILLUMINA";
			my $PU=$self->{$lib}{'fq1'}{0}{'PU'} if (exists $self->{$lib}{'fq1'}{0}{'PU'});
			my $rg="--rg-id $ID --rg-platform $PL --rg-library $LB --rg-sample $SM";
			$rg.=" --rg-platform-unit $PU" if (defined $PU);
			$tophat_cmd .= "\${tophat} \${tophatpara} $rg -o $lib\_pe\_tophat \$REFERENCE $fq1 $fq2\n";
			$tophat_cmd .= "\${samtools} index $lib\_pe\_tophat/accepted_hits.bam\n";
			push @bam,$self->{'-workdir'}."/tophat/$lib\_pe\_tophat/accepted_hits.bam";
		}
		if (exists $fq{0} && @{$fq{0}}>0){
			my @fq1=@{$fq{0}};
			my $fq1=join ",",@fq1;
			my $ID=(exists $self->{$lib}{'fq'}{0}{'ID'})?$self->{$lib}{'fq'}{0}{'ID'}:$lib;
			my $SM=(exists $self->{$lib}{'fq'}{0}{'SM'})?$self->{$lib}{'fq'}{0}{'SM'}:$lib;
			my $LB=(exists $self->{$lib}{'fq'}{0}{'LB'})?$self->{$lib}{'fq'}{0}{'LB'}:$lib;
			my $PL=(exists $self->{$lib}{'fq'}{0}{'PL'})?$self->{$lib}{'fq'}{0}{'PL'}:"ILLUMINA";
			my $PU=$self->{$lib}{'fq'}{0}{'PU'} if (exists $self->{$lib}{'fq'}{0}{'PU'});
			my $rg="--rg-id $ID --rg-platform $PL --rg-library $LB --rg-sample $SM";
			$rg.=" --rg-platform-unit $PU" if (defined $PU);
			$tophat_cmd .= "\${tophat} \${tophatpara} $rg -o $lib\_se\_tophat \$REFERENCE $fq1\n";
			$tophat_cmd .= "\${samtools} index $lib\_se\_tophat/accepted_hits.bam\n";
			push @bam,$self->{'-workdir'}."/tophat/$lib\_se\_tophat/accepted_hits.bam";
		}
		push @bam, @{$self->{$lib}{"$ref-bwabam"}} if (exists $self->{$lib}{"$ref-bwabam"});
		if (@bam>1)
		{
			if (exists $self->{"software:picard"}) {
				my $MergeSamFiles="$1/MergeSamFiles.jar" if ($self->{"software:picard"}=~/(.*)\/[^\/\s]+$/);
				$MergeSamFiles=qq(java -Xmx4g -Djava.io.tmpdir=./tmp_merge -jar $MergeSamFiles) if ($MergeSamFiles!~/^java/ && $MergeSamFiles !~ /\-jar/);
				my $MergeSamFilesPara=$self->{"setting:MergeSamFiles"};
				my $merge_bam=join " INPUT=",@bam;
				$tophat_cmd .= "$MergeSamFiles INPUT=$merge_bam $MergeSamFilesPara OUTPUT=$lib.tophat.merge.bam\n";
				$tophat_cmd .= "\${samtools} index $lib.merge.bam\n";
				@{$self->{$lib}{"tophatbam"}}=();
				$self->{$lib}{"tophatbam"}=$self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.tophat.merge.bam";
				if (exists $self->{$lib}{"$ref-bwabam"})
				{
					@{$self->{$lib}{"$ref-bwabam"}}=();
					push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.tophat.merge.bam";
				}
			}
			else {
				my $merge_bam=join " ",@bam;
				$tophat_cmd .= qq(\${samtools} view -H $bam[0] |grep -v "^\@RG" | grep -v "^\@PG" >> $lib.inh.sam\n);
				foreach my $tophatbam(@bam)
				{
					$tophat_cmd .=  qq(\${samtools} view -H $tophatbam |grep RG >> $lib.inh.sam\n);
				}
				$tophat_cmd .=  qq(\${samtools} view -H $bam[0] |grep PG >> $lib.inh.sam\n);
				$tophat_cmd .=  "\${samtools} merge -f -nr -h $lib.inh.sam $lib.tophat.merge.bam $merge_bam\n";
				$self->{$lib}{"tophatbam"}=$self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.tophat.merge.bam";
				if (exists $self->{$lib}{"$ref-bwabam"})
				{
					@{$self->{$lib}{"$ref-bwabam"}}=();
					push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"setting:aln_outdir"}."/$lib/$lib.tophat.merge.bam";
				}
			}
		}
		else
		{
			$self->{$lib}{"tophatbam"}=$bam[0];
		}
	}
	$tophat_cmd .= "cd ..\n";
	return ($tophat_cmd);
}

sub runGSNAP ($) {
	
}

sub runGMAP ($) {
	
}

sub runSTAR ($) {
	
}

sub runBLAST ($) {
	
}

#psLayout version 3
#match  mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts  tStarts
#       match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
sub runBLAT ($) {
	my $self = shift;
	my $ref = shift;
	$ref ||= 'ref';
	my $blat = GATE::Error::checkPath($self->{"software:blat"});
	if (!defiined $blat) {
		return "";
	}
	my $blat_cmd = qq(echo `date`; echo "run blat"\n);
	$blat_cmd .= qq(export PATH="$self->{"setting:PATH"}":\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	my $reference = GATE::Error::checkPath($self->{"database:$ref"});
	$blat_cmd .= qq(export REFERENCE="$reference"\n);
	$blat_cmd .= qq(export blat="$blat"\n);
	my $blatpara = $self->{'setting:blat'};
	$blat_cmd .= qq(export blatpara="$blatpara"\n);
	my $workdir=GATE::Error::checkPath($self->{"-workdir"});
	$blat_cmd .= qq(export workdir="$workdir"\n);
	$blat_cmd .= qq(cd \${workdir}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
}

sub runLAST($) {
	
}

sub runLASTZ($) {
	
}

#########################################################
#                                                       #
#                      Genotyping                       #
#                                                       #
#########################################################

## GATK v2.7-2-g6bda569
## Picard v1.97
sub runGATK ($$) {
	my $self=shift;
	my $ref=shift;
	$ref ||= 'ref';
	my $gatk_cmd = qq(echo `date`; echo "run Samtools-Picard-GATK"\n);
	$gatk_cmd .= qq(export PATH="$self->{"setting:PATH"}":\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	my $reference=GATE::Error::checkPath($self->{"database:$ref"});
	$gatk_cmd .= qq(export REFERENCE="$reference"\n);
	my $heap=$self->{"setting:heap"};
	$gatk_cmd .= qq(export heap="$heap"\n);
	my $multi_t=$self->{"setting:multithreads"};
	$gatk_cmd .= qq(export multithreads=$multi_t\n);
	my $samtools = GATE::Error::checkPath($self->{"software:samtools"});
	$gatk_cmd .= qq(export samtools="$samtools"\n);
	$gatk_cmd .= "\${samtools} faidx \$REFERENCE\n" unless (GATE::Error::checkIndex('samtools',$reference)==1);
	my $samtools_path=$1 if ($samtools =~ /(\S+)\/samtools/);
	my ($picard,$picardpath,$MergeSamFiles,$MarkDuplicates,$bamtools);
	if (defined $samtools_path)
	{
		$gatk_cmd .= qq(export samtools_path="$samtools_path"\n);
		$gatk_cmd .= qq (export PATH="\${samtools_path}":"\${samtools_path}/bcftools":\$PATH\n);
	}
	if (exists $self->{"software:picard"})
	{
		$picard=$self->{"software:picard"};
		$picardpath=$1 if ($self->{"software:picard"}=~/(.*)\/[^\/\s]+$/);
		$gatk_cmd .= qq(export picard="$picard"\n);
		$picard = correctJavaCmd($picard,"\${heap}");
		$MergeSamFiles=(exists $self->{"software:MergeSamFiles"})?$self->{"software:MergeSamFiles"}:qq($picardpath/MergeSamFiles.jar);
		$gatk_cmd .= qq(export MergeSamFiles="$MergeSamFiles"\n);
		$MergeSamFiles=correctJavaCmd($MergeSamFiles,"\${heap}","./tmp_merge");
		$MarkDuplicates=(exists $self->{"software:MarkDuplicates"})?$self->{"software:MarkDuplicates"}:"$picardpath/MarkDuplicates.jar";
		$gatk_cmd .= qq(export MarkDuplicates="$MarkDuplicates"\n);
		$MarkDuplicates=correctJavaCmd($MarkDuplicates,"\${heap}","./tmp_rmdup");
		if (defined $picard) {
			my $dict=GATE::Error::checkIndex('picard',$reference);
			if ($dict==0) {
				my $CreateSequenceDictionary="";
				my $picardpath=$1 if ($self->{"software:picard"}=~/(.*)\/[^\/\s]+$/);
				$CreateSequenceDictionary=(exists $self->{"software:CreateSequenceDictionary"})?$self->{"software:CreateSequenceDictionary"}:qq($picardpath/CreateSequenceDictionary.jar);
				$CreateSequenceDictionary=correctJavaCmd($CreateSequenceDictionary,"\${heap}");
				my $ref_prefix=$1 if ($reference=~/([^\s]+)\.fa/i);
				$gatk_cmd .= "$CreateSequenceDictionary R=\$REFERENCE O=$ref_prefix\.dict CREATE_INDEX=true CREATE_MD5_FILE=true TRUNCATE_NAMES_AT_WHITESPACE=true VALIDATION_STRINGENCY=STRICT\n";
			}
		}
	}
	if (defined $self->{"software:bamtools"})
	{
		$bamtools=GATE::Error::checkPath($self->{"software:bamtools"});
		$gatk_cmd .= qq(export bamtools\="$bamtools"\n);
	}
	my $gatk=GATE::Error::checkPath($self->{"software:gatk"}) if (exists $self->{"software:gatk"});
	$gatk=correctJavaCmd($gatk,$heap,"./tmp_gatk");
	$gatk_cmd .= qq(export gatk="$gatk"\n);
	if (exists $self->{"database:dbSNP"})
	{
		$gatk_cmd .=qq(export dbSNP="$self->{"database:dbSNP"}"\n);
	}
	my $workdir=GATE::Error::checkPath($self->{"-workdir"});
	my $vardir=GATE::Error::checkPath($self->{"setting:var_outdir"});
	$gatk_cmd .= qq(export workdir="$workdir"\n);
	$gatk_cmd .= qq(export vardir="$vardir"\n);
	$gatk_cmd .= qq(cd \${workdir}\n);
	$gatk_cmd .= qq([[ -d \${vardir} ]] || mkdir -p \${vardir}\n) if (!-d qq($self->{"-workdir"}/$vardir));
	$gatk_cmd .= qq(cd \${vardir}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	my $multi=0;
	foreach my $lib(@libraries) {
		my $bam="";
		$gatk_cmd .= qq(echo `date`; echo "$lib"\n);
		if (!-e qq($self->{"-workdir"}/$self->{"setting:var_outdir"}/$lib))
		{
			$gatk_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) ;
		}
		$gatk_cmd .= "cd $lib\n";
		if (!exists $self->{$lib}{"$ref-bam"}) {
			if (exists $self->{$lib}{"$ref-bwabam"} && @{$self->{$lib}{"$ref-bwabam"}}>1) {
				my $merge_bam="";
				if (defined $MergeSamFiles && $MergeSamFiles ne "") {
## Merge BAM 
					#$gatk_cmd .= qq(mkdir -p ./tmp_merge);
					#$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
					#$gatk_cmd .= qq(export tmp_merge="./tmpmerge");
					#$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
					my $MergeSamFilesPara=(exists $self->{"setting:MergeSamFiles"})?$self->{"setting:MergeSamFiles"}:'USE_THREADING=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT';
					$merge_bam=join " INPUT\=",@{$self->{$lib}{"ref-bwabam"}};
					$gatk_cmd .= "$MergeSamFiles INPUT\=$merge_bam $MergeSamFilesPara OUTPUT\=$lib.merge.bam";
					$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
					if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1 )
					{
						$gatk_cmd .= qq(rm -rf ./tmp_merge);
						$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
					}
					$bam="$lib.merge.bam";
				} else {
					$merge_bam=join " ",@{$self->{$lib}{"ref-bwabam"}};
					$gatk_cmd .=  qq(\${samtools} view -H ${$self->{$lib}{"ref-bwabam"}}[0] |grep -v "^\@RG" | grep -v "^\@PG" >> $lib.inh.sam);
					$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
					foreach my $bwabam(@{$self->{$lib}{"ref-bwabam"}})
					{
						$gatk_cmd .=  qq(\${samtools} view -H $bwabam |grep "^\@RG" >> $lib.inh.sam);
						$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
					}
					$gatk_cmd .=  qq(\${samtools} view -H ${$self->{$lib}{"ref-bwabam"}}[0] |grep "^\@PG" >> $lib.inh.sam);
					$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
					$gatk_cmd .=  "\${samtools} merge -f -nr -h $lib.inh.sam $lib.merge.bam $merge_bam";
					$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
					$bam="$lib.merge.bam";
				}
				$self->{$lib}{"$ref-bam"}=qq($self->{"-workdir"}/$self->{"setting:var_outdir"}/$lib/$lib.merge.bam);
			}
			elsif (exists $self->{$lib}{"$ref-bwabam"}) {
				$bam=${$self->{$lib}{"$ref-bwabam"}}[0];
				#$gatk_cmd .= "\${samtools} mpileup -ugf \$REFERENCE $bam | bcftools view -bvcg - | bcftools view -cg - > $lib.var.vcf && vcfutils.pl varFilter -D100 $lib.var.vcf > $lib.var.flt.vcf";
				#$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
				$self->{$lib}{"$ref-bam"}=${$self->{$lib}{"$ref-bwabam"}}[0];
			}
		} else {
			$bam=$self->{$lib}{"$ref-bam"};
		}

## Filter
#${bamtools} filter \
#  -isMapped true \
#  -isPaired true \
#  -isProperPair true \
#  -in  $PWDS/${subjectID}.fxmt.bam \
#  -out $PWDS/${subjectID}.fxmt.flt.bam

#${bamtools} stats \
#  -insert \
#  -in $PWDS/${subjectID}.fxmt.flt.bam \
#  > $PWDS/${subjectID}.fxmt.flt.stats
		if (defined $bamtools) {
			my $fltbam="$1.flt.bam" if ($bam=~/([^\/\s]+)\.bam$/);
			$gatk_cmd .= qq(\${bamtools} filter -isMapped true -isPaired true -in $bam -out $fltbam);
			$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			$bam="$lib.merge.flt.bam";
		}
## Remove duplicates
#java -Xmx${heap}m -Djava.io.tmpdir\=${tmp_folder}_rmdup \
#  -jar ${picard}/MarkDuplicates.jar \
#  I\=$PWDS/${subjectID}.fxmt.flt.bam \
#  O\=$PWDS/${subjectID}.rmdup.bam \
#  M\=$PWDS/${subjectID}.duplicate_report.txt \
#  VALIDATION_STRINGENCY\=SILENT \
#  REMOVE_DUPLICATES\=true
		if (!exists $self->{'rule:rmdup'} || GATE::Error::boolean($self->{'rule:rmdup'})==1) {
			if (defined $MarkDuplicates && $MarkDuplicates ne "") {
				#$gatk_cmd .= qq(mkdir -p ./tmp_rmdup);
				#$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
				#$gatk_cmd .= qq(export tmp_rmdup="./tmp_rmdup");
				#$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
				my $MarkDuplicatesPara=(exists $self->{"setting:MarkDuplicates"})?$self->{"setting:MarkDuplicates"}:'VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true ASSUME_SORTED=true';
				my $rmdupbam="$1.rmdup.bam" if ($bam=~/([^\/\s]+)\.bam/);
				$gatk_cmd .= qq($MarkDuplicates INPUT=$bam OUTPUT=$rmdupbam M=$lib.duplicate_report.txt $MarkDuplicatesPara);
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
				$gatk_cmd .= qq(\${samtools} index $rmdupbam);
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
				$bam=$rmdupbam;
				if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1 )
				{
					$gatk_cmd .= qq(rm -rf ./tmp_rmdup);
					$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
				}
			}else{
				my $rmdupbam="$1.rmdup.sort" if ($bam=~/([^\/\s]+)\.bam/);
				$gatk_cmd .= "\${samtools} rmdup $bam - | \${samtools} rmdup -S - - | \${samtools} sort -m 3000000000 - $rmdupbam";
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
				$gatk_cmd .= "\${samtools} index $rmdupbam.bam";
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
				$bam="$rmdupbam.bam";
			}
		}
## Stat
#${bamtools} stats \
#  -insert \
#  -in $PWDS/${subjectID}.rmdup.bam \
#  > $PWDS/${subjectID}.rmdup.stats

		if (defined $self->{"software:bamtools"}) {
				my $stats="$1.stats" if ($bam=~/([^\/\s]+)\.bam/);
				$gatk_cmd .= qq(\${bamtools} stats -insert -in $bam > $stats\n);
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
		} 

		if (defined $gatk) {
##Local Realignment
#------------------
#Samtools calls short indels with local realignment, but it does not write a modified BAM file after the realignment.
#The GATK though provides such a tool that realigns reads in regions with suspected indel artifacts and generates a BAM with cleaned alignments.
			#$gatk_cmd .= "mkdir ./tmp_realign";
			#$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			$gatk_cmd .= qq($gatk -T RealignerTargetCreator -R \$REFERENCE -I $bam -o $lib.gatk.intervals);
			if (exists $self->{"database:dbSNP"})
			{
				$gatk_cmd .= qq( --dbsnp \$dbSNP);
			}
			$gatk_cmd .= (exists $self->{"rule:multimode"}) ? qq( -nt $self->{"setting:multithreads"} && ) : "\n";
			my $realignedbam="$1.realigned.bam" if ($bam=~/([^\/\s]+)\.bam/);
			$gatk_cmd .= qq($gatk -T IndelRealigner -R \$REFERENCE -I $bam -o $realignedbam -targetIntervals $lib.gatk.intervals -LOD 0.4 -compress 6 -l INFO);
			if (exists $self->{"database:dbSNP"})
			{
				$gatk_cmd .= qq( --dbsnp \$dbSNP);
			}
			$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			if (defined $self->{"software:bamtools"})
			{
				my $stats="$1.stats" if ($realignedbam=~/([^\/\s]+)\.bam/);
				$gatk_cmd .= qq(\${bamtools} stats -insert -in $realignedbam  > $stats\n);
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			}
			$bam=$realignedbam;

## Analysis of covariatesDetermine the covariates affecting base quality scores in the BAM file.
# Previously, it was:
#java -Xmx${heap}m -Djava.io.tmpdir\=${tmp_folder}_covar \
#-jar $gatk \
#-R $REF \
#-l INFO \
#-I $PWDS/${subjectID}.realigned.srt.bam \
#-knownSites $DBSNP \
#-T CountCovariates \
#-nt $CPUs \
#-cov ReadGroupCovariate \
#-cov QualityScoreCovariate \
#-cov CycleCovariate \
#-cov DinucCovariate \
#-recalFile $PWDS/${subjectID}.flt.recal_v1.csv  \
#-L $ExonFile
#
# but now, it's rather
#java -Xmx4g -jar GenomeAnalysisTK.jar -l INFO -R hg19.fa -knowSites dbsnp132.txt -I input.marked.realigned.fixed.bam -T BaseRecalibrator -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate --out input.recal_data.csv
#From GATK doc:
#java -Xmx4g -jar GenomeAnalysisTK.jar \
#-T BaseRecalibrator \
#-I my_reads.bam \
#-R resources/Homo_sapiens_assembly18.fasta \
#-knownSites bundle/hg18/dbsnp_132.hg18.vcf \
#-knownSites another/optional/setOfSitesToMask.vcf \
#-o recal_data.grp
			#$gatk_cmd .= "mkdir ./tmp_covar";
			#$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			if (exists $self->{"database:dbSNP"})
			{
				$gatk_cmd .= qq($gatk -T BaseRecalibrator -R \$REFERENCE -knowSites \$dbSNP -l INFO -I $bam -o $lib.recalibration_report.grp -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate);
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";

## Generate AnlyzeCovariates Plots
#
#java -Xmx${heap}m -Djava.io.tmpdir\=${tmp_folder}_covar \
#-jar ${gatk_analyzecovar} \
#-recalFile $PWDS/${subjectID}.flt.recal_v1.csv  \
#-outputDir $PWDS/analyzeCovar_v1  \
#-ignoreQ 3
				my $AnalyzeCovariates;
				if (exists $self->{"software:AnalyzeCovariates"})
				{
					$AnalyzeCovariates=correctJavaCmd($AnalyzeCovariates,$heap,"./tmp_covar");
				}
				else
				{
					$AnalyzeCovariates=GATE::Error::checkPath("$1/resource/AnalyzeCovariates.jar") if ($self->{"software:gatk"}=~/(\S+)\/[^\/\s]+$/);
					if (defined $AnalyzeCovariates && -f $AnalyzeCovariates)
					{
						$AnalyzeCovariates=correctJavaCmd($AnalyzeCovariates,$heap,"./tmp_covar");
					}
				}
				if (defined $AnalyzeCovariates)
				{
					$gatk_cmd .= qq(mkdir analyzeCovar_v1);
					$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
					$gatk_cmd .= qq($AnalyzeCovariates -recalFile $lib.flt.recal_v1.csv -outputDrir analyzeCovar_v1 -ignoreQ 3);
					$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
				}

##Base Quality Recalibration
#---------------------------
#Base quality recalibration is optional, but it helps to improve the SNP quality and is preferred whenever possible.
#This tool recalibrates base quality scores of sequencing-by-synthesis reads in an aligned BAM file.
#Often base quality scores as reported by the sequencing machines can be inaccurate or uninformative.
#After recalibration, the quality scores in the QUAL field in each read in the output BAM are more accurate in that
#the reported quality score is much closer to its actual probability of mismatching the reference genome. Moreover,
#the recalibration tool attempts to correct for variation in quality with machine cycle and sequence context, and by
#doing so provides not only more accurate quality scores but also more widely dispersed (and thus more informative) ones.

## BaseRecalibrator
#First pass of the base quality score recalibration -- Generates recalibration table based on various user-specified
#covariates (such as read group, reported quality score, machine cycle, and nucleotide context).
			#$gatk_cmd .= qq($gatk -I $bam -R \$REFERENCE  -T BaseRecalibrator -l INFO \\\n);
			#$gatk_cmd .= qq(-recalFile recal_data.csv -cov ReadGroupCovariate -cov QualityScoreCovariate \\\n);
			#$gatk_cmd .= qq(-cov CycleCovariate -cov DinucCovariate -cov TileCovariate\n);
			#$gatk_cmd .= qq($gatk -I $bam -R \$REFERENCE -T TableRecalibration -l INFO \\\n);
			#$gatk_cmd .= qq(-recalFile recal_data.csv --output_bam $lib.merge.sort.alnRecal.bam\n);
			#$gatk_cmd .= "$samtools index $lib.merge.sort.alnRecal.bam\n";
			#$bam="$lib.merge.sort.alnRecal.bam";
#Creating a recalibrated BAM
#java -Xmx${heap}m -Djava.io.tmpdir\=${tmp_folder}_recal \
#-jar $gatk \
#-l INFO \
#-R $REF \
#-I $PWDS/${subjectID}.realigned.srt.bam \
#-T TableRecalibration \
#--default_platform Illumina \
# --default_read_group MP1  \
#--out $PWDS/${subjectID}.realigned.recal.bam \
#-recalFile $PWDS/${subjectID}.flt.recal_v1.csv  \
#-L $ExonFile

#To create a recalibrated BAM you can use GATK's PrintReads with the engine on-the-fly recalibration capability.
#Here is a typical command line to do so:
#java -jar GenomeAnalysisTK.jar \
#   -T PrintReads \
#   -R reference.fasta \
#   -I input.bam \
#   -BQSR recalibration_report.grp \
#   -o output.bam
				my $recalbam = "$1.recal.bam" if ($bam=~/([^\/\s]+)\.bam$/);
				$gatk_cmd .= qq($gatk -T PrintReads -R \$REFERENCE -I $bam -BQSR $lib.recalibration_report.grp -o $recalbam && $samtools index $lib.recal.bam);
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
				$bam = $recalbam;

## Re-analysis of covariatesDetermine the covariates affecting base quality scores
## in the realigned recalibrated BAM file for the comparison. This is an optional step. 
#java -Xmx${heap}m -Djava.io.tmpdir\=${tmp_folder}_covar \
# -jar $gatk \
# -R $REF \
# -I $PWDS/${subjectID}.realigned.recal.bam \
# -knownSites ${DBSNP} \
# -T CountCovariates \
# -nt $CPUs \
# -cov ReadGroupCovariate \
# -cov QualityScoreCovariate \
# -cov CycleCovariate \
# -cov DinucCovariate \
# -recalFile $PWDS/${subjectID}.flt.recal_v2.csv  \
# -L $ExonFile
				$gatk_cmd .= qq($gatk -T BaseRecalibrator -R \$REFERENCE -knowSites \$dbSNP -I $bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -o $lib.flt.recal_v2.csv);
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? qq( -nt $self->{"setting:multithreads"} && ) : "\n";
			}
			#else
			#{
			#	$gatk_cmd .= qq($gatk -T BaseRecalibrator -I $bam -R \$REFERENCE -run_without_dbsnp_potentially_ruining_quality -o $lib.recalibration_report.grp --intermediate_csv_file $lib.recal.csv --plot_pdf_file $lib.comp.pdf);
			#	$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			#	my $recalbam = "$1.recal.bam" if ($bam=~/([^\/\s]+)\.bam$/);
			#	$gatk_cmd .= qq($gatk -T PrintReads -R \$REFERENCE -I $bam -BQSR $lib.recalibration_report.grp -o $recalbam && $samtools index $lib.recal.bam);
			#	$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			#	$bam = $recalbam;
			#}

## BAQ calmd
			my $baqbam = "$1.baq.bam" if ($bam=~/([^\/\s]+)\.bam$/);
			$gatk_cmd .= "\${samtools} calmd -Abr $bam \$REFERENCE > $baqbam && \${samtools} index $baqbam";
			$bam = $baqbam;
			$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";

## Genotyping calling
			$gatk_cmd .= "$gatk -T UnifiedGenotyper -R \$REFERENCE -I $bam -baq CALCULATE_AS_NECESSARY -o $lib.gatk.var.vcf -U -S SILENT -rf BadCigar";
			$gatk_cmd .= (exists $self->{"rule:multimode"}) ? qq( -nt $self->{"setting:multithreads"} -nct $self->{"setting:multithreads"} && ) : "";
			$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			if (exists $self->{"software:tabix"} || $self->{"software:bgzip"})
			{
				my $bgzip;
				if (exists $self->{"software:bgzip"})
				{
					$bgzip=GATE::Error::checkPath($self->{"software:bgzip"});
				}
				else
				{
					my $tabix=GATE::Error::checkPath($self->{"software:tabix"}) if (exists $self->{"software:tabix"});
					$bgzip=$tabix;
					$bgzip=~s/tabix$/bgzip/;
				}
				$gatk_cmd .= "$bgzip $lib.gatk.var.vcf";
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			}

##coverage
#-------------------
#The GATK provides a tool to assess the depth of coverage over any number of intervals in a BAM file.
#The coverage will be displayed in the output for every single position plus as an average over each interval specified.
			if (exists $self->{"setting:TargetIntervalList"})
			{
				my $TL=checkPaht($self->{"setting:TargetIntervalList"});
				$gatk_cmd .= "$gatk -T DepthOfCoverage -I $bam -R \$REFERENCE -o $lib.coverage.depth -L $TL -ct 0 -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -ct 25 -ct 30 -pt readgroup --omitDepthOutputAtEachBase --omitIntervalStatistics --omitLocusTable";
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			} else {
				#my $output_prefix = $1 if ($bam=~/([^\/\s]+)\.bam/i);
				$gatk_cmd .= "$gatk -T DepthOfCoverage -I $bam -R \$REFERENCE -l INFO -o $lib.coverage.depth -ct 0 -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -ct 25 -ct 30 -pt readgroup --omitDepthOutputAtEachBase --omitIntervalStatistics --omitLocusTable";
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			}
			#$bam=$self->{'-workdir'}."/".$self->{"setting:var_outdir"}."/$lib/$lib.realigned.baq.bam";
			$self->{$lib}{"$ref-bam"}=$self->{'-workdir'}."/".$self->{"setting:var_outdir"}."/$lib/$bam";
			$self->{$lib}{"$ref-vcf"}=$self->{'-workdir'}."/".$self->{"setting:var_outdir"}."/$lib/$lib.gatk.var.vcf";
			if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1 )
			{
				$gatk_cmd .= "rm -rf $lib.merge.sort.bam $lib.realigned.bam tmp_realign tmp_gatk";
				$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			}
		} else {
			$gatk_cmd .= "\${samtools} mpileup -ugf \$REFERENCE $bam | bcftools view -bvcg - | bcftools view -cg - > $lib.var.vcf &&  vcfutils.pl varFilter -D100 $lib.var.vcf > $lib.var.flt.vcf";
			$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			$self->{$lib}{"$ref-bam"}=$self->{'-workdir'}."/".$self->{"setting:var_outdir"}."/$lib/$bam";
		}
		if (exists $self->{"software:sam2reads"}) {
			my $name=$1 if ($bam=~/([^\/\s]+).bam/);
			my $sam2reads=GATE::Error::checkPath($self->{"software:sam2reads"});
			$gatk_cmd .= "\${samtools} view -F 4 $bam | $sam2reads -R1 $name\_mapped.R1.fastq -R2 $name\_mapped.R2.fastq -R3 $name\_mapped.R3.fastq -f fq";
			$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
			$gatk_cmd .= "\${samtools} view -f 4 $bam | $sam2reads -R1 $name\_unmapped.R1.fastq -R2 $name\_unmapped.R2.fastq -R3 $name\_unmapped.R3.fastq -f fq";
			$gatk_cmd .= (exists $self->{"rule:multimode"}) ? " && " : "\n";
		}
		if (exists $self->{"overlap"} && exists $self->{"sam2bed"} && exists $self->{"msort"}) {
			$gatk_cmd .= $self->stat_mappedreads("bam",$bam,"lib",$lib);
		}
		if (exists $self->{"sam2bed"}){
			my $bed="$1.bed" if ($bam=~/([^\/\s]+)\.bam$/);
			$gatk_cmd .= qq( [[ -f $bed ]] || \${samtools} view -F4 $bam | $self->{"sam2bed"} -u -n > $bed\n);
			push @{$self->{'LIB'}{$lib}{'INPUT'}{'BED'}},$self->{'-workdir'}."/".$self->{"setting:var_outdir"}."/$lib/$bed";
		}
		$multi++;
		$gatk_cmd  =~ s/\&*\s*$//;
		$gatk_cmd .= ((exists $self->{"setting:multithreads"}) && ($multi % $self->{"setting:multithreads"}!=0) && ($multi<@libraries)) ? " &\n" : "\n";
		$gatk_cmd .= "cd ..\n";
	}
	$gatk_cmd  =~ s/\&*\s*$/\n/;
	$gatk_cmd .= "cd ..\n";
	return $gatk_cmd;
}

sub runDindel ($) {
#/usr/local/bin/dindel --analysis getCIGARindels --bamFile realigned.baq.bam --outputFile dindel_output --ref reference.fa
#python /usr/local/bin/makeWindows.py --inputVarFile dindel_output.variants.txt --windowFilePrefix realign_windows --numWindowsPerFile 1000
#perl -e 'my @f=glob("realign_windows.*.txt");foreach (@f){my $prefix="dindel_stage2_output_windows.$1" if ($_=~/windows\.(\d+)\./);system "/usr/local/bin/dindel --analysis indels --doDiploid --bamFile realigned.baq.bam --ref reference.fa --varFile $_ --libFile dindel_output.libraries.txt --outputFile $prefix";}'
#ls dindel_stage2_output_windows.*.glf.txt > dindel_stage2_outputfiles.txt
#python /usr/local/bin/mergeOutputDiploid.py --inputFiles dindel_stage2_outputfiles.txt --outputFile variantCalls.VCF --ref reference.fa 
	my $self=shift;
	my $ref=shift;
	$ref ||= 'ref';
	if (!exists $self->{"software:dindel"}) {
		return "";
	}
	my $dindel=GATE::Error::checkPath($self->{"software:dindel"});
	my $dindel_cmd = qq(echo `date`; echo "run dindel"\n);
	$dindel_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$dindel_cmd .= qq(export dindel=$dindel\n);
	my $para = $self->{"setting:dindel"};
	$dindel_cmd .= qq(export para="$para"\n);
	my $multithreads = $self->{"setting:multithreads"};
	my $multirun = $self->{"software:multithreads-run"};
	$dindel_cmd .= qq(export multirun=$multirun\n);
	my $reference = $self->{"database:$ref"};
	$dindel_cmd .= qq(export REFERENCE=$reference\n);
	my $dindel_output=(exists $self->{"setting:dindel_outdir"})?$self->{"setting:dindel_outdir"}:"dindel";
	$dindel_cmd .= "[[ -d $dindel_output ]] || mkdir $dindel_output\n";
	$dindel_cmd .= "cd $dindel_output\n";
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		if (exists $self->{$lib}{"$ref-bam"}) {
			$dindel_cmd .= "[[ -d $lib ]] || mkdir $lib\n";
			$dindel_cmd .= "cd $lib\n";
			my $bam=$self->{$lib}{"$ref-bam"};
			$dindel_cmd .= qq(\${dindel} --anlysis getCIGARindels --bamFile $bam --outputFile $lib.dindel_output --ref \${REFERENCE}\n);
			$dindel_cmd .= qq(makeWindows.py --inputVarFile $lib.dindel_output.variants.txt --windowFilePrefix $lib.realign_windows --numWindowsPerFile 1000\n);
			#/usr/local/bin/dindel --analysis indels --doDiploid --bamFile realigned.baq.bam --ref reference.fa --varFile $_ --libFile dindel_output.libraries.txt --outputFile $prefix
			$dindel_cmd .= qq(perl -e \'my \@f=glob\("$lib.realign_windows.*.txt"\););
			$dindel_cmd .= qq(for\(my \$i=0;\$i<\@f;\$i+=$multithreads\){my \@cmdary=();foreach my \$j\(\$i..\(\$i+$multithreads-1\)\){my \$prefix="$lib.dindel_stage2_output_windows.\$1" if \(\$f[\$i]=~/windows\\.\(\\d+\)\\./\););
			$dindel_cmd .= qq(push \@cmdary,qq\(\"$dindel --analysis indels --doDiploid --bamFile $bam --ref $reference --varFile \$f[\$i] --libFile $lib.dindel_output.libraries.txt --outputFile \$prefix\"\);}my \$cmd="$multirun ".join " ",\@cmdary;system \$cmd}\'\n);
			$dindel_cmd .= qq(ls $lib.dindel_stage2_output_windows.*.glf.txt > $lib.dindel_stage2_outputfiles.txt\n);
			$dindel_cmd .= qq(mergeOutputDiploid.py --inputFiles $lib.dindel_stage2_outputfiles.txt --outputFile $lib.variantCalls.VCF --ref \${REFERENCE}\n);
			$dindel_cmd .= qq(rm $lib.realign_windows.*.txt $lib.dindel_stage2_output_windows.*.glf.txt\n) if (exists $self->{"rule:Clean"} && GATE::Error::boolean($self->{"rule:Clean"})==1 );
			$dindel_cmd .= qq(cd ..\n);
		}
	}
	return $dindel_cmd;
}

sub runPindel ($) {
#mkdir output
#./pindel -f demo/hs_ref_chr20.fa -p demo/COLO-829_20-p_ok.txt -c 20 -o output/ref
#./pindel -f demo/simulated_reference.fa -i demo/simulated_config.txt -c ALL -o output/simulated
#./pindel -f <reference.fa> -p <pindel_input> [and/or -i bam_configuration_file] -c <chromosome_name> -o <prefix_for_output_files>
	
}

#http://vcftools.sourceforge.net/options.html
sub runVCFtools ($) {
	
}

# iSAAC
sub runiSAAC ($) {
	
}

sub runSNVer ($) {
##a) For individual sequencing data
# /path/to/java -jar /path/to/SNVer-0.2.0/SNVerIndividual.jar \
# -i pe.sorted.dedup.bam -o prefix_of_output -r ref.fasta -l target.bed

## b) For pooled Sequencing data
# /path/to/java -jar /path/to/SNVer-0.2.0/SNVerPool.jar -c pool.info \
# -i input_bam -o prefix_of_output -r ref.fasta -l target.bed
## or
# /path/to/java -jar /path/to/SNVer-0.2.0/SNVerPool.jar -n 96 \
# -i input_bam -o prefix_of_output -r ref.fasta -l target.bed
## Annotation
#/path/to/annovar/convert2annovar.pl -format vcf4 pe.vcf > input
#/path/to/annovar/summarize_annovar.pl --verdbsnp 132 --buildver hg19 \
#--outfile sum input /path/to/humandb 	
}

sub runSnpEff ($) {
#"A program for annotating and predicting the effects of single nucleotide polymorphisms,
#SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.", 
#Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 
#2012 Apr-Jun;6(2):80-92. PMID: 22728672 [PubMed - in process]
	my $self=shift;
}

sub runSnpSift ($) {
	my $self=shift;
}

sub runSOAPsnp($) {
#soapsnp -B aln.bam -d reference.fa -o cns -r 0.0005 -e 0.001 -u -L 150 -2 -Q J -s dbSNP -T region.out -m
#The dbSNP file consist of a lot of lines like this one:
#	chr1    201979756       1       1       0       0.161   0       0       0.839   rs568
#	The columns from left to right are: name of chromosome, coordinate on the chromosome, whether 
#	the SNP	has allele frequency information (1 is true, 0 is false), whether the SNP is validated 
#	by experiment (1 is true, 0 is false), whether the SNP is actually an indel (1 is true, 0 is false),
#	frequency of A, frequency of C, frequency of T, frequency of G, SNP id. For known SNP sites that do
#	not have allele frequency information, the frequency information can be arbitrarily determined as 
#	any positive values, which only imply what alleles have already been deposited in the database.
	my $self=shift;
}

sub checkdbSNP {
	my $dbSNP = shift;
	my $count=0;
	open (IN,$dbSNP) || die $!;
	while(<IN>){
		next if (/^\#/ || /^\s+/);
		my @t=split "\t+", $_;
		for (my $i=1;$i<$#t-1;$i++){
			$count++ if (/^\d+/);
		}
		last;
	}
	close IN;
	if ($count==8){
		return 1;
	} else {
		return 0;
	}
}

sub runSOAPsnv ($) {
	
}

sub runSOAPindel ($) {
	
}

sub runSOAPpopindel ($) {
	
}

sub runSOAPsv ($) {
	
}

sub runBreakDancer ($) {
	
}

sub runCRISP ($) {
	my $self = shift;
}

sub runPlink ($) {
	
}

sub runCNVNator ($) {
	
}

#########################################################
#                                                       #
#                       Assembly                        #
#                                                       #
#########################################################

sub runCufflinks($) {
	my ($self,%attrs)=@_;
	if (!exists $self->{"software:cufflinks"} || !defined $self->{"software:cufflinks"}) {
		return "";
	}
	my $cufflinks_cmd = qq(echo `date`; echo "run Cufflinks"\n);
	$cufflinks_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	my $cufflinks=GATE::Error::checkPath($self->{'software:cufflinks'});
	$cufflinks_cmd .= qq(export cufflinks=$cufflinks\n);
	my $gffread=GATE::Error::checkPath($self->{'software:gffread'}) if (exists $self->{'software:gffread'});
	$cufflinks_cmd .= qq(export gffread=$gffread\n) if (defined $gffread);
	my $bowtie=GATE::Error::checkPath($self->{'software:bowtie'});
	$cufflinks_cmd .= qq(export bowtie=$bowtie\n);
	my $reference=GATE::Error::checkPath($self->{'database:ref'});
	$cufflinks_cmd .= qq(export reference=$reference\n);
	my $refGene=GATE::Error::checkPath($self->{'database:refGene'});
	if (defined $gffread)
	{
		if ($refGene=~/\.gff/i)
		{
			my $gff=$refGene;
			$gff=~s/gff/cufflinks\.gff/;
			if (!-f $gff)
			{
				$cufflinks_cmd .= qq(\${gffread} -E $refGene -o- > $gff 2> gffread.log\n);
			}
			$refGene=$gff;
		}
		elsif ($refGene=~/\.gtf/i)
		{
			my $gtf=$reference;
			$gtf=~s/gff/cufflinks\.gff/;
			if (-f $gtf)
			{
				$cufflinks_cmd .= qq(\${gffread}  -E $refGene -T -o- > $gtf 2> gffread.log\n);
			}
			$refGene=$gtf;
		}
	}
	$cufflinks_cmd .= qq(export refGene=$refGene\n);
	my $para=$self->{'setting:cufflinks'};
	$cufflinks_cmd .= qq(export para="$para"\n);
	my $workdir=$self->{"-workdir"};
	$cufflinks_cmd .= qq(export workdir=$workdir\n);
	$cufflinks_cmd .= qq(cd \${workdir}\n);
	
	$cufflinks_cmd .= "[[ -d cufflinks ]] || mkdir cufflinks\n" if (!-d qq($self->{"-workdir"}/cufflinks));
	$cufflinks_cmd .= "cd cufflinks\n";
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		if (exists $attrs{'-ref'} && ($attrs{'-ref'} eq "yes" || $attrs{'-ref'}==1)) {
			$cufflinks_cmd .= qq(\${cufflinks} \${para} -g \$refGene -b \$reference -o $lib\_cufflinks $self->{$lib}{"tophatbam"}\n);
		} else {
			$cufflinks_cmd .= qq(\${cufflinks} \${para} -b \$reference -o $lib\_cufflinks $self->{$lib}{"tophatbam"}\n);
		}
		$cufflinks_cmd .= qq(\${gffread} -w transcripts.fa -g \$reference transcripts.gtf\n);
		$self->{$lib}{'cufflinks'}=$self->{'-workdir'}."/cufflinks/$lib\_cufflinks";
		$self->{$lib}{'transcript-gtf'}=$self->{'-workdir'}."/cufflinks/$lib\_cufflinks/transcripts.gtf";
		push @{$self->{'transcript-gtf'}},$self->{'-workdir'}."/cufflinks/$lib\_cufflinks/transcripts.gtf";
	}
	$cufflinks_cmd .= "cd ..\n";
	return ($cufflinks_cmd);
}

sub runCuffMerge($) {
	my $self=shift;
	my $ref=shift;
	my $gene=shift;
	$ref ||= 'ref';
	$gene ||= 'refGene';
	if (!exists $self->{"software:cuffmerge"} || !defined $self->{"software:cuffmerge"}) {
		return "";
	}
	#cuffmerge -s ../05.db/w14_v7.2.fasta -p 4 assembly_GTF_list.txt
	my $reference=GATE::Error::checkPath($self->{"database:$ref"});
	my $para=$self->{'setting:cuffmerge'};
	my $cuffmerge="";
	if (exists $self->{'software:cuffmerge'}) {
		$cuffmerge=GATE::Error::checkPath($self->{'software:cuffmerge'});
	} else {
		$cuffmerge=GATE::Error::checkPath($self->{'software:cufflinks'});
		$cuffmerge=~s/cufflinks$/cuffmerge/;
	}
	my $cuffmerge_cmd = qq(echo `date`; echo "run Cuffmerge"\n);
	$cuffmerge_cmd .= qq(cd $self->{"-workdir"}\n);
	$cuffmerge_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$cuffmerge_cmd .= qq(export REFERENCE=$reference\n);
	$cuffmerge_cmd .= "[[ -d cuffmerge ]] || mkdir cuffmerge\n" if (!-d qq($self->{"-workdir"}/cuffmerge));
	$cuffmerge_cmd .= "cd cuffmerge\n";
	if (exists $self->{'database:$gene'}) {
		my $refGene=GATE::Error::checkPath($self->{'database:$gene'});
		my $all_transcripts=join "\\n",@{$self->{'transcript-gtf'}};
		$cuffmerge_cmd.=qq(perl -e 'print "$all_transcripts\\n"' > GTF_list.txt\n);
		$cuffmerge_cmd.="$cuffmerge $para -g $refGene -s \$REFERENCE -p 10 GTF_list.txt\n";
	} else {
		my $all_transcripts=join "\\n",@{$self->{'transcript-gtf'}};
		$cuffmerge_cmd.=qq(perl -e 'print "$all_transcripts\\n"' > GTF_list.txt\n);
		$cuffmerge_cmd.=qq($cuffmerge $para -s \$REFERENCE -p $self->{"setting:multithreads"} GTF_list.txt\n);
	}
	$self->{'merged-gtf'}=$self->{'-workdir'}."/cuffmerge/merged_asm/merged.gtf";
	$self->{'merged-transcripts-gtf'}=$self->{'-workdir'}."/cuffmerge/merged_asm/transcripts.gtf";
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		$cuffmerge_cmd .= qq(echo `date`; echo "$lib"\n);
		if (exists $self->{$lib}{'transcript-gtf'}) {
			$cuffmerge_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/cuffmerge/$lib));
			$cuffmerge_cmd .= "cd $lib\n";
			if (exists $self->{'database:$gene'}) {
				my $refGene=GATE::Error::checkPath($self->{'database:$gene'});
				$cuffmerge_cmd.=qq(echo "$self->{$lib}{'transcript-gtf'}" > $lib\_GTF_list.txt\n);
				$cuffmerge_cmd.="$cuffmerge $para -g $refGene -s \$REFERENCE -p 10 $lib\_GTF_list.txt\n";
			} else {
				$cuffmerge_cmd.=qq(echo "$self->{$lib}{'transcript-gtf'}" > $lib\_GTF_list.txt\n);
				$cuffmerge_cmd.="$cuffmerge $para -s \$REFERENCE -p 10 $lib\_GTF_list.txt\n";
			}
			$self->{$lib}{'merged-gtf'}=$self->{'-workdir'}."/cuffmerge/merged_asm/merged.gtf";
			$self->{$lib}{'merged-transcripts-gtf'}=$self->{'-workdir'}."/cuffmerge/merged_asm/transcripts.gtf";
			$cuffmerge_cmd .= "cd ..\n";
		}
	}
	$cuffmerge_cmd .= "cd ..\n";
	return ($cuffmerge_cmd);
}

sub runCuffCompare($) {
#cuffcompare -i input_gtf_list -r ../reference/genes.gtf -o compare_out 2> cuffcompare.log &
#cuffcompare -s ../05.db/w14_v7.2.fasta -p 4 assembly_GTF_list.txt
	my $self=shift;

	if (!exists $self->{"software:cuffcompare"} || !defined $self->{"software:cuffcompare"}) {
		return "";
	}
	my $cuffcompare_cmd = qq(echo `date`; echo "run Cuffcompare"\n);;
	$cuffcompare_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	my $reference=GATE::Error::checkPath($self->{'database:ref'});
	$cuffcompare_cmd .= qq(export REFERENCE\="$reference"\n);
	if (exists $self->{'database:$gene'}) {
		my $refGene=$self->{'database:$gene'};
		$cuffcompare_cmd.=qq(export refGene="$refGene"\n);
	}
	my $para=$self->{'setting:cuffcompare'};
	$cuffcompare_cmd .= qq(export para="$para"\n);
	my $cuffcompare="";
	if (exists $self->{'software:cuffcompare'}) {
		$cuffcompare=GATE::Error::checkPath($self->{'software:cuffcompare'});
	} else {
		$cuffcompare=GATE::Error::checkPath($self->{'software:cufflinks'});
		$cuffcompare=~s/cufflinks$/cuffcompare/;
	}
	$cuffcompare_cmd=qq(export cuffcompare=$cuffcompare\n);
	my $workdir=$self->{"-workdir"};
	$cuffcompare_cmd=qq(export workdir=$workdir\n);
	$cuffcompare_cmd=qq(cd \${workdir}\n);
	$cuffcompare_cmd .= "[[ -d cuffcompare ]] || mkdir cuffcompare\n" if (!-d qq($self->{"-workdir"}/cuffcompare));
	$cuffcompare_cmd .= "cd cuffcompare\n";

	if (exists $self->{'database:refGene'}) {
		my $refGene=$self->{'database:refGene'};
		$cuffcompare_cmd.=qq(export refGene="$refGene"\n);
		my $all_transcripts=join "\\n",@{$self->{'transcript-gtf'}};
		$cuffcompare_cmd.=qq(perl -e 'print "$all_transcripts\\n"' > GTF_list.txt\n);
		$cuffcompare_cmd.="\${cuffcompare} \${para} -r \${refGene} -i GTF_list.txt -o cuffcompare\n";
	} else {
		my $all_transcripts=join "\\n",@{$self->{'transcript-gtf'}};
		$cuffcompare_cmd.=qq(perl -e 'print "$all_transcripts\\n"' > GTF_list.txt\n);
		$cuffcompare_cmd.="\${cuffcompare} \${para} -i GTF_list.txt -o cuffcompare\n";
	}
	$self->{'compare-transcript-gtf'}=$self->{'-workdir'}."/cuffcompare/cuffcompare.combined.gtf";
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		$cuffcompare_cmd .= qq(echo `date`; echo "$lib"\n);
		if (exists $self->{$lib}{'transcript-gtf'}) {
			$cuffcompare_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/cuffcompare/$lib));
			$cuffcompare_cmd .= "cd $lib\n";
			if (exists $self->{'database:$gene'}) {
				$cuffcompare_cmd.=qq(echo "$self->{$lib}{'transcript-gtf'}" > $lib\_GTF_list.txt\n);
				$cuffcompare_cmd.="\${cuffcompare} \${para} -r \${refGene} -i GTF_list.txt -o $lib\n";
			} else {
				$cuffcompare_cmd.=qq(echo "$self->{$lib}{'transcript-gtf'}" > $lib\_GTF_list.txt\n);
				$cuffcompare_cmd.="\${cuffcompare} \${para} -i GTF_list.txt -o $lib\n";
			}
			$self->{$lib}{'compare-transcript-gtf'}=$self->{'-workdir'}."/cuffcompare/$lib/$lib.combined.gtf";
			$cuffcompare_cmd .= "cd ..\n";
		}
	}
	$cuffcompare_cmd .= "cd ..\n";
	return ($cuffcompare_cmd);
}

sub runTrinity($) {
#Trinity.pl --seqType fq --JM 100G --left reads_1.fq  --right reads_2.fq --CPU 6
#TRINITY_RNASEQ_ROOT/util/alignReads.pl --left left.fq --right right.fq --seqType fq --target Trinity.fasta --aligner bowtie
#TRINITY_RNASEQ_ROOT/util/RSEM_util/run_RSEM.pl --transcripts Trinity.fasta --name_sorted_bam bowtie_out.nameSorted.sam.+.sam.PropMapPairsForRSEM.bam --paired
	my $self=shift;
	my $trinity=GATE::Error::checkPath($self->{"software:trinity"});
	my $trinity_cmd = qq(echo `date`; echo "run Trinity"\n);
	$trinity_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	if (exists $self->{'setting:TRINITY_RNASEQ_ROOT'}) {
		$trinity_cmd .= qq(export TRINITY_RNASEQ_ROOT="$self->{'setting:TRINITY_RNASEQ_ROOT'}"\n);
	} else {
		my $trinity_root=$1 if ($trinity=~/(\S+)\/[^\/]$/);
		$trinity_cmd .= qq(export TRINITY_RNASEQ_ROOT="$trinity_root"\n);
	}
	my $para=(exists $self->{"setting:trinity"})?$self->{"setting:trinity"}:qq(--seqType fq --JM 100G --CPU $self->{"setting:multithreads"});
	$trinity_cmd .= qq(cd $self->{"-workdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		$trinity_cmd .= qq(echo `date`; echo "$lib"\n);
		$trinity_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$lib));
		$trinity_cmd .= qq(cd $lib\n);
		$trinity_cmd .= qq([[ -d trinity_asm ]] || mkdir trinity_asm\n) unless (-d qq($self->{"-workdir"}/$lib/trinity_asm));
		$trinity_cmd .= qq(cd trinity_asm\n);
		my %read=getlibInput($self->{"LIB"}{$lib});
		if (exists $read{1} && exists $read{2}) {
			if (@{$read{1}}>1 && @{$read{2}}>1){
				my ($reads1,$reads2);
				if (@{$read{1}}==@{$read{2}}) {
					$reads1=join " ",@{$read{1}};
					$reads2=join " ",@{$read{2}};
				} else {
					if (@{$read{1}}<@{$read{2}}) {
						$reads1=join " ",@{$read{1}}[0..(@{$read{1}}-1)];
						$reads2=join " ",@{$read{2}}[0..(@{$read{1}}-1)];
						push @{$read{0}},@{$read{1}}[@{$read{1}}..(@{$read{2}}-1)];
					} else {
						$reads1=join " ",@{$read{1}}[0..(@{$read{2}}-1)];
						$reads2=join " ",@{$read{2}}[0..(@{$read{2}}-1)];
						push @{$read{0}},@{$read{2}}[@{$read{2}}..(@{$read{1}}-1)];
					}
				}
				my $single;
				if (exists $read{0}) {
					$single=join " ",@{$read{0}}[0..(@{$read{0}}-1)];
				}
				if (defined $reads1 && defined $reads2)
				{
					$trinity_cmd .= qq(cat $reads1 > reads_1.fq && cat $reads2 > reads_2.fq\n);
					if (defined $single)
					{
						if (@{$read{0}}>1) {
							$trinity_cmd .= qq(cat $single > single_reads.fq\n);
							$single = "single_reads.fq";
						}
					}
					$trinity_cmd .= qq($trinity --left reads_1.fq --right reads_2.fq $para);
					$trinity_cmd .= qq( --single $single);
					$trinity_cmd .= qq(\n);
					$trinity_cmd .= qq(\${TRINITY_RNASEQ_ROOT}/util/alignReads.pl --left reads_1.fq --right reads_2.fq --seqType fq --target Trinity.fasta --aligner bowtie);
					$trinity_cmd .= qq( --single $single);
					$trinity_cmd .= qq(\n);
					$trinity_cmd .= qq(\${TRINITY_RNASEQ_ROOT}/util/RSEM_util/run_RSEM.pl --transcripts Trinity.fasta --name_sorted_bam bowtie_out/bowtie_out.nameSorted.bam --paired\n);
				}
			} else {
				my $single;
				if (exists $read{0}) {
					$single=join " ",@{$read{0}}[0..(@{$read{0}}-1)];
				}
				if (defined $single)
				{
					if (@{$read{0}}>1) {
						$trinity_cmd .= qq(cat $single > single_reads.fq\n);
						$single = "single_reads.fq";
					}
				}
				$trinity_cmd .= qq($trinity --left ${$read{1}}[0] --right ${$read{2}}[0] $para);
				$trinity_cmd .= qq( --single $single);
				$trinity_cmd .= qq(\n);
				$trinity_cmd .= qq(\${TRINITY_RNASEQ_ROOT}/util/alignReads.pl --left reads_1.fq --right reads_2.fq --seqType fq --target Trinity.fasta --aligner bowtie);
				$trinity_cmd .= qq( --single $single);
				$trinity_cmd .= qq(\n);
				$trinity_cmd .= qq(\${TRINITY_RNASEQ_ROOT}/util/RSEM_util/run_RSEM.pl --transcripts Trinity.fasta --name_sorted_bam bowtie_out/bowtie_out.nameSorted.bam --paired\n);
			}
		} elsif (exists $read{1}) {
			$trinity_cmd .= qq($trinity --single ${$read{1}}[0] $para\n);
			$trinity_cmd .= qq(\${TRINITY_RNASEQ_ROOT}/util/alignReads.pl --single ${$read{1}}[0] --seqType fq --target Trinity.fasta --aligner bowtie\n);
			$trinity_cmd .= qq(\${TRINITY_RNASEQ_ROOT}/util/RSEM_util/run_RSEM.pl --transcripts Trinity.fasta --name_sorted_bam bowtie_out/bowtie_out.nameSorted.bam --SS_lib_type F --thread_count 8\n);
		} elsif (exists $read{2}) {
			$trinity_cmd .= qq($trinity --single ${$read{1}}[2] $para\n);
			$trinity_cmd .= qq(\${TRINITY_RNASEQ_ROOT}/util/alignReads.pl --single ${$read{2}}[0] --seqType fq --target Trinity.fasta --aligner bowtie\n);
			$trinity_cmd .= qq(\${TRINITY_RNASEQ_ROOT}/util/RSEM_util/run_RSEM.pl --transcripts Trinity.fasta --name_sorted_bam bowtie_out/bowtie_out.nameSorted.bam--SS_lib_type F --thread_count 8\n);
		} elsif (exists $read{0}) {
			$trinity_cmd .= qq($trinity --single ${$read{0}}[0] $para\n);
			$trinity_cmd .= qq(\${TRINITY_RNASEQ_ROOT}/util/alignReads.pl --single ${$read{0}}[0] --seqType fq --target Trinity.fasta --aligner bowtie\n);
			$trinity_cmd .= qq(\${TRINITY_RNASEQ_ROOT}/util/RSEM_util/run_RSEM.pl --transcripts Trinity.fasta --name_sorted_bam bowtie_out/bowtie_out.nameSorted.bam --SS_lib_type F --thread_count 8\n);
		}
		push @{$self->{"LIB"}{$lib}{'INPUT'}{'denovo'}},qq($self->{"-workdir"}/trinity_asm/$lib/Trinity.fasta);
	}
	return $trinity_cmd;
}

sub runVelvetOases ($) {
	my $self=shift;
	my $velveth=GATE::Error::checkPath($self->{"software:velveth"});
	my $velvet_cmd = qq(echo `date`; echo "run Velvet-Oases"\n);
	$velvet_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$velvet_cmd .= qq(export velvet="$velveth"\n);
	my $velvetg="";
	if (exists $self->{"software:velvetg"})
	{
		$velvetg=GATE::Error::checkPath($self->{"software:velvetg"});
		$velvet_cmd .= qq(export velveg="$velvetg"\n);
	}
	else
	{
		$velvetg=$self->{"software:velveth"};
		$velvetg=~s/velveth$/velvetg/;
		$velvet_cmd .= qq(export velveg="$velvetg"\n);
	}
	my $oases=GATE::Error::checkPath($self->{"software:oascs"}) if (exists $self->{"software:oascs"});
	$velvet_cmd .= qq(export oases="$oases"\n) if (defined $oases);
	my $velvethpara = $self->{'setting:velveth'};
	$velvet_cmd .= qq(export velvethpara="$velvethpara"\n);
	my $velvetgpara = $self->{'setting:velvetg'};
	$velvet_cmd .= qq(export velvetgpara="$velvetgpara"\n);
	my $oasespara = $self->{'setting:oases'} if (exists $self->{'setting:oases'});
	$velvet_cmd .= qq(export oasespara="$oasespara"\n) if (defined $oasespara);
	my @libraries=sort keys %{$self->{'LIB'}};
	$velvet_cmd .= qq(cd $self->{"-workdir"}\n);
	$velvet_cmd .= qq([[ -d velvet ]] || mkdir velvet\n);
	foreach my $lib(@libraries) {
		$velvet_cmd .= qq(echo `date`; echo "$lib"\n);
		$velvet_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$lib));
		$velvet_cmd .= qq(cd $lib\n);
		my %read=getlibInput($self->{"LIB"}{$lib});
		my $input="";
		if (exists $read{0}) {
			for (my $i=0;$i<@{$read{0}};$i++) {
				$input.="-short --".$self->checkFileFormat(${$read{0}}[$i])." ${$read{0}}[$i]";
			}
		}
		if (exists $read{1} && exists $read{2} && @{$read{2}}==@{$read{1}}) {
			for (my $i=0;$i<@{$read{2}};$i++) {
				my $reads1=${$read{1}}[$i];
				my $reads2=${$read{2}}[$i];
				my $format=$self->checkFileFormat($reads1);
				my $reads=shuffleSequences($reads1,$reads2,$format);
				$input.="-shortPaired --".$self->checkFileFormat($reads)." $reads";
			}
		}
		$velvet_cmd .= qq(\${velveth} $lib \${velvethpara} $input\n);
		$velvet_cmd .= qq(\${velvetg} $lib \${velveghpara});
		$velvet_cmd .= qq( -read_trkg yes) if ($velvetgpara!~/read_trkg/);
		$velvet_cmd .= "\n";
		if (defined $oases)
		{
			$velvet_cmd .= qq(\${oases} $lib \${oasespara}\n) ;
			push @{$self->{"LIB"}{$lib}{'INPUT'}{'denovo'}},qq($self->{"-workdir"}/velvet/$lib/transcripts.fa);
		} else {
			push @{$self->{"LIB"}{$lib}{'INPUT'}{'denovo'}},qq($self->{"-workdir"}/velvet/$lib/contigs.fa);
		}
	}
	return $velvet_cmd;
}

sub shuffleSequences ($) {
	
}

sub runABySS ($) {
	my $self=shift;
}

sub runTransABySS ($) {
	my $self=shift;
}

sub runSOAPdenovo ($) {
#SOAPdenovo-63mer all -s WGS_soapdenovo.config -K 25 -F -p 32 -a 92 -m 51 -M 1 -E -k 31 -F -L 100 -V -o WGS_25mer 1> 25mer.log 2> 25mer.log
	my $self=shift;
	if (!exists $self->{"software:soapdenovo"})
	{
		return "";
	}
	my $soapdenovo=GATE::Error::checkPath($self->{"software:soapdenovo"}) if (exists $self->{"software:soapdenovo"});
	my $soapdenovo_cmd = qq(echo `date`; echo "run SOAPdenvo|SOAPdenovo-Trans"\n);
	$soapdenovo_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$soapdenovo_cmd .= qq(export soapdenvo="$soapdenovo"\n);
	my $config=GATE::Error::checkPath($self->{"setting:soapdenovo_config"}) if (exists $self->{"setting:soapdenovo_config"});
	my $para=$self->{"setting:soapdenovo"};
	$soapdenovo_cmd .= qq(export soapdenvo_para="$para"\n);
	$config=$1 if ($para=~/\-s\s+(\S+)/);
	$soapdenovo_cmd .= qq(cd $self->{"-workdir"}\n);
	$soapdenovo_cmd .= qq([[ -d soapdenovo] || mkdir soapdenvo\n);
	if (defined $config) {
		$soapdenovo_cmd .= qq(\${soapdenvo} \${soapdenvo_para}});
		$soapdenovo_cmd .= qq( -s $config) if ($para !~/\-s/);
		$soapdenovo_cmd .= "\n";
		my $prefix = "";
		if ($para=~/\-o\s+(\S+)/) {
			$prefix = $1;
		}
		else {
			$prefix = "soapdenovo";
		}
		$self->{"denovo_genomics"}=qq($self->{"-workdir"}/soapdenvo/$prefix.scafSeq);
	} else {
		print STDERR "no SOAPdenovo.config file found!\n";
		my @libraries=sort keys %{$self->{'LIB'}};
		foreach my $lib(@libraries) {
			my $makeconfig=$self->make_config("program","SOAPdenovo","outdir",qq($self->{"-workdir"}/$soapdenovo),"lib",$lib);
			$soapdenovo_cmd .= qq(\${soapdenvo} \${soapdenvo_para}} -s $makeconfig -o $lib\n);
			$self->{"denovo_genomics"}=qq($self->{"-workdir"}/soapdenvo/$lib.scafSeq);
		}
	}
	return $soapdenovo_cmd;
}

sub runSOAPdenvoTrans ($) {
	
}

sub make_config {
	my ($self,%attrs)=@_;
	my $program = $attrs{"program"};
	my $outdir  = $attrs{"outdir"};
	my $lib     = $attrs{"lib"};
	my $config  = "$outdir/$lib.$program.config";
	if ($program =~ /SOAPdenovo/i && exists $self->{"LIB"}{$lib}) {
		if (-f $config) {
			die "$config is existing!\n";
		}
		open(OUT,">$config") or die $!;
		my %read=getlibInput($self->{"LIB"}{$lib});
		if (exists $read{1} && exists $read{2}) {
			print OUT "[LIB]\n";
			for (my $i=0;$i<@{$read{2}};$i++) {
				my $type=(checkFileFormat(${$read{1}}[$i])=~/fastq/i)?"fq":"fa";
				my $R1=$type."1";
				my $ins=(exists $self->{$lib}{$R1}{$i}{'PI'}) ? $self->{$lib}{$R1}{$i}{'PI'} : 200;
				print OUT "avg_ins=$ins\n";
				my $reverse=($ins>600)?1:0;
				print OUT "reverse_seq=$reverse\n";
				my $rank=1;
				if ($ins<=500) {
					$rank=2;
				} elsif ($ins<1000) {
					$rank=3;
				}
				else {
					$rank=int($ins/1000)+3;
				}
				print OUT "rank=$rank\n";
				print OUT "pair_num_cutoff=3\n";
				print OUT "asm_flag=3\n";
				my $head=($type=~/fq/) ? "q" : "f";
				print OUT $head.qq(1=${$read{1}}[$i]\n);
				print OUT $head.qq(2=${$read{2}}[$i]\n);
			}
		}
		if (exists $read{0}) {
			print OUT "[LIB]\n";
			for (my $i=0;$i<@{$read{0}};$i++) {
				print OUT "rank=1\nasm_flags=1\n";
				my $head=(checkFileFormat(${$read{0}}[$i])=~/fastq/i)?"q":"f";
				print OUT qq($head=${$read{0}}[$i]\n);
			}
		}
		
		close OUT;
		return $config;
	}
}


sub runGapCloser ($) {
	
}

sub runGapFiller ($) {
	
}

sub runPhrap ($) {
#phrap seq.fas -new_ace -revise_greedy -shatter_greedy -forcelevel 0 -repeat_stringency 0.95 > phrap.out
	
}

sub runPhusion ($) {
	
}

sub runPAGIT ($) {
	
}

sub runREAPR ($) {
	
}

sub runCAP3 ($) {
	
}

sub runPCAP ($) {
	
}

sub runNewbler ($) {
#/opt/454/bin/newAssembly BAC1
#/opt/454/bin/addRun BAC1 ~/01.data/W14/454/BAC/5_BAC_454/BAC_1.sff 
#/opt/454/bin/runProject -cpu 8 BAC1
	my $self = shift;
	if (!exists $self->{"setting::Newbler"} || !-d $self->{"setting::Newbler"}) {
		return "";
	}
	my $newbler_cmd = qq(echo `date`; echo "run Newbler"\n);
	$newbler_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$newbler_cmd .= qq(export NEWBLER="$self->{"setting::Newbler"}"\n);
	my $para = $self->{'setting:Newbler'} if (exists $self->{'setting:Newbler'});
	$newbler_cmd .= qq(export para=$para) if (defined $para);
	my $workdir .= GATE::Error::checkPath($self->{"-workdir"});
	$newbler_cmd .= qq(export workdir=$workdir\n);
	$newbler_cmd .= qq(cd \${workdir}\n);
	$newbler_cmd .= "[[ -d newbler }] || mkdir newbler\n" if (!-d qq($self->{"-workdir"}/newbler));
	$newbler_cmd .= qq(cd newbler\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries)
	{
		my %read=getlibInput($self->{"LIB"}{$lib});
		$newbler_cmd .= qq(\${NEWBLER}/newAssembly $lib\n);
		if (exists $read{1} && exists $read{2}) {
			for (my $i=0;$i<@{$read{2}};$i++) {
				$newbler_cmd .= qq(\${NEWBLER}/addRun $lib -p ${$read{1}}[$i]\n);
				$newbler_cmd .= qq(\${NEWBLER}/addRun $lib -p ${$read{2}}[$i]\n);
			}
		}
		if (exists $read{0}) {
			for (my $i=0;$i<@{$read{0}};$i++) {
				$newbler_cmd .= qq(\${NEWBLER}/addRun $lib ${$read{0}}[$i]\n);
			}
		}
		$newbler_cmd .= qq(\${NEWBLER}/runProject \${para}\n);
	}
	return $newbler_cmd;
}

sub runMSR_CA($) {
#==>sr_config.txt<==
#PATHS
#JELLYFISH_PATH=/usr/local/genome/MSR-CA-1.6.1/bin
#SR_PATH=/usr/local/genome/MSR-CA-1.6.1/bin
#CA_PATH=/usr/local/genome/MSR-CA-1.6.1/CA/Linux-amd64/bin
#END
#
#DATA
#PE= pe 400 40  400.20110831.s_2_uniq1.fastq 400.20110831.s_2_uniq2.fastq
#JUMP= mp 1000 100 1.0k.20110518.s_5_uniq2.fastq 1.0k.20110518.s_5_uniq1.fastq
#END
#
#PARAMETERS
#LIMIT_JUMP_COVERAGE = 60
#CA_PARAMETERS = ovlMerSize=30 cgwErrorRate=0.15 ovlMemory=4GB
#KMER_COUNT_THRESHOLD = 1
#EXTEND_JUMP_READS=1
#NUM_THREADS= 16
#JF_SIZE=500000000
#END
	my $self = shift;
}

#########################################################
#                                                       #
#              Quantitative analysis                    #
#                                                       #
#########################################################

sub runCuffdiff($) {
	my $self=shift;
	my $ref = shift;
	my $gene = shift;
	$ref ||= 'ref';
	$gene ||= 'refGene';
#cuffdiff -p 12 -o UF_BF_ST_diff -b ../reference/hg19.fa -u ../reference/genes.gtf
#-L UF,BF,ST UF1_tophat/accepted_hits.bam,UF2_tophat/accepted_hits.bam,UF3_tophat/accepted_hits.bam
#BF1_tophat/accepted_hits.bam,BF2_tophat/accepted_hits.bam,BF3_tophat/accepted_hits.bam
#ST1_tophat/accepted_hits.bam,ST2_tophat/accepted_hits.bam,ST3_tophat/accepted_hits.bam
	my $cuffdiff_cmd = qq(echo `date`; echo "run Cuffdiff"\n);
	$cuffdiff_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	my $reference=GATE::Error::checkPath($self->{"database:$ref"});
	$cuffdiff_cmd .= qq(export reference=$reference\n);
	my $refGene=(exists $self->{"database:$gene"})?GATE::Error::checkPath($self->{"database:$gene"}):$self->{$gene};
	$cuffdiff_cmd .= qq(export refGene=$refGene\n);
	my $cuffdiff="";
	if (exists $self->{'software:cuffdiff'}) {
		$cuffdiff=GATE::Error::checkPath($self->{'software:cuffdiff'});
	} else {
		$cuffdiff=GATE::Error::checkPath($self->{'software:cufflinks'});
		$cuffdiff=~s/cufflinks$/cuffdiff/;
	}
	$cuffdiff_cmd .= qq(export cuffdiff=$cuffdiff\n);
	my $para = $self->{'setting:cuffdiff'};
	$cuffdiff_cmd .= qq(export para="$para"\n);
	my $workdir .= GATE::Error::checkPath($self->{"-workdir"});
	$cuffdiff_cmd .= qq(export workdir=$workdir\n);
	$cuffdiff_cmd .= qq(cd \${workdir}\n);
	$cuffdiff_cmd .= "[[ -d cuffdiff ]] || mkdir cuffdiff\n" if (!-d qq($self->{"-workdir"}/cuffdiff));
	$cuffdiff_cmd .= "cd cuffdiff\n";
	my @libraries=sort keys %{$self->{'LIB'}};
	my ($label,$bam)=("","");
	foreach my $sm(@libraries) {
		$label.="$sm,";
		$bam.=$self->{$sm}{"tophatbam"}." ";
	}
	$label=~s/\,$//;
	$bam=~s/[\,\s]+$//;
	$cuffdiff_cmd .= qq(\${cuffdiff} \${para} -o expdiff -b \$reference -L $label \$refGene $bam\n);
	$self->{expdiff}=$self->{'-workdir'}."/cuffdiff/expdiff";
	$cuffdiff_cmd .= "cd ..\n";
	return ($cuffdiff_cmd);
}

sub runcummeRbund ($) {
	my $self = shift;
}


#########################################################
#                                                       #
#              Alternative Splicing                     #
#                                                       #
#########################################################

sub runAlternativeSplicing($) {
	my $self = shift;
	if (!exists $self->{"software:AlternativeSplicing"} || !-f $self->{"software:AlternativeSplicing"})
	{
		return "";
	}
	my $samtools = GATE::Error::checkPath($self->{"software:samtools"});
	my $SamToFastq=GATE::Error::checkPath($self->{"software:SamToFastq"});
	my $AlternativeSplicing=GATE::Error::checkPath($self->{"software:AlternativeSplicing"});
	my $bwa=GATE::Error::checkPath($self->{"software:bwa"});
	my $getAlnGene = GATE::Error::checkPath($self->{"software:getAlnGene"});
	my $reference=GATE::Error::checkPath($self->{"database:ref"});
	my $refGene=GATE::Error::checkPath($self->{"database:refGene"});
	
	my @libraries=sort keys %{$self->{'LIB'}};
	my $runAS_cmd = qq(echo `date`; echo "run AlternativeSplicing"\n);
	$runAS_cmd .= qq(cd $self->{"-workdir"}\n);
	$runAS_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$runAS_cmd .= qq(export REFERENCE=$reference\n);
	$runAS_cmd .= qq(export refGene=$refGene\n);
	$runAS_cmd .= qq(mkdir $self->{"setting:as_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:as_outdir"}));
	$runAS_cmd .= "cd as\n";
	foreach my $lib(@libraries) {
		$runAS_cmd .= qq(echo `date`; echo "$lib"\n);
		$runAS_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) if (!-d qq($self->{"-workdir"}/as/$lib));
		$runAS_cmd .= "\ncd $lib\n";
		my $bam="";
		if (!exists $self->{$lib}{"ref-bam"}) {
			if (exists $self->{$lib}{"ref-bwabam"} && @{$self->{$lib}{"ref-bwabam"}}>1) {
				my $merge_bam=join " ",@{$self->{$lib}{"ref-bwabam"}};
				$runAS_cmd .=  qq($samtools view -H ${$self->{$lib}{"ref-bwabam"}}[0] > $lib.inh.sam\n);
				$runAS_cmd .=  "$samtools merge -nr -h $lib.inh.sam $lib.merge.bam $merge_bam\n";
				$runAS_cmd .= "$samtools sort $lib.merge.bam $lib.merge.sort\n";
				$runAS_cmd .= "$samtools index $lib.merge.bam $lib.merge.sort.bam\n";
				$bam="$lib.merge.sort.bam";
				$self->{$lib}{"ref-bam"}=$self->{'-workdir'}."/".$self->{"setting:as_outdir"}."/$lib/$lib.merge.sort.bam";
			} else {
				$bam=${@{$self->{$lib}{"ref-bwabam"}}}[0];
				$self->{$lib}{"ref-bam"}=${@{$self->{$lib}{"ref-bwabam"}}}[0];
			}
		} else {
			$bam=$self->{$lib}{"ref-bam"};
		}
		$runAS_cmd .= "## get no Mapped reads\n";
		my $noMappedbam=(split /\//,$bam)[-1];
		
		$noMappedbam=~s/\.bam$/\.noMappedToGenome\.bam/;
		$runAS_cmd .= "$samtools view -f 4 -bh $bam > $noMappedbam\n";
		my $noMappedfq=$noMappedbam;
		$noMappedfq=~s/\.bam$/\.fq/;
		$runAS_cmd .= "## convert to fastq\n";
		$runAS_cmd .= "$SamToFastq VALIDATION_STRINGENCY=SILENT I=$noMappedbam F=$noMappedfq\n";
		$runAS_cmd .= "## check aln refgene\n";
		$runAS_cmd .= "$getAlnGene $bam \${refGene} > $lib\_aln\_gene.gtf\n";
		my $alnGene = "$lib\_aln\_gene.gtf";
		$runAS_cmd .= "## get splicing junction reads\n";
		$runAS_cmd .= "$AlternativeSplicing -sp -i $bam -r \$REFERENCE -g $alnGene -p 10\n";
		$runAS_cmd .= "$bwa index -a bwtsw SPresult/spliceJunctions.fa\n";
		$runAS_cmd .= "## mapping with BWA\n";
		my $noMappedsai=$noMappedfq;
		$noMappedsai=~s/fq$/sai/;
		$runAS_cmd .= "$bwa aln -t 12 -n 2 SPresult/spliceJunctions.fa $noMappedfq > $noMappedsai\n";
		$runAS_cmd .= "$bwa samse SPresult/spliceJunctions.fa $noMappedsai $noMappedfq > mappingToJunctions.sam\n";
		$runAS_cmd .= "## Annotate the patterns of alternative splicing patterns\n";
		$runAS_cmd .= "$AlternativeSplicing -as -i mappingToJunctions.sam -g $alnGene -s SPresult/spliceSites.txt -t SPresult/transShortRun.txt\n";
		$runAS_cmd .= "$samtools view -b -h -S ASresult/JunctionRecord.sam > ASresult/JunctionRecord.bam\n";
		$runAS_cmd .= "$samtools sort ASresult/JunctionRecord.bam ASresult/JunctionRecord.sort\n";
		$runAS_cmd .= "$samtools index ASresult/JunctionRecord.sort.bam\n";
		$runAS_cmd .= "cd ..\n";
		$self->{$lib}{"AS_result-gff"}=$self->{'-workdir'}."/".$self->{"setting:as_outdir"}."/$lib/ASresult/AS_result.gff";
		$self->{$lib}{"Junctions-bed"}=$self->{'-workdir'}."/".$self->{"setting:as_outdir"}."/$lib/ASresult/Junctions.bed";
	}
	$runAS_cmd .= "cd ..\n";
	return ($runAS_cmd);
}

sub runASAP ($) {
	my $self = shift;
	my $samtools = GATE::Error::checkPath($self->{"software:samtools"});
	my $asap = GATE::Error::checkPath($self->{"software:asap"});
	my $para = (exists $self->{"setting:asap"}) ? $self->{"setting:asap"} : "";
	my $reference=GATE::Error::checkPath($self->{"database:ref"});
	my $refGene=GATE::Error::checkPath($self->{"database:refGene"});
	
	my @libraries=sort keys %{$self->{'LIB'}};
	my $runAS_cmd = qq(echo `date`; echo "run ASAP"\n);
	$runAS_cmd .= qq(cd $self->{"-workdir"}\n);
	$runAS_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$runAS_cmd .= qq(export REFERENCE=$reference\n);
	$runAS_cmd .= qq(export refGene=$refGene\n);
	$runAS_cmd .= qq(mkdir $self->{"setting:as_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:as_outdir"}));
	$runAS_cmd .= "cd as\n";
	foreach my $lib(@libraries) {
		$runAS_cmd .= qq(echo `date`; echo "$lib"\n);
		$runAS_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) if (!-d qq($self->{"-workdir"}/as/$lib));
		$runAS_cmd .= "cd $lib\n";
		if (exists $self->{$lib}{"tophatbam"}) {
			my $bam=$self->{$lib}{"tophatbam"};
			if (exists $self->{$lib}{'transcript-gtf'})
			{
				my $transcript=$self->{$lib}{'transcript-gtf'};
				$runAS_cmd .= "\${asap} -ref \${refGene} -tran $transcript -bam $bam $para > $lib\_AS.gff 2> $lib.error.log\n";
			}
		}
		$runAS_cmd .= "cd ..\n";
	}
	$runAS_cmd .= "cd ..\n";
	return ($runAS_cmd);
}

sub runGenePlot ($) {
	my $self = shift;
	my $ref = shift;
	my $gene = shift;
	$ref ||= 'ref';
	$gene ||= 'refGene';
#perl PlotGeneSplicing.pl --list gene_list.xls --ref refGene.gff \
#--tr TS1:S1_transctripts.gtf,TS2:S2_transcripts.gtf,TS3:S3_transcripts.gtf\
#--as AS1:AS1_result.gff,AS2:AS2_result.gff,AS3:AS3_result.gff \
#--bam BS1:BS1.bam,BS2:BS2.bam,BS3:BS3.bam \
#--jun JS1:JS1.bed,JS2:JS2.bed,JS3:JS3.bed
	if (!exists $self->{"software:GenePlot"} || !defined $self->{"software:GenePlot"}) {
		return "";
	}
	my $GenePlot = GATE::Error::checkPath($self->{"software:GenePlot"});
	my $plotas_cmd = qq(echo `date`; echo "run GenePlot"\n);
	$plotas_cmd .= qq(cd $self->{"-workdir"}\n);
	$plotas_cmd .= $self->runCuffcompare() if (!exists $self->{'compare-transcript-gtf'});
	#print qq(genelist=$self->{"database:list"}\n);
	my $genelist = GATE::Error::checkPath($self->{"database:genelist"});
	my $refGene = GATE::Error::checkPath($self->{"database:$gene"});
	
	$plotas_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$plotas_cmd .= "cd as\n";
	$plotas_cmd .= "mkdir AS_Figure\n" if (!-d qq($self->{"-workdir"}/as/AS_Fgiure/));
	$plotas_cmd .= "cd AS_Figure\n";
	$plotas_cmd .= qq($GenePlot --list $genelist --ref $refGene \\\n);

	my @libraries=sort keys %{$self->{'LIB'}};
	my ($tr,$as,$bam,$jun)=("","","","");
	foreach my $sm(@libraries) {
		$tr.=qq($sm:$self->{$sm}{'compare-transcript-gtf'},);
		$as.=qq($sm:$self->{$sm}{'AS_result-gff'},);
		$bam.=qq($sm:$self->{$sm}{"$ref-bam"},);
		$jun.=qq($sm:$self->{$sm}{'Junctions-bed'},);
	}
	$tr=~s/\,$//;
	$as=~s/\,$//;
	$bam=~s/\,$//;
	$jun=~s/\,$//;
	$plotas_cmd .= "-tr $tr \\\n" if ($tr ne "");
	$plotas_cmd .= "-as $as \\\n" if ($as ne "");
	$plotas_cmd .= "-bam $bam \\\n" if ($bam ne "");
	$plotas_cmd .= "-jun $jun \\\n" if ($jun ne "");
	$plotas_cmd .= "cd ../../\n";
	if (defined $GenePlot) {
		return ($plotas_cmd);
	} else {
		return "";
	}
}

sub runScripture ($) {
#
#Parameters 
# -alignment <Alignment file in BAM, SAM or Alignemnt format> 
# -maskFileDir <Mask File directory> 
# -out <Output file name>
# -sizeFile <Chromosome size file> 
# -chr <Chromsomosome to segment> 
# -chrSequence <Necessary to filter spliced reads by splice site information. Notice that this is only compatible with region files that contain regions of only one chromosome> 
# Optional arguments: 
# -windows <Comma separated list of windows to evaluate defaults to contiguous regions of coverage>
# -trim <Include this flag if trimming of the ends of windows based on read coverage  is desired this is expensive> 
# -alpha <Desired FDR>
#  -dontFilterCanonicalSplice
# -start <To segment only a subregion of the chromosome include its start> -end <To segment only a subregion of the chromosome include its end>
# -minSpliceSupport <Minimum count to support splice reads, default is 1> 
# -pairedEnd <Paired end alignment files> -strandSpecificReads <Strand specific alignment file> -scoreRegions <Full BED to score> -upWeightSplices -lambda <If a prior background expectation for number of reads per base exists> -exons <BED file of exons> -introns <Introns and counts>
#
#Task: AddPairs -  Uses a paired end alignment to tune graph 
#	-in <Graph in .dot format. Standard input is assumed> 
#	-pairedEnd <Paired end information (as in previous task), in single line BED format>
#	 -maskFileDir <Directory containing mask files for the genome> 
#	-chr <Chromosome (only a chromosome at a time is supported at this point)> 
#	-sizeFile <Chromosome size file> 
#	-out <Output file name>
#
#Task: fastScore -  Computes several expression related scores for a set of annotations using a the graph .dot file 
#	-in <chr.dot file> 
#	-annotations <BED file with annotation to score>
#	 -chr <chr e.g: chrZ>
#	 -alpha <optional>
#	 -out
#
#
#Task: score -  Computes several expression related scores for a set of annotations -in <Full BED file with annotations to score> 
#	-alignment <Alignment file in BAM, SAM or Alignemnt format> 
#	-sizeFile <Chromosome size file> 
#	-out <Output file name> 
#	 -maskFileDir <Mask File directory>
#
#Task: extractDot - Extracts a graph for the specified region 
#	-in <Dot file from a previous Scripture run 
#	-chr <Chromosome> 
#	-start <Start of region> 
#	-end<End of region> 
#	-out <output file>
#
#Task: getIdenticalGappedReadsTranscripts -  Report all transcripts that ALL their introns are spanned by identical gapped reads -in <Full BED file with annotations to score> 
#	-alignment <Alignment file in BAM, SAM or Alignemnt format> 
#	-sizeFile <Chromosome size file> 
#	-out <Output file name> 
#	 -maskFileDir <Mask File directory>
#
#Task: makePairedFile Makes a paired end alignment file from two sets of independtly aligned left and right ends, ideally the files should be name-sorted 
#	-pair1 <First pair alignments> 
#	-pair2 <Second pair alignments> 
#	-out <output consolidated paired end alignment> 
#	-sorted  <Include this flag if the data is already read name sorted, ideally both input files should be sorted by read name using unix sort for example> 
#	-sortBam <Include this flag if the input includes a single BAM file; a temporary SAM file sorted by name would be generated and than removed  > 
#	-sizeFile <Chromosome size file>
#
#Task: chipScan - Segment the genome assuming contiguous data. Similar to the default task but optimized for contiguous data. 
# -alignment <Alignment file in BAM, SAM or Alignemnt format> 
# -maskFileDir <Mask File directory> 
# -out <Output file name>
# -chr <Chromosome to segment>
# -sizeFile <Chromosome size file> 
# -windows <Comma separated list of windows to evaluate defaults to contiguous regions of coverage> 
# Optional arguments:
# -findMaxContiguous <Each significant window is trimmed by finding the contiguous sub region with coverage over a predefined threshold> -trim <Include this flag if trimming of the ends of windows based on read coverage  is desired this is expensive> 
# -alpha <Desired FDR>
#
#Task: trim -  Trims end of transcripts by removing all bases whose coverage is below the specified quantile of transcript expression -in <Full BED file with annotations to trim> 
#	-alignment <Alignment file in BAM, SAM or Alignemnt format> 
#	-sizeFile <Chromosome size file> 
#	-out <Output file name> 
#	 -maskFileDir <Mask File directory>
#	-quantile <Coverage quantile below which end bases should be trimmed>

## Make paired file task
# java -Xmx2000m jar scripture.jar -task makePairedFile <Mandatory parameters> <Options>
## Segmentation task
# java -Xmx2000m jar scripture.jar <Mandatory parameters> <optional parameter>
## Add pairs task
# java -Xmx2000m jar scripture.jar -task addpairs <Mandatory parameters>
## Score task
# java -Xmx2000m -jar scripture.jar -task score <Mandatory parameters>

#dot -Tps chr1.segments.bed.dot -o chr1.ps
	my $self = shift;
	my $ref = shift;
	$ref ||= 'ref';
	if (exists $self->{"software:scripture"}){
		return "";
	}
	my $scripture_cmd = $self->runTopHat('ref');
	$scripture_cmd .= qq(echo `date`; echo "run Scripture"\n);
	$scripture_cmd .= qq(export PATH="$self->{"setting:PATH"}":\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	my $reference=GATE::Error::checkPath($self->{"database:$ref"});
	$scripture_cmd .= qq(export REFERENCE="$reference"\n);
	my $heap = $self->{"setting:heap"};
	$scripture_cmd .= qq(export heap="$heap"\n);
	my $scripture = GATE::Error::checkPath($self->{"software:scripture"}) if (exists $self->{"software:scripture"});
	$scripture=correctJavaCmd($scripture,$heap,"./tmp_scripture");
	my $para = $self->{"setting:scripture"} if (exists $self->{"setting:scripture"});
	$scripture_cmd .= qq(export scripture="$scripture"\n);
	$scripture_cmd .= qq(export para="$para"\n);
	my $multi_t=$self->{"setting:multithreads"};
	$scripture_cmd .= qq(export multithreads=$multi_t\n);
	my $samtools = GATE::Error::checkPath($self->{"software:samtools"});
	$scripture_cmd .= qq(export samtools="$samtools"\n);
	my $which=`which igvtools`;chomp $which;
	my $igvtools=(-f $which && -e $which) ? $which : GATE::Error::checkPath($self->{"software:igvtools"});
	$scripture_cmd .= qq(export igvtools="$igvtools"\n);
	my $dot=`which dot`;chomp $dot;
	$dot=GATE::Error::checkPath($self->{"software:dot"}) if (exists $self->{"software:dot"});
	$scripture_cmd .= qq(export dot="$dot"\n);
	my $workdir=GATE::Error::checkPath($self->{"-workdir"});
	$scripture_cmd .= qq(export workdir="$workdir"\n);
	$scripture_cmd .= qq(cd \${workdir}\n);
	my $scripture_outdir = (exists $self->{"setting:scripture_outdir"}) ? GATE::Error::checkPath($self->{"setting:scripture_outdir"}) : "./scripture_out";
	$scripture_cmd .= qq(export scripture_outdir="$scripture_outdir"\n);
	$scripture_cmd .= qq([[ -d \${scripture_outdir} ]] || mkdir -p \${scripture_outdir}\n) if (!-d qq($self->{"-workdir"}/$scripture_outdir));
	$scripture_cmd .= qq(cd \${scripture_outdir}\n);
	if (-f "$reference.fai" || -f "$reference.dict") {
		if (-f "$reference.fai") {
			$scripture_cmd .= "cut -f 1,2 $reference.fai > $reference.size";
		}
		elsif (-f "$reference.dict") {
			$scripture_cmd .= qq(cut -f 1,2 $reference.dict | sed 's\/SN\\:\/\/;s\/LN\\:\/\/;' > $reference.size);
		}
	}
	my @libraries=sort keys %{$self->{'LIB'}};
	my $multi=0;
	foreach my $lib (@libraries) {
		#my %fq=getlibInput($self->{"LIB"}{$lib});
		$scripture_cmd .= qq([[ -d $lib ]] || mkdir -p $lib\n);
		$scripture_cmd .= qq(cd $lib\n);
		my $tophatbam = $self->{$lib}{"tophatbam"};
		my $tophatsam = "$1.sam" if ($tophatbam=~/(\S+)\.bam$/);
		$scripture_cmd .= qq(\${samtools} view $tophatbam > $tophatsam && );
		$scripture_cmd .= qq(\${igvtools} index $tophatsam && );
		$scripture_cmd .= qq(\${scripture} -alignment $tophatsam -out $lib.segments.txt -sizeFile $reference.size -chrSequences \${REFERENCE});
		$scripture_cmd .= qq( \${para}) if (defined $para);
		if ($multi % $multi_t==0) {
			$scripture_cmd .= "\n";
		}
		else {
			$scripture_cmd .= (exists $self->{"rule:multimode"}) ? " &\n " : "\n";
		}
		$scripture_cmd .= qq(cd ../\n);
		$multi++;
	}
	$scripture_cmd .= qq(cd ../\n);
	return $scripture_cmd;
}

sub runSOAPsplic ($) {
	
}


#########################################################
#                                                       #
#                 ChiP-seq/MeDIP-seq                    #
#                                                       #
#########################################################

sub runMACS ($) {
	#macs14 -t test.bowtie2.sam -c control.bowtie2.sam --bw 180 --gsize mm --verbose 3 --diag  -B -S --nomodel --name test_vs_control --mfold 10,30 2>macs.log
	my $self =shift;
	my $ref = shift;
	$ref ||= 'ref';
	if (!exists $self->{"software:macs"}) {
		return "";
	}
	my $macs=GATE::Error::checkPath($self->{"software:macs"});
	my $para=(exists $self->{'setting:macs'}) ? $self->{'setting:macs'} : " --bw=180 --verbose=3 --diag  -B -S --nomodel --mfold 10,30";
	my $macs_cmd = qq(echo `date`; echo "run MACS"\n);
	$macs_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$macs_cmd .= qq(export macs=$macs\n);
	$macs_cmd .= qq(export para="$para"\n);
	$macs_cmd .= qq(cd $self->{"-workdir"}\n);
	my $macs_outdir = (exists $self->{"setting:macs_outdir"}) ? $self->{"setting:macs_outdir"} : "macs";
	$macs_cmd .= qq([[ -d $macs_outdir ]] || mkdir $macs_outdir\n);
	$macs_cmd .= qq(cd $macs_outdir\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		$macs_cmd .= qq(echo `date`; echo "$lib"\n);
		$macs_cmd .= qq([[ -d $lib ]] || mkdir $lib\n);
		$macs_cmd .= qq(cd $lib\n);
		if (exists $self->{'LIB'}{$lib}{'INPUT'}{'CFILE'} && exists $self->{'LIB'}{$lib}{'INPUT'}{'TFILE'}) {
			my $tfile=$self->{'LIB'}{$lib}{'INPUT'}{'TFILE'};
			my $cfile=$self->{'LIB'}{$lib}{'INPUT'}{'CFILE'};
			$macs_cmd .= qq(\${macs} -t $tfile -c $cfile \${para} --name $lib\n);
		}
		$macs_cmd .= "cd ..\n";
	}
	return $macs_cmd;
}

#Fseq v1.84
sub runFseq ($) {
	#macs14 -t test.bowtie2.sam -c control.bowtie2.sam --bw 180 --gsize mm --verbose 3 --diag  -B -S --nomodel --name test_vs_control --mfold 10,30 2>macs.log
	my $self =shift;
	my $ref = shift;
	$ref ||= 'ref';
	if (!exists $self->{"software:macs"}) {
		return "";
	}
	my $fseq=GATE::Error::checkPath($self->{"software:fseq"});
	my $para=(exists $self->{'setting:fseq'}) ? $self->{'setting:fseq'} : " -of wig -t 8.0 -v ";
	my $fseq_cmd = qq(echo `date`; echo "run Fseq"\n);
	$fseq_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$fseq_cmd .= qq(export fseq=$fseq\n);
	$fseq_cmd .= qq(export para="$para"\n);
	$fseq_cmd .= qq(cd $self->{"-workdir"}\n);
	my $fseq_outdir = (exists $self->{"setting:fseq_outdir"}) ? $self->{"setting:fseq_outdir"} : "fseq";
	$fseq_cmd .= qq([[ -d $fseq_outdir ]] || mkdir $fseq_outdir\n);
	$fseq_cmd .= qq(cd $fseq_outdir\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		$fseq_cmd .= qq(echo `date`; echo "$lib"\n);
		$fseq_cmd .= qq([[ -d $lib ]] || mkdir $lib\n);
		$fseq_cmd .= qq(cd $lib\n);
 #fseq -v -of wig chr1.bed chr2.bed
 		if (exists $self->{'LIB'}{$lib}{'INPUT'}{'BED'}){
 			$fseq_cmd .= qq(\${fseq} $para );
 			$fseq_cmd .= join " ",@{$self->{'LIB'}{$lib}{'INPUT'}{'BED'}};
 		} 
		$fseq_cmd .= "cd ..\n";
	}
	return $fseq_cmd;
}

sub runRUM ($) {
	
}

sub runCisGenome ($) {
	
}

#########################################################
#                                                       #
#           Gene Prediction/Annotation                  #
#                                                       #
#########################################################

sub runGenomeThreader ($) {
#gth -showintronmaxlen 6000 -gcmaxgapwidth 6000 -o MW_v2c.gth.xml -xmlout -maskpolyatails -paralogs \
#-genomic ../../00.db/MW_v2c.fasta -cdna ../../00.db/cassava_cdna.fa \
#-protein ../../00.db/Euphorbiaceae_protein_sequence.fasta > 1.log 2> 2.log
	my $self=shift;
	if (!exists $self->{"software:gth"} || !-e $self->{"software:gth"})
	{
		return "";
	}
	my $gth=GATE::Error::checkPath($self->{"software:gth"});
	my $gth_cmd = qq(echo `date`; echo "run GenomeThreader"\n);
	$gth_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$gth_cmd .= qq(export gth="$gth"\n);
	my $denovo_genomics="";
	$denovo_genomics=$self->{"database:denovo"} if (exists $self->{"database:denovo"});
	if (!defined $denovo_genomics) {
		$gth_cmd .= $self->runTrinity();
		$denovo_genomics=$self->{"denovo_genomics"} if (exists $self->{"denovo_genomics"});
	}
	my $cdna=$self->{"database:cdna"} if (exists $self->{"database:cdna"});
	$gth_cmd .= qq(export cdna="$cdna"\n);
	my $protein=$self->{"database:protein"} if (exists $self->{"database:protein"});
	$gth_cmd .= qq(export protein="$protein"\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	$gth_cmd .= qq(cd $self->{"-workdir"}\n);
	$gth_cmd .= qq(mkdir gth\n) unless (-d qq($self->{"-workdir"}/gth));
	my $k=0;
	foreach my $lib(@libraries) {
		$gth_cmd .= qq(echo `date`; echo "$lib"\n);
		$gth_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$lib));
		$gth_cmd .= qq(cd $lib\n);
		if (exists $self->{"LIB"}{$lib}{'INPUT'}{'denovo'}) {
			$denovo_genomics = $self->{"LIB"}{$lib}{'INPUT'}{'denovo'};
		}
		elsif (exists $self->{"LIB"}{"$lib"}{'denovo'}) {
			$denovo_genomics = ${$self->{"LIB"}{$lib}{'denovo'}}[0];
		}
		elsif (exists $self->{"LIB"}{"$lib"}{'fa'}) {
			$denovo_genomics = ${$self->{"LIB"}{$lib}{'fa'}}[0];
		} else {
			$k++;
			next;
		}
		$gth_cmd .= qq(\${gth} -showintronmaxlen 6000 -gcmaxgapwidth 6000 -o gth.out.xml -xmlout -maskpolyatails -paralogs \\\n);
		$gth_cmd .= qq(-genomic $denovo_genomics\\\n);
		$gth_cmd .= qq(-cdns \${cdna}\\\n) if (defined $cdna);
		$gth_cmd .= qq(-protein \${protein}\\\n) if (defined $protein);
		$gth_cmd .= qq( > gth.log 2>&1\n);
		$gth_cmd .= qq(cd ../\n);
	}
	if ($k==@libraries) {
		$gth_cmd .= qq(\${gth} -showintronmaxlen 6000 -gcmaxgapwidth 6000 -o gth.out.xml -xmlout -maskpolyatails -paralogs \\\n);
		$gth_cmd .= qq(-genomic $denovo_genomics\\\n);
		$gth_cmd .= qq(-cdns \${cdna}\\\n) if (defined $cdna);
		$gth_cmd .= qq(-protein \${protein}\\\n) if (defined $protein);
		$gth_cmd .= qq( > gth.log 2>&1\n);
		$gth_cmd .= qq(cd ../\n);
	}
	return $gth_cmd;

}

#gm_es.pl denovo_genomics.fa
sub runGeneMark ($) {
	my $self=shift;
	if (!exists $self->{"software:gm"} || !-e $self->{"software:gm"})
	{
		return "";
	}
	my $gm=GATE::Error::checkPath($self->{"software:gm"});
	my $gm_cmd = qq(echo `date`; echo "run GeneMark"\n);
	$gm_cmd = qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$gm_cmd .= qq(export gm="$gm"\n);
	my $denovo_genomics="";
	$denovo_genomics=$self->{"database:denovo"} if (exists $self->{"database:denovo"});
	$denovo_genomics=$self->{"denovo_genomics"} if (exists $self->{"denovo_genomics"});
	if (!defined $denovo_genomics) {
		$gm_cmd .= $self->runTrinity();
		$denovo_genomics=$self->{"denovo_genomics"} if (exists $self->{"denovo_genomics"});
	}
	my @libraries=sort keys %{$self->{'LIB'}};
	my $k=0;
	foreach my $lib(@libraries) {
		$gm_cmd .= qq(echo `date`; echo "$lib"\n);
		$gm_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$lib));
		$gm_cmd .= qq(cd $lib\n);
		if (exists $self->{"LIB"}{$lib}{'INPUT'}{'denovo'}) {
			$denovo_genomics = $self->{"LIB"}{$lib}{'INPUT'}{'denovo'};
		}
		elsif (exists $self->{"LIB"}{"$lib"}{'denovo'}) {
			$denovo_genomics = ${$self->{"LIB"}{$lib}{'denovo'}}[0];
		}
		elsif (exists $self->{"LIB"}{"$lib"}{'fa'}) {
			$denovo_genomics = ${$self->{"LIB"}{$lib}{'fa'}}[0];
		} else {
			$k++;
			next;
		}
		$gm_cmd .= qq(\{gm} $denovo_genomics\n);
	}
	if ($k==@libraries) {
		$gm_cmd .= qq(\{gm} $denovo_genomics\n);
	}
	if (defined $denovo_genomics) {
		return $gm_cmd;
	}
}

sub runAUGUSTUS ($) {
	
}

sub runSNAP ($) {
	
}

sub runGenScan ($) {
	
}

sub runPASA ($) {
#PASA in the Context of a Complete Eukaryotic Annotation Pipeline
#(A) ab initio gene finding using a selection of the following software tools: GeneMarkHMM, FGENESH, Augustus, and SNAP, GlimmerHMM.
#(B) protein homology detection and intron resolution using the GeneWise software and the uniref90 non-redundant protein database.
#(C) alignment of known ESTs, full-length cDNAs, and most recently, Trinity RNA-Seq assemblies to the genome using GMAP
#(D) PASA alignment assemblies based on overlapping transcript alignments from step ( C)
#(E) use of EVidenceModeler (EVM) to compute weighted consensus gene structure annotations based on the above (A, B, C, D)
#(F) use of PASA to update the EVM consensus predictions, adding UTR annotations and models for alternatively spliced isoforms (leveraging D and E).
#(G) limited manual refinement of genome annotations (F) using Argo or Apollo
	my $self=shift;
}

sub runEVM ($) {
	my $self=shift;
}

#########################################################
#                                                       #
#                    Annotation                         #
#                                                       #
#########################################################

sub runBlast2GO ($) {
	my $self=shift;
}

#########################################################
#                                                       #
#                    Find ncRNA                         #
#                                                       #
#########################################################
sub runtRNAScan ($) {
	my $self=shift;
	my $tRNAscanSE=GATE::Error::checkPath($self->{"software:tRNAscan-SE"});
	if (!defined $tRNAscanSE || $tRNAscanSE eq "") {
		return "";
	}
	my $tRNA_cmd = qq(echo `date`; echo "run tRNAscan-SE"\n);
	$tRNA_cmd = qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$tRNA_cmd .= qq(export tRNAscanSE="$tRNAscanSE"\n);
	my $tRNApara=$self->{"setting:tRNAscan-SE"};
	$tRNA_cmd .= qq(export tRNApara="$tRNApara"\n);
	my $denovo_genomics="";
	$denovo_genomics=$self->{"database:denovo"} if (exists $self->{"database:denovo"});
	$denovo_genomics=$self->{"denovo_genomics"} if (exists $self->{"denovo_genomics"});
	if (!defined $denovo_genomics) {
		$tRNA_cmd .= $self->runTrinity();
		$denovo_genomics=$self->{"denovo_genomics"} if (exists $self->{"denovo_genomics"});
	}
	my @libraries=sort keys %{$self->{'LIB'}};
	my $k=0;
	foreach my $lib(@libraries) {
		$tRNA_cmd .= qq(echo `date`; echo "$lib"\n);
		$tRNA_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$lib));
		$tRNA_cmd .= qq(cd $lib\n);
		if (exists $self->{"LIB"}{$lib}{'INPUT'}{'denovo'}) {
			$denovo_genomics = $self->{"LIB"}{$lib}{'INPUT'}{'denovo'};
		}
		elsif (exists $self->{"LIB"}{"$lib"}{'denovo'}) {
			$denovo_genomics = ${$self->{"LIB"}{$lib}{'denovo'}}[0];
		}
		elsif (exists $self->{"LIB"}{"$lib"}{'fa'}) {
			$denovo_genomics = ${$self->{"LIB"}{$lib}{'fa'}}[0];
		} else {
			$k++;
			next;
		}
		$tRNA_cmd .= qq(\{tRNAscanSE} \${tRNApara} -o $denovo_genomics\.tRNA -f $denovo_genomics\.structure $denovo_genomics\n);
	}
	if ($k==@libraries) {
		$tRNA_cmd .= qq(\{tRNAscanSE} \${tRNApara} -o $denovo_genomics\.tRNA -f $denovo_genomics\.structure $denovo_genomics\n);
	}
	if (defined $denovo_genomics) {
		return $tRNA_cmd;
	}
}

sub runInferal ($) {
	my $self=shift;
}

#########################################################
#                                                       #
#                    Gene Fusion                        #
#                                                       #
#########################################################

sub runTophatFusion ($) {
	my $self = shift;
}

sub runSOAPfusion ($) {
	my $self = shift;
}

sub runSOAPfuse ($) {
	my $self = shift;
}

sub runCRAC($) {
	my $self = shift;
	my $ref = shift;
	$ref ||= "ref";
	if (!exists $self->{"software:crac"} || !defined $self->{"software:crac"}){
		return "";
	}
	my $crac = $self->{"software:crac"};
	my $cracpara = $self->{"setting:crac"};
	$cracpara .= qq( -k 22 \n) if ($cracpara!~/\-k/);
	my $reference=GATE::Error::checkPath($self->{"database:$ref"});
	my $crac_cmd = qq(echo `date`; echo "run CRAC"\n);;
	$crac_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$crac_cmd .= qq(export REFERENCE="$reference"\n);
	$crac_cmd .= qq(export crac=$crac\n);
	$crac_cmd .= qq(export seecerpara="$cracpara"\n);
	$crac_cmd .= qq(export qc_outdir=$self->{"setting:qc_outdir"}\n);
	$crac_cmd .= qq(cd $self->{"-workdir"}\n);
	$crac_cmd .= qq([[ -d \${qc_outdir} || mkdir \${qc_outdir}\n) if (!-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}));
	$crac_cmd .= qq(cd $self->{"setting:qc_outdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		$crac_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"setting:qc_outdir"}/$lib));
		$crac_cmd .= "cd $lib\n";
		my %fq=getlibInput($self->{"LIB"}{$lib});
		if (exists $fq{1} && exists $fq{2}){
			for (my $i=0;$i<@{$fq{2}};$i++) {
				my $fq1=${$fq{1}}[$i];
				my $fq2=${$fq{2}}[$i];
				my $bam=(@{$fq{2}}>1)?"$lib.pair.".($i+1).".sam":"$lib.pair.bam";
				my $chimera=(@{$fq{2}}>1)?"$lib.pair.".($i+1).".chimera":"$lib.pair.chimera";
				$crac_cmd .= qq(\${crac} \${cracpara} -i \$REFERENCE -r $fq1 $fq2 -o $bam );
				$crac_cmd .= qq( --paired-end-chimera $chimera ) if ($cracpara=~/chimera/);
				$crac_cmd .= qq(--nb-threads $self->{'setting:multithreads'}) if (exists $self->{'setting:multithreads'});
				$crac_cmd .= qq(\n);
			}
		}
		if (exists $fq{0} && @{$fq{0}}>0){
			for (my $i=0;$i<@{$fq{0}};$i++) {
				my $fq=${$fq{0}}[$i];
				my $bam=(@{$fq{0}}>1)?"$lib.single.".($i+1).".sam":"$lib.single.bam";
				my $chimera=(@{$fq{0}}>1)?"$lib.single.".($i+1).".chimera":"$lib.single.chimera";
				$crac_cmd .= qq(\${crac} \${cracpara} -i \$REFERENCE -r $fq -o $bam );
				$crac_cmd .= qq( --chimera $chimera ) if ($cracpara=~/chimera/);
				$crac_cmd .= qq(\n);
			}
		}
		$crac_cmd .= "cd ..\n";
	}
	return $crac_cmd;
}

#########################################################
#                                                       #
#                    Gene Editing                       #
#                                                       #
#########################################################

sub runRepeatMasker ($) {
	
}

sub runTRF ($) {
	
}


#########################################################
#                                                       #
#                    Phylogenetic                       #
#                                                       #
#########################################################

sub runClustalW2($) {
	my $self=shift;
	my $which=`which clustalw2`;chomp $which;
	my $clustalw2=(-f $which && -e $which) ? $which : GATE::Error::checkPath($self->{"software:clustalw2"});
	if (!defined $clustalw2)
	{
		return "";
	}
#clustalw2 -INFILE=chloroplast.ext -ALIGN -TREE -PIM -TYPE=DNA -OUTPUT=NEXUS -STATS=chloroplast_stats.log > chloroplast.log 2>&1 &
	my $clustalw2_cmd=qq(echo `date`; echo "run ClustalW2"\n);;
	$clustalw2_cmd.= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$clustalw2_cmd.=qq(export clustalw2=$clustalw2\n);
	my $para=(exists $self->{"setting:clustalw2"})?$self->{"setting:clustalw2"}:"-ALIGN -TREE -PIM -TYPE=DNA -OUTPUT=NEXUS -STATS=stats.log";
	$clustalw2_cmd.=qq(export clustalw2_para=$para\n);
	my $phylogen_outdir = (exists $self->{"setting:phylogen_outdir"}) ? $self->{"setting:phylogen_outdir"} : "phylogen";
	my $workdir=GATE::Error::checkPath($self->{"-workdir"});
	$clustalw2_cmd.=qq(export workdir=$workdir\n);
	$clustalw2_cmd.=qq(export phylgen_outdir=$phylogen_outdir\n);
	$clustalw2_cmd.=qq(cd \${wordir}\n);
	$clustalw2_cmd.=qq(mkdir \${phylgen_outdir}\n) if (!-d qq($self->{"-workdir"}/$phylogen_outdir));
	$clustalw2_cmd.=qq(cd \${phylgen_outdir}\n);
	if (exists $self->{clustalw2_INFILE})
	{
		my $infile=GATE::Error::checkPath($self->{clustalw2_INFILE});
		$clustalw2_cmd.=qq(\${clustalw2} \${clustalw2_para} -INFILE=$infile\n);
		my $prefix=$1 if ($infile=~/([^\/\s]\S+)\.[^\.\s]+$/);
		$self->{"NEXUS"}="$workdir/$phylogen_outdir/$prefix.nxs";
	}
	else
	{
		my @libraries=sort keys %{$self->{'LIB'}};
		foreach my $lib(@libraries) {
			$clustalw2_cmd.="mkdir $lib\n";
			$clustalw2_cmd.="cd $lib\n";
			my @file=();
			push @file,@{$self->{"LIB"}{$lib}{'fa'}} if (exists $self->{"LIB"}{$lib}{'fa'});
			push @file,@{$self->{"LIB"}{$lib}{'fasta'}} if (exists $self->{"LIB"}{$lib}{'fasta'});
			push @file,@{$self->{"LIB"}{$lib}{'INFILE'}} if (exists $self->{"LIB"}{$lib}{'INFILE'});
			if (@file>0) {
				for (my $i=0;$i<@file;$i++)
				{
					$clustalw2_cmd.="\${clustalw2} \${clustalw2_para} -INFILE=$file[$i]\n";
					my $prefix=$1 if ($file[$i]=~/([^\/\s]\S+)\.[^\.\s]+$/);
					push @{$self->{$lib}{"NEXUS"}},"$workdir/$phylogen_outdir/$prefix.nxs";
				}
			}
			$clustalw2_cmd.="cd ..\n";
		}
	}
	$clustalw2_cmd.="cd ..\n";
	return $clustalw2_cmd;
}

sub runMrBayes
{
	my $self=shift;
	my $which = `which mb`;chomp $which;
	my $mb = (-f $which && -e $which) ? $which : GATE::Error::checkPath($self->{"software:mb"});
	if (!defined $mb) {
		return "";
	}
	my $mb_cmd=qq(echo `date`; echo "run MrBayes"\n);
	$mb_cmd .= qq(export PATH=$self->{"setting:PATH"}:\$PATH\n) if (exists $self->{"setting:PATH"} && $self->{"setting:PATH"}!~/\/usr\/local\/bin/);
	$mb_cmd.=qq(export mb=$mb\n);
	my $conf=(exists $self->{"setting:mb"})?$self->{"setting:mb"}:"\tlset nst=6 rates=invgamma;\n\tmcmc ngen=20000000;\n\tsump relburnin=yes burninfrac=0.25;\n\tsumt relburnin=yes burninfrac=0.25;\n";
	$mb_cmd.=qq(export mb_conf=$conf\n);
	my $phylogen_outdir = (exists $self->{"setting:phylogen_outdir"}) ? $self->{"setting:phylogen_outdir"} : "phylogen";
	my $workdir=GATE::Error::checkPath($self->{"-workdir"});
	$mb_cmd.=qq(export workdir=$workdir\n);
	$mb_cmd.=qq(export phylgen_outdir=$phylogen_outdir\n);
	$mb_cmd.=qq(cd \${wordir}\n);
	$mb_cmd.=qq(mkdir \${phylgen_outdir}\n) if (!-d qq($self->{"-workdir"}/$phylogen_outdir));
	$mb_cmd.=qq(cd \${phylgen_outdir}\n);
	if (exists $self->{"NEXUS"}) {
		$mb_cmd.=qq(echo "begin mrbayes;\n\texec $self->{"NEXUS"}\n\${mb_conf}\nend;\n" > batch.txt\n);
		$mb_cmd.=qq(\${mb} batch.txt\n);
	}else {
		my @libraries=sort keys %{$self->{'LIB'}};
		foreach my $lib(@libraries) {
			$mb_cmd.="mkdir $lib\n";
			$mb_cmd.="cd $lib\n";
			for (my $i=0;$i<@{$self->{$lib}{"NEXUS"}};$i++)
			{
				$mb_cmd.=qq(echo "begin mrbayes;\\texec ${$self->{$lib}{"NEXUS"}}[$i]\n\${mb_conf}\nend;\n" > $lib-$i.batch.txt\n);
				$mb_cmd.=qq(\${mb} $lib-$i.batch.txt\n);
			}
			$mb_cmd.="cd ..\n";
		}
	}
	$mb_cmd.="cd ..\n";
	return $mb_cmd;
}

sub runCluster ($) {
#cluster -f tabfile -l -cg m -ca a -na -e 7	
}

sub runCirCOS ($) {
	
}

sub runANFO ($) {
	
}

sub mapDamage ($) {
	
}

1;
