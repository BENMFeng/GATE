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
$VERSION = "0.9i,05-08-2013";
package GATE::Element;
#package GATE::Extension;
use FindBin qw($Bin $Script);
use lib "$FindBin::Bin/../lib";
use File::Basename qw(basename dirname);
use Cwd;


sub parseConfig($) {
	my $self = shift;
	my $file=checkPath($self->{'-config'});
	my ($name,$index,$barcode,$lib,$mergepe)=("","","","",""); 
	my %RG;
	my %countSample;
	open (IN,$file) || die "Can't open config file:$file for reading\n";
	while(<IN>) {
		chomp;
		next if ($_ eq "" || $_=~/^\#/);
		if (/\[(.*)\]/) {
			$name=$1;
			if ($name eq "LIB")
			{
				($index,$barcode,$lib,$mergepe)=("","","","");
				delete @RG{keys %RG};
				%RG=();
			}
		}
		if ($name eq "INPUT") {
			if (/^([A-Z]{2})\=([^\=]+)/) {
				my ($rg,$info)=($1,$2);
				$RG{$rg}=$info;
				$lib=$1 if (/LB\=([^\=]+)/);
			} elsif (/([^\=]+)\=([^\=]+)/) {
				my ($lb,$path)=($1,$2);
				$self->{'LIB'}{$lib}{$name}{$lb}++;
				$self->{$lib}{$lb}[$self->{'LIB'}{$lib}{$name}{$lb}-1]=checkPath($path);
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
					${$self->{$name}{$lib}{$type}}[$countSample{"$name $lib $type"}-1]=checkPath($path);
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
		} else {
			if (/([^\=]+)\=([^\=]+)/) {
				my ($lib,$path)=($1,$2);
				$path=checkPath($path) if ($name=~/database/ || $name=~ /software/);
				$self->{"$name:$lib"} = $path;
			}
		}
	}
	close IN;
	return bless($self);
}

sub checkPath{
	my $path=shift;
	#print STDERR "path=$path\n";
	if (defined $path && $path !~ /^\//) {
		if($path=~/\s([^\/\s]+\/\S+)/) {
			my $d=$1;
			my ($dir,$file)=($1,$2) if ($d=~/(\S+\/)([^\/\s]+)$/);
			if (-d $dir && -f $d) {
				my $rawdir=`pwd`;
				chomp $rawdir;
				chdir($rawdir);
				chdir($dir);
				my $targetdir=`pwdf`;
				chomp $targetdir;
				chdir($rawdir);
				$path="$targetdir/$file";
			} else {
				print STDERR get_time()."\t\t$path error: folder $dir or file $d doesn't exist\n";
				die();
			}
		}
	} else {
		if (defined $path &&  $path=~/^(\/\S+)/) {
			if (!-d $1 && !-f $1) {
				print STDERR get_time()."\t\t$path error: folder or file doesn't exist\n";
				die();
			}
		}
	}
	return $path;
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

#########################################################
#                                                       #
#                    Data Processing                    #
#                                                       #
#########################################################

sub selectIdxFastq ($) {
	my $self = shift;
	if (!exists $self->{"software:selectIdxFastq"} || !defined $self->{"software:selectIdxFastq"} || (!exists $self->{'idx'} && !exists $self->{'bar'})) {
		return "";
	}
	my $selectIdxFastq=$self->{"software:selectIdxFastq"};
	my $para=(exists $self->{"CustomSetting:selectIdxFastq"})?$self->{"CustomSetting:selectIdxFastq"}:"-mis 1 -qual 30";
	$para=~s/\-.+prefix\s+\S+//;
	my $Idx_cmd = qq(echo `date`; echo "run selectIdxFastq"\n);
	$Idx_cmd .= qq(cd $self->{"-workdir"}\n);
	$Idx_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$Idx_cmd .= qq(export selectIdxFastq=$selectIdxFastq\n);
	$Idx_cmd .= qq(export selectIdxFastq_para="$para"\n);
	$Idx_cmd .= qq(mkdir $self->{"CustomSetting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}));
	$Idx_cmd .= qq(cd $self->{"CustomSetting:qc_outdir"}\n);
	my $Idx_cmd_multi_head=$Idx_cmd;
	my $Idx_cmd_multi="";
	my $withIdx=0;
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		next if (!exists $self->{'idx'}{$lib} && !exists $self->{'bar'}{$lib});
		$Idx_cmd_multi .= $Idx_cmd_multi_head;
		$Idx_cmd .= qq(echo `date`; echo "$lib"\n);
		$Idx_cmd .= qq([[ -d $lib ]] || mkdir $lib\n)  unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
		$Idx_cmd_multi .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
		$Idx_cmd .= "cd $lib\n";
		$Idx_cmd_multi .= "cd $lib && ";
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
						die qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out2 is existent!\n) if (-f qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out1));
						die qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out2 is existent!\n) if (-f qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out2));
						${$self->{"LIB"}{$lib}{$lbmark}}[$k]=qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out1);
						${$self->{"LIB"}{$lib}{"$lb$j"}}[$k]=qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out2);
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
							die qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out1 is existent!\n) if (-f qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out1));
							${$self->{"LIB"}{$lib}{$lbmark}}[$k]=qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out1);
							$withIdx++;
						}
					}
				}
			}
		}
		$Idx_cmd .= "cd ..\n";
		$Idx_cmd_multi .= "cd .. && ";
		$Idx_cmd_multi .= "cd .. ";
		if ($withIdx % $self->{"CustomSetting:multithreads"}!=0)
		{
			$Idx_cmd_multi .= " &\n";
		} else {
			$Idx_cmd_multi .= "\n";
		}
	}
	$Idx_cmd .= "cd ..\n";
	if ($withIdx>0){
		if (exists $self->{"CustomSetting:multimode"}) {
			chomp $Idx_cmd_multi;
			$Idx_cmd_multi=~s/\&+$//;
			$Idx_cmd_multi.="\n";
			$Idx_cmd_multi.=print_check_process('selectIdxFastq');
			return $Idx_cmd_multi;
		}else{
			return $Idx_cmd;
		}
	}
}

sub mergeOverlapPE($) {
	my $self = shift;
	if ( (!exists $self->{"software:mergeOverlapPE"} && (!exists $self->{"software:bwa-pemerge"})) || (!exists $self->{'pemerge'}) ) {
		return "";
	}
	if (exists $self->{"software:mergeOverlapPE"}) {
		my $mergeOverlapPE=$self->{"software:mergeOverlapPE"};
		my $moppara = $self->{"CustomSetting:mergeOverlapPE"};
		my $mop_cmd = qq(echo `date`; echo "run mergeOverlapPE"\n);
		$mop_cmd .= qq(cd $self->{"-workdir"}\n);
		$mop_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
		$mop_cmd .= qq(mkdir $self->{"CustomSetting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}));
		$mop_cmd .= qq(cd $self->{"CustomSetting:qc_outdir"}\n);
		my $mop_cmd_multi_head=$mop_cmd;
		my $mop_cmd_multi="";
		my $withPE=0;
		my @libraries=sort keys %{$self->{'LIB'}};
		foreach my $lib(@libraries) {
			next if (!exists $self->{'pemerge'}{$lib});
			$mop_cmd_multi.= $mop_cmd_multi_head;
			$mop_cmd_multi .= qq(echo `date`; echo "$lib"\n);
			$mop_cmd .= qq(echo `date`; echo "$lib"\n);
			$mop_cmd .= qq([[ -d $lib ]] || mkdir $lib\n)  unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
			$mop_cmd_multi .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
			$mop_cmd .= "cd $lib\n";
			$mop_cmd_multi .= "cd $lib && ";
			my %fq=getlibSeq($self->{"LIB"}{$lib});
			my @Reads1=();
			@Reads1=@{$fq{1}} if (exists $fq{1});
			my @Reads2=();
			@Reads2=@{$fq{2}} if (exists $fq{2});
			if (@Reads1==0 || @Reads2==0)
			{
				$mop_cmd .= "cd ../\n";
				$mop_cmd_multi .= "cd ../";
				if ($withPE % $self->{"CustomSetting:multithreads"}!=0)
				{
					$mop_cmd_multi .= " &\n";
				} else {
					$mop_cmd_multi .= "\n";
				}
				next;
			}
			for (my $i=0;$i<@Reads1;$i++)
			{
				if (exists $self->{$lib}{'fq1'}{$i}{"MergePE"} && ($self->{$lib}{'fq1'}{$i}{"MergePE"} =~ /TRUE/i || $self->{$lib}{'fq1'}{$i}{"MergePE"} =~ /Yes/i)) {
					my $prefix=(@Reads1>1)?"$lib-".($i+1):$lib;
					$mop_cmd .= "$mergeOverlapPE $Reads1[$i] $Reads2[$i] -prefix $prefix $moppara\n";
					$mop_cmd_multi .= "$mergeOverlapPE $Reads1[$i] $Reads2[$i] -prefix $prefix $moppara && ";
					${$self->{"LIB"}{$lib}{'fq1'}}[$i]=qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$prefix\_R1.fastq);
					${$self->{"LIB"}{$lib}{'fq2'}}[$i]=qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$prefix\_R2.fastq);
					push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$prefix\_merged.fastq);
					foreach my $k(keys %{$self->{$lib}{'fq1'}{$i}})
					{
						$self->{$lib}{'fq'}{@{$self->{"LIB"}{$lib}{'fq'}}-1}{$k}=$self->{$lib}{'fq1'}{$i}{$k};
					}
					$withPE++;
				}
			}
			$mop_cmd .= "cd ..\n";
			$mop_cmd_multi .= "cd .. && ";
			$mop_cmd_multi .= "cd ..";
			if ($withPE % $self->{"CustomSetting:multithreads"} != 0)
			{
				$mop_cmd_multi .= " &\n";
			} else {
				$mop_cmd_multi .= "\n";
			}
		}
		$mop_cmd .= "cd ..\n";
		if ($withPE>0){
			if (exists $self->{"CustomSetting:multimode"}) {
				chomp $mop_cmd_multi;
				$mop_cmd_multi =~ s/\&+$//;
				$mop_cmd_multi .= "\n";
				$mop_cmd_multi .= print_check_process('mergeOverlapPE');
				return $mop_cmd_multi;
			}else{
				return $mop_cmd;
			}
		}
	}elsif (exists $self->{"software:bwa-pemerge"}) {
		my $bwapemerge=$self->{"software:bwa-pemerge"};
		my $moppara = (exists $self->{"CustomSetting:bwa-pemerge"}) ? $self->{"CustomSetting:bwa-pemerge"} : "";
		my $mop_cmd = qq(echo `date`; echo "run bwa-pemerge"\n);
		$mop_cmd .= qq(cd $self->{"-workdir"}\n);
		$mop_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
		$mop_cmd .= qq(mkdir $self->{"CustomSetting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}));
		$mop_cmd .= qq(cd $self->{"CustomSetting:qc_outdir"}\n);
		my $mop_cmd_multi_head=$mop_cmd;
		my $mop_cmd_multi="";
		my $withPE=0;
		my @libraries=sort keys %{$self->{'LIB'}};
		foreach my $lib(@libraries) {
			next if (!exists $self->{'pemerge'}{$lib});
			$mop_cmd_multi.= $mop_cmd_multi_head;
			$mop_cmd_multi .= qq(echo `date`; echo "$lib"\n);
			$mop_cmd .= qq(echo `date`; echo "$lib"\n);
			$mop_cmd .= qq([[ -d $lib ]] || mkdir $lib\n)  unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
			$mop_cmd_multi .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
			$mop_cmd .= "cd $lib\n";
			$mop_cmd_multi .= "cd $lib && ";
			my %fq=getlibSeq($self->{"LIB"}{$lib});
			my @Reads1=();
			@Reads1=@{$fq{1}} if (exists $fq{1});
			my @Reads2=();
			@Reads2=@{$fq{2}} if (exists $fq{2});
			if (@Reads1==0 || @Reads2==0)
			{
				$mop_cmd .= "cd ../\n";
				$mop_cmd_multi .= "cd ../";
				if ($withPE % $self->{"CustomSetting:multithreads"}!=0)
				{
					$mop_cmd_multi .= " &\n";
				} else {
					$mop_cmd_multi .= "\n";
				}
				next;
			}
			for (my $i=0;$i<@Reads1;$i++)
			{
				if (exists $self->{$lib}{'fq1'}{$i}{"MergePE"} && ($self->{$lib}{'fq1'}{$i}{"MergePE"} =~ /TRUE/i || $self->{$lib}{'fq1'}{$i}{"MergePE"} =~ /Yes/i)) {
					my $prefix=(@Reads1>1)?"$lib-".($i+1):$lib;
					$mop_cmd .= "$bwapemerge $moppara -m $Reads1[$i] $Reads2[$i] > $prefix\_merged.fastq\n";
					$mop_cmd_multi .= qq($bwapemerge $moppara -m -t $self->{"CustomSetting:multithreads"} $Reads1[$i] $Reads2[$i] > $prefix\_merged.fastq && );
					$mop_cmd .= qq($bwapemerge $moppara -u $Reads1[$i] $Reads2[$i] | awk '\{if \(NR\%8<4\)\{print \$0 > "$prefix\_R1.fastq"\}else\{print \$0 > "$prefix\_R2.fastq"\}\}'\n);
					$mop_cmd_multi .= qq($bwapemerge $moppara -u -t $self->{"CustomSetting:multithreads"} $Reads1[$i] $Reads2[$i] | awk '\{if \(NR\%8<4\)\{print \$0 > "$prefix\_R1.fastq"\}else\{print \$0 > "$prefix\_R2.fastq"\}\}' && );
					${$self->{"LIB"}{$lib}{'fq1'}}[$i]=qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$prefix\_R1.fastq);
					${$self->{"LIB"}{$lib}{'fq2'}}[$i]=qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$prefix\_R2.fastq);
					push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$prefix\_merged.fastq);
					foreach my $k(keys %{$self->{$lib}{'fq1'}{$i}})
					{
						$self->{$lib}{'fq'}{@{$self->{"LIB"}{$lib}{'fq'}}-1}{$k}=$self->{$lib}{'fq1'}{$i}{$k};
					}
					$withPE++;
				}
			}
			$mop_cmd .= "cd ..\n";
			$mop_cmd_multi .= "cd .. && ";
			$mop_cmd_multi .= "cd ..";
			if ($withPE % $self->{"CustomSetting:multithreads"} != 0)
			{
				$mop_cmd_multi .= " &\n";
			} else {
				$mop_cmd_multi .= "\n";
			}
		}
		$mop_cmd .= "cd ..\n";
		if ($withPE>0){
			if (exists $self->{"CustomSetting:multimode"}) {
				chomp $mop_cmd_multi;
				$mop_cmd_multi =~ s/\&+$//;
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
		my $SolexaQA=checkPath($self->{"software:SolexaQA"});
		my $para=$self->{"CustomSetting:SolexaQA"};
		#SolexaQA.pl ../BF3.1.fq -v -m -s 10000 -b -sanger -d ./ 
		my $qa_cmd = qq(echo `date`; echo "run SolexaQA"\n);
		$qa_cmd .= qq(cd $self->{"-workdir"}\n);
		$qa_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
		$qa_cmd .= qq(mkdir $self->{"CustomSetting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}));
		$qa_cmd .= qq(cd $self->{"CustomSetting:qc_outdir"}\n);
		my $qa_cmd_multi_head=$qa_cmd;
		my $qa_cmd_multi="";
		my $multi=0;
		my @libraries=sort keys %{$self->{'LIB'}};
		foreach my $lib(@libraries){
			$qa_cmd_multi.= $qa_cmd_multi_head;
			$qa_cmd_multi .= qq(echo `date`; echo "$lib"\n);
			$qa_cmd .= qq(echo `date`; echo "$lib"\n);
			$qa_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
			$qa_cmd_multi .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
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
							if (!-f qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/QA/$fq))
							{
								$qa_cmd .= "gzip -cd $reads > $fq\n";
								$qa_cmd_multi .= "gzip -cd $fq && ";
							}
							$qa_cmd .= "$SolexaQA $fq ";
							$qa_cmd .= ($para=~/\-d/) ? " $para\n" : "$para -d QA\n";
							$qa_cmd_multi .= "$SolexaQA $fq ";
							$qa_cmd_multi .= ($para=~/\-d/) ? " $para && " : "$para -d QA &&";
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
			$qa_cmd_multi .= "cd .. && ";
			$qa_cmd_multi .= "cd .. ";
			$multi++;
			if ($multi % $self->{"CustomSetting:multithreads"}!=0){
				$qa_cmd_multi .= "&\n";
			}else{
				$qa_cmd_multi .= "\n";
			}
		}
		$qa_cmd .= "cd ..\n";
		if (exists $self->{"CustomSetting:multimode"} && $self->{"CustomSetting:multimode"} =~ /y/i) {
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
		my $check_fastq=checkPath($self->{"software:check_fastq"});
		my $para=$self->{"CustomSetting:check_fastq"};
		my $distribute_fqcheck=$self->{"CustomSetting:distribute_fqcheck"} if (exists $self->{"CustomSetting:distribute_fqcheck"});
		#SolexaQA.pl ../BF3.1.fq -v -m -s 10000 -b -sanger -d ./ 
		my $qa_cmd = qq(echo `date`; echo "run check_fastq"\n);
		$qa_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
		$qa_cmd .= qq(export check_fastq="$check_fastq"\n);
		$qa_cmd .= qq(cd $self->{"-workdir"}\n);
		$qa_cmd .= qq(mkdir $self->{"CustomSetting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}));
		$qa_cmd .= qq(cd $self->{"CustomSetting:qc_outdir"}\n);
		my $qa_cmd_multi_head=$qa_cmd;
		my $qa_cmd_multi="";
		my $multi=0;
		my @libraries=sort keys %{$self->{'LIB'}};
		foreach my $lib(@libraries){
			$qa_cmd_multi.= $qa_cmd_multi_head;
			$qa_cmd_multi .= qq(echo `date`; echo "$lib"\n);
			$qa_cmd .= qq(echo `date`; echo "$lib"\n);
			$qa_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
			$qa_cmd_multi .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
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
							if (!-f qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/QA/$fq))
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
			if ($multi % $self->{"CustomSetting:multithreads"}!=0){
				$qa_cmd_multi .= "&\n";
			}else{
				$qa_cmd_multi .= "\n";
			}
		}
		$qa_cmd .= "cd ..\n";
		if (exists $self->{"CustomSetting:multimode"} && $self->{"CustomSetting:multimode"} =~ /y/i) {
			chomp $qa_cmd_multi;
			$qa_cmd_multi=~s/\&+$//;
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
	my $filterPCRdup=checkPath($self->{"software:filterPCRdup"});
	my $para=(exists $self->{"CustomSetting:filterPCRdup"})? $self->{"CustomSetting:filterPCRdup"} : "";
	my $fltpcr_cmd = qq(echo `date`; echo "run filterPCRdup"\n);
	$fltpcr_cmd .= qq(cd $self->{"-workdir"}\n);
	$fltpcr_cmd .= qq(mkdir $self->{"CustomSetting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}));
	$fltpcr_cmd .= qq(cd $self->{"CustomSetting:qc_outdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib (@libraries) {
		$fltpcr_cmd .= qq(echo `date`; echo "$lib"\n);
		$fltpcr_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
		$fltpcr_cmd .= "cd $lib\n";
		my %fq=getlibSeq($self->{"LIB"}{$lib});
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
				if (!-f qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$fq1))
				{
					$fltpcr_cmd .= "gzip -cd $reads1 > $fq1\n" if ($reads1=~/gz$/);
				}
			}
			my $fqname2=(split /\//,$reads2)[-1] if (defined $reads2);
			my $fq2=$reads2 if (defined $reads2);
			if ( (defined $reads2) && ( ($reads2 !~ /fastq$/i && $reads2 !~ /fq$/i) || ($reads2 =~ /gz$/) ) ) {
				$fq2=$1 if ($fqname2=~/(\S+)\.gz/);
				if (!-f qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$fq2))
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
				${$self->{"LIB"}{$lib}{'fq1'}}[$j]=qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out1);
				${$self->{"LIB"}{$lib}{'fq2'}}[$j]=qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out2);
			}
			elsif (@Reads1>0)
			{
				${$self->{"LIB"}{$lib}{'fq1'}}[$j]=qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out1);
			}
		}
		for (my $j=0;$j<@Reads3;$j++) {
			my $reads3=$Reads3[$j];
			my $fqname3=(split /\//,$reads3)[-1];
			my $fq3=$reads3;
			if (($reads3 !~ /fastq$/i && $reads3 !~ /fq$/i) || ($reads3 =~ /gz$/)) {
				$fq3=$1 if ($fqname3=~/(\S+)\.gz/);
				if (!-f qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$fq3))
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
			${$self->{"LIB"}{$lib}{'fq'}}[$j]=qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$out3);
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
	my $scanAP=checkPath($self->{"software:scanAP"});
	my $fastqcut=$self->{"software:fastqcut"} if (exists $self->{"software:fastqcut"});
	my $fastqcut_para = $self->{"CustomSetting:fastqcut"} if (exists $self->{"CustomSetting:fastqcut"});
	my $align_matrix=$self->{"CustomSetting:align.mat"} if (exists $self->{"CustomSetting:align.mat"});
	my $para=(exists $self->{"CustomSetting:scanAP"})? $self->{"CustomSetting:scanAP"} : "";
	my $AP=checkPath($self->{"database:AP"});
	my $trim_seq=checkPath($self->{"software:trim_seq"}) if (exists $self->{"software:trim_seq"});
	my $trim_seq_para=$self->{"CustomSetting:trim_seq"};
	my $fltap_cmd = qq(echo `date`; echo "run scanAP"\n);
	$fltap_cmd .=qq(cd $self->{"-workdir"}\n);
	$fltap_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$fltap_cmd .= qq(export AP=$AP\n);
	$fltap_cmd .= qq(mkdir $self->{"CustomSetting:qc_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}));
	$fltap_cmd .= qq(cd $self->{"CustomSetting:qc_outdir"}\n);
	my $fltap_cmd_multi_head=$fltap_cmd;
	my $fltap_cmd_multi="";
	my $multi=0;
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib (@libraries) {
		$fltap_cmd_multi .= $fltap_cmd_multi_head;
		$fltap_cmd_multi .= qq(echo `date`; echo "$lib"\n);
		$fltap_cmd .= qq(echo `date`; echo "$lib"\n);
		$fltap_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
		$fltap_cmd_multi .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
		$fltap_cmd .= "cd $lib\n";
		$fltap_cmd_multi .= "cd $lib && ";
		$fltap_cmd .= qq(cp $align_matrix ./\n) unless (!defined $align_matrix || -f qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/align.mat));
		$fltap_cmd_multi .=  qq(cp $align_matrix ./ && ) unless (!defined $align_matrix || -f qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/align.mat));
		
		my %Reads;
		foreach my $i(sort keys %{$self->{"LIB"}{$lib}}) {
			my ($reads1,$reads2,$filter1,$filter2)=("","","","");
			foreach my $reads(@{$self->{"LIB"}{$lib}{$i}}) {
				my $fqname=(split /\//,$reads)[-1];
				my $fq=$reads;
				if (($reads !~ /fastq$/i && $reads !~ /fq$/i) || ($reads =~ /gz$/)) {
					$fq=$1 if ($fqname=~/(\S+)\.gz/);
					$fqname = $fq;
					if (!-f qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$fq))
					{
						$fltap_cmd .= "gzip -cd $reads > $fq\n" if ($reads=~/gz$/);
						$fltap_cmd_multi .= "gzip -cd $reads > $fq && " if ($reads=~/gz$/);
					}
				}
				
				my ($stat,$detail)=("$fqname.stat","$fqname.detail");
				if ($fqname=~/(\S+)\.[^\.\s]+$/) {
					$stat="$1.stat";
					$detail="$1.detail";
				}
				if ( ( (!exists $self->{"CustomSetting:reuse"}) || (exists $self->{"CustomSetting:reuse"} && $self->{"CustomSetting:reuse"} =~/N/i) ) && !-f qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib/$detail))
				{
					$fltap_cmd .= "$scanAP -i $fq -a \$AP -s $stat -d $detail $para\n";
					$fltap_cmd_multi .= "$scanAP -i $fq -a \$AP -s $stat -d $detail $para && ";
				}
				$fltap_cmd .= "$trim_seq $trim_seq_para -trim_detail $detail $fq\n" if (exists $self->{"software:trim_seq"});
				$fltap_cmd_multi .= "$trim_seq $trim_seq_para -trim_detail $detail $fq  && " if (exists $self->{"software:trim_seq"});
				my $prefix=$1 if ($fqname=~/(\S+)\.[^\.\s]+$/);
				if (defined $fastqcut)
				{
					$fltap_cmd .= "$fastqcut $fqname.trim.out $fastqcut_para -prefix $fqname > $prefix.clean.fastq\n";
					$fltap_cmd_multi .= " $fastqcut $fqname.trim.out $fastqcut_para -prefix $fqname > $prefix.clean.fastq && ";
					push @{$Reads{$i}},["$prefix.clean.fastq","$fqname.filter.out $fqname.qcut.out"];
				}
				else
				{
					$fltap_cmd .= "mv $fqname.trim.out $prefix.clean.fastq\n";
					$fltap_cmd_multi .= " mv $fqname.trim.out $prefix.clean.fastq && ";
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
						my $fltfastq2pe=checkPath($self->{"software:fltfastq2pe"});
						my ($pair1,$pair2)=($reads1,$reads2);
						$pair1=~s/fastq$/pair.fastq/i;
						$pair2=~s/fastq$/pair.fastq/i;
						my ($single1,$single2)=($reads1,$reads2);
						$single1=~s/fastq$/single.fastq/i;
						$single2=~s/fastq$/single.fastq/i;
						$fltap_cmd .= "$fltfastq2pe -fastq1 $reads1 -fastq2 $reads2 $filter1 $filter2\n";
						$fltap_cmd .= "rm $reads1 $reads2\n" if (exists $self->{"CustomSetting:Clean"} && ($self->{"CustomSetting:Clean"}=~/y/i || $self->{"CustomSetting:Clean"}=~/TRUE/i) );
						$fltap_cmd_multi .= "$fltfastq2pe -fastq1 $reads1 -fastq2 $reads2 $filter1 $filter2 && ";
						$fltap_cmd_multi .= "rm $reads1 $reads2 && " if (exists $self->{"CustomSetting:Clean"} && ($self->{"CustomSetting:Clean"}=~/y/i || $self->{"CustomSetting:Clean"}=~/TRUE/i) );
						${$self->{"LIB"}{$lib}{"$lb$i"}}[$p]=$self->{"-workdir"}."/".$self->{"CustomSetting:qc_outdir"}."/$lib/$pair1";
						${$self->{"LIB"}{$lib}{"$lb$j"}}[$p]=$self->{"-workdir"}."/".$self->{"CustomSetting:qc_outdir"}."/$lib/$pair2";
						push @{$self->{"LIB"}{$lib}{"fq"}},$self->{"-workdir"}."/".$self->{"CustomSetting:qc_outdir"}."/$lib/$single1";
						foreach my $k(keys %{$self->{$lib}{"$lb$i"}{$p}})
						{
							$self->{$lib}{'fq'}{@{$self->{"LIB"}{$lib}{"fq"}}-1}{$k}=$self->{$lib}{"$lb$i"}{$p}{$k};
						}
						push @{$self->{"LIB"}{$lib}{"fq"}},$self->{"-workdir"}."/".$self->{"CustomSetting:qc_outdir"}."/$lib/$single2";
						foreach my $k(keys %{$self->{$lib}{"$lb$i"}{$p}})
						{
							$self->{$lib}{'fq'}{@{$self->{"LIB"}{$lib}{"fq"}}-1}{$k}=$self->{$lib}{"$lb$i"}{$p}{$k};
						}
						$p++;
					} else {
						${$self->{"LIB"}{$lib}{"$lb$i"}}[$s1]=$self->{"-workdir"}."/".$self->{"CustomSetting:qc_outdir"}."/$lib/$reads1";
						$s1++;
						${$self->{"LIB"}{$lib}{"$lb$j"}}[$s2]=$self->{"-workdir"}."/".$self->{"CustomSetting:qc_outdir"}."/$lib/$reads2";
						$s2++;
					}
				} elsif (defined $reads1 && $reads1 ne "") {
					${$self->{"LIB"}{$lib}{$lbmark}}[$s1]=$self->{"-workdir"}."/".$self->{"CustomSetting:qc_outdir"}."/$lib/$reads1";
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
					${$self->{"LIB"}{$lib}{$lbmark}}[$s2]=$self->{"-workdir"}."/".$self->{"CustomSetting:qc_outdir"}."/$lib/$reads2";
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
		$fltap_cmd .= "rm *.out\n" if (exists $self->{"CustomSetting:Clean"} && ($self->{"CustomSetting:Clean"}=~/y/i || $self->{"CustomSetting:Clean"}=~/TRUE/i) );
		$fltap_cmd_multi .= "rm *.out && " if (exists $self->{"CustomSetting:Clean"} && ($self->{"CustomSetting:Clean"}=~/y/i || $self->{"CustomSetting:Clean"}=~/TRUE/i) );
		$fltap_cmd .= qq(rm align.mat\n) unless (!defined $align_matrix || -f "align.mat");
		$fltap_cmd_multi .= qq(rm align.mat && ) unless (!defined $align_matrix || -f "align.mat");
		$fltap_cmd .= "cd ..\n";
		$fltap_cmd_multi .= "cd .. && ";
		$fltap_cmd_multi .= "cd .. ";
		$multi++;
		if ($multi % $self->{"CustomSetting:multithreads"}!=0){
			$fltap_cmd_multi .= "&\n";
		} else {
			$fltap_cmd_multi .= "\n";
		}
	}
	
	$fltap_cmd .= "cd ..\n";
	if (exists $self->{"CustomSetting:multimode"}) {
		chomp $fltap_cmd_multi;
		$fltap_cmd_multi=~s/\&+$//;
		$fltap_cmd_multi.="\n";
		$fltap_cmd_multi.=print_check_process('scanAP','trim_seq','fastqcut','fltfastq2pe');
		return $fltap_cmd_multi;
	} else {
		return $fltap_cmd;
	}
}

sub runRSeQC ($) {
	my $self = shift;
	if (!exists $self->{"CustomSetting:RSeQCPATH"} || !defined $self->{"CustomSetting:RSeQCPATH"}){
		return "";
	}
	my $rseqc_path = $self->{"CustomSetting:RSeQCPATH"};
	my $refseqbed = $self->{"database:refseq-bed"};
	
	my $rseqc_cmd = qq(echo `date`; echo "run RSeQC"\n);
	$rseqc_cmd .= qq(cd $self->{"-workdir"}\n);
	$rseqc_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$rseqc_cmd .= qq(export PYTHONPATH=$self->{"CustomSetting:PYTHONPATH"}:\$PYTHONPATH\n) if (exists $self->{"CustomSetting:PYTHONPATH"});
	$rseqc_cmd .= qq(export PATH=$rseqc_path:\$PATH\n);
	unless (exists $self->{'cmd'}{'aln'} && $self->{'cmd'}{'aln'}==0) {
		$rseqc_cmd .= $self->runBWA("ref");
		$self->{'cmd'}{'aln'} = 0;
	}
	$self->{'cmd'}{'aln'}=0;
	$rseqc_cmd .= qq(export qc_outdir=$self->{"CustomSetting:qc_outdir"}\n);
	$rseqc_cmd .= qq([[ -d \${qc_outdir} ]] || mkdir \${qc_outdir}\n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}));
	$rseqc_cmd .= qq(cd $self->{"CustomSetting:qc_outdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib (@libraries) {
		$rseqc_cmd .= qq(echo `date`; echo "$lib"\n);
		$rseqc_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
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
	my $reference=checkPath($self->{"database:ref"});
	my $genecode=checkPath($self->{"database:gencode-gtf"});
	my $gcfile=checkPath($self->{"database:gencode-gc"});
	my $para=$self->{'CustomSetting:RNA-SeQC'};
	
	my $rnaseqqc_cmd = qq(echo `date`; echo "run RNA-SeQC"\n);
	$rnaseqqc_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$rnaseqqc_cmd .= qq(export qc_outdir=$self->{"CustomSetting:qc_outdir"}\n);
	$rnaseqqc_cmd .= qq(cd $self->{"-workdir"}\n);
	$rnaseqqc_cmd .= qq([[ -d \${qc_outdir} ]] || mkdir \${qc_outdir}\n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}));
	unless (exists $self->{'cmd'}{'aln'} && $self->{'cmd'}{'aln'}==0) {
		$rnaseqqc_cmd .= $self->runBWA("ref");
		$self->{'cmd'}{'aln'}=0;
	}
	$rnaseqqc_cmd .= qq(cd $self->{"CustomSetting:qc_outdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		$rnaseqqc_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
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
	my $seecerpara = $self->{"CustomSetting:SEECER"};
	my $seecer_cmd = qq(echo `date`; echo "run SEECER"\n);;
	$seecer_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$seecer_cmd .= qq(export seecer=$seecer\n);
	$seecer_cmd .= qq(export jellyfish=$jellyfish\n);
	$seecer_cmd .= qq(export seecerpara="$seecerpara"\n);
	$seecer_cmd .= qq(export qc_outdir=$self->{"CustomSetting:qc_outdir"}\n);
	$seecer_cmd .= qq(cd $self->{"-workdir"}\n);
	$seecer_cmd .= qq([[ -d \${qc_outdir} ]] || mkdir \${qc_outdir}\n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}));
	$seecer_cmd .= qq(cd $self->{"CustomSetting:qc_outdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		$seecer_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
		$seecer_cmd .= "cd $lib\n";
		my %fq=getlibSeq($self->{"LIB"}{$lib});
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

sub runBWA($$) {
	my $self=shift;
	my $ref=shift;
	if (!exists $self->{"software:bwa"} || !defined $self->{"software:bwa"} || !-e $self->{"software:bwa"}){
		return "";
	}
	$ref ||= 'ref';
	my $bwa_cmd;
	$bwa_cmd = qq(echo `date`; echo "run BWA"\n);
	$bwa_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	my $bwa=checkPath($self->{"software:bwa"});
	$bwa_cmd .= qq(export bwa="$bwa"\n);
	my $samtools=checkPath($self->{"software:samtools"});
	$bwa_cmd .= qq(export samtools="$samtools"\n);
	my $reference=checkPath($self->{"database:$ref"});
	$bwa_cmd .= qq(export REFERENCE="$reference"\n);
	my $alnpara=$self->{'CustomSetting:bwaaln'} if (exists $self->{'CustomSetting:bwaaln'});
	my $mempara=$self->{'CustomSetting:bwamem'} if (exists $self->{'CustomSetting:bwamem'});
	if ($alnpara!~/\-t\s+\d+/ && exists $self->{'CustomSetting:multithreads'})
	{
		$alnpara.=" -t ".$self->{'CustomSetting:multithreads'};
	}
	$bwa_cmd .= qq(export alnpara="$alnpara"\n) if (defined $alnpara);
	$bwa_cmd .= qq(export mempara="$mempara"\n) if (defined $mempara);
	my $sampepara="";
	$sampepara=$self->{'CustomSetting:sampe'} if (exists $self->{'CustomSetting:sampe'});
	my $samsepara="";
	$samsepara=$self->{'CustomSetting:samse'} if (exists $self->{'CustomSetting:samse'});
	$bwa_cmd .= qq(export heap="$self->{'CustomSetting:heap'}"\n);
	my ($picard,$FixMateInformation)=("","");
	if (exists $self->{"software:picard"})
	{
		$picard=checkPath($self->{"software:picard"});
		my $picardpath=$1 if ($self->{"software:picard"}=~/(.*)\/[^\/\s]+$/);
		$bwa_cmd .= qq(export picard="$picard"\n);
		$picard = correctJavaCmd($picard,"\${heap}");
		$FixMateInformation=(exists $self->{"software:FixMateInformation"})?$self->{"software:FixMateInformation"}:qq($picardpath/FixMateInformation.jar);
	}
	
	my $workdir=checkPath($self->{"-workdir"});
	$bwa_cmd .= qq(export workdir="$workdir"\n);
	my $alndir=checkPath($self->{"CustomSetting:aln_outdir"});
	$bwa_cmd .= qq(export alndir="$alndir"\n);
	$bwa_cmd .= qq(cd \${workdir}\n);
	$bwa_cmd .= qq([[ -d \${alndir} ]] || mkdir \${alndir}\n) if (!-d qq($self->{"-workdir"}/$alndir));;
	$bwa_cmd .= qq(cd \${alndir}\n);
	$bwa_cmd .= "\${bwa} index -a bwtsw \$REFERENCE\n" unless (checkIndex('bwa',$reference)==1);
	my $db=$1 if ($reference =~ /([^\/\.]+)\.fa/);
	my @libraries=sort keys %{$self->{'LIB'}};
	my $skip = 0;
	foreach my $lib(@libraries) {
		$bwa_cmd .= qq(echo `date`; echo "$lib"\n);
		$bwa_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) if (!-e qq($self->{"-workdir"}/$self->{"CustomSetting:aln_outdir"}/$lib));
		$bwa_cmd .= "cd $lib\n";
		if (exists $self->{$lib}{"$ref-bwabam"}) {
			$skip++;
			$bwa_cmd .= "cd ../\n";
			next;
		}
		my %fq=getlibSeq($self->{"LIB"}{$lib});
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
						$bwa_cmd .= "\${samtools} rmdup $lib.pair.$k.bam $lib.pair.$k.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.pair.$k.rmdup.bam $lib.pair.$k.rmdup.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.pair.$k.rmdup.sort.bam\n";
						$bam="$lib.pair.$k.rmdup.sort.bam";
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
						$bwa_cmd .= "\${samtools} rmdup $lib.pair.bam $lib.pair.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.pair.rmdup.bam $lib.pair.rmdup.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.pair.rmdup.sort.bam\n";
						$bam="$lib.pair.rmdup.sort.bam";
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
						$bwa_cmd .= "\${samtools} rmdup $lib.pair.$k.bam $lib.pair.$k.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.pair.$k.rmdup.bam $lib.pair.$k.rmdup.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.pair.$k.rmdup.sort.bam\n";
						$bam="$lib.pair.$k.rmdup.sort.bam";
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
						$bwa_cmd .= "\${samtools} rmdup $lib.pair.bam $lib.pair.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.pair.rmdup.bam $lib.pair.rmdup.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.pair.rmdup.sort.bam\n";
						$bam="$lib.pair.rmdup.sort.bam";
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
					my $FixMateInformationPara=(exists $self->{"CustomSetting:FixMateInformation"})? $self->{"CustomSetting:FixMateInformation"}:"SO=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT";
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
						my $bamtools=checkPath($self->{"software:bamtools"});
						$bwa_cmd .= qq(export bamtools\="$bamtools"\n);
						my $stats="$1.stats" if ($fxmtbam=~/(\S+)\.bam/);
						$bwa_cmd .= qq(\${bamtools} stats -insert -in $fxmtbam > $stats\n);
					}
					$bwa_cmd .= qq(rm -rf tmp_fixmate\n);
					$bam=$fxmtbam;
				}
				push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$bam";
				$bwa_cmd .= "rm *.sai *.pair.bam *.rmdup.bam\n" if (exists $self->{"CustomSetting:Clean"} && ($self->{"CustomSetting:Clean"}=~/y/i || $self->{"CustomSetting:Clean"}=~/TRUE/i) );
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
						$bwa_cmd .= "\${samtools} rmdup $lib.single.$k.bam $lib.single.$k.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.single.$k.rmdup.bam $lib.single.$k.rmdup.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.single.$k.rmdup.sort.bam\n";
						push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.single.$k.rmdup.sort.bam";
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
						$bwa_cmd .= "\${samtools} rmdup $lib.single.bam $lib.single.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.single.rmdup.bam $lib.single.rmdup.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.single.rmdup.sort.bam\n";
						push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.single.rmdup.sort.bam";
					}
					$bwa_cmd .= "rm *.single.bam *.single.?.bam *.rmdup.bam\n" if (exists $self->{"CustomSetting:Clean"});
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
						$bwa_cmd .= "\${samtools} rmdup $lib.single.$k.bam $lib.single.$k.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000  $lib.single.$k.rmdup.bam $lib.single.$k.rmdup.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.single.$k.rmdup.sort.bam\n";
						push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.single.$k.rmdup.sort.bam";
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
						$bwa_cmd .= "\${samtools} rmdup $lib.single.bam $lib.single.rmdup.bam\n";
						$bwa_cmd .= "\${samtools} sort -m 3000000000 $lib.single.rmdup.bam $lib.single.rmdup.sort\n";
						$bwa_cmd .= "\${samtools} index $lib.single.rmdup.sort.bam\n";
						push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.single.rmdup.sort.bam";
					}
					$bwa_cmd .= "rm *.sai *.single.bam *.single.?.bam *.rmdup.bam\n" if (exists $self->{"CustomSetting:Clean"});
				}
			}
		}
		if (exists $self->{'CustomSetting:bwamerge'} && exists $self->{"software:picard"}) {
			my $MergeSamFiles="$1/MergeSamFiles.jar" if ($self->{"software:picard"}=~/(.*)\/[^\/\s]+$/);
			$MergeSamFiles=qq(java -Xmx\${heap} -Djava.io.tmpdir=./tmp_merge -jar $MergeSamFiles) if ($MergeSamFiles!~/^java/ && $MergeSamFiles !~ /\-jar/);
			my $MergeSamFilesPara=$self->{"CustomSetting:MergeSamFiles"};
			my $merge_bam=join " INPUT=",@{$self->{$lib}{"ref-bwabam"}};
			$bwa_cmd .= "$MergeSamFiles INPUT=$merge_bam $MergeSamFilesPara OUTPUT=$lib.merge.bam\n";
			$bwa_cmd .= "\${samtools} index $lib.merge.bam\n";
			$bwa_cmd .= "\${samtools} rmdup $lib.merge.bam - | \${samtools} rmdup -S - - | \${samtools} sort - $lib.merge.rmdup.sort\n";
			$bwa_cmd .= "\${samtools} index $lib.merge.rmdup.sort.bam\n";
			$bwa_cmd .= "rm -rf ./tmp_merge $lib.merge.bam*\n" if (exists $self->{"CustomSetting:Clean"});
			@{$self->{$lib}{"$ref-bwabam"}}=();
			my $bam=$self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.merge.rmdup.sort.bam";
			push @{$self->{$lib}{"$ref-bwabam"}},$bam;
			if (exists $self->{"overlap"} && exists $self->{"sam2bed"} && exists $self->{"msort"}) {
				$bwa_cmd .= $self->stat_mappedreads("bam",$bam,"lib",$lib);
			}
		}
		if (exists $self->{'CustomSetting:UnmappedRealign'})
		{
			if (exists $self->{$lib}{"$ref-bwabam"} && exists $self->{$lib}{"software:sam2reads"} && exists $self->{$lib}{"software:fltfastq2pe"})
			{
				my $sam2reads=checkPath($self->{$lib}{"software:sam2reads"});
				$bwa_cmd .= qq(export sam2reads="$sam2reads"\n);
				my $fltfastq2pe=checkPath($self->{$lib}{"software:fltfastq2pe"});
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
						push @{$self->{"LIB"}{$lib}{'fq1'}},qq($self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.$i-unmpped.R1.pair.fastq);
						push @{$self->{"LIB"}{$lib}{'fq2'}},qq($self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.$i-unmpped.R2.pair.fastq);
						push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.$i-unmpped.R1.single.fastq);
						push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.$i-unmpped.R2.single.fastq);
						push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.$i-unmpped.R3.fastq);
					}
					else
					{
						$bwa_cmd .= qq(\${samtools} -f 4 $bwabam | \${sam2reads} -R1 $lib.unmpped.R1.fastq -R2 $lib.unmpped.R2.fastq -R3 $lib.unmpped.R3.fastq -f fq\n);
						$bwa_cmd .= qq(\${fltfastq2pe} -fastq1 $lib.unmpped.R1.fastq -fastq2 $lib.unmpped.R2.fastq\n);
						push @{$self->{"LIB"}{$lib}{'fq1'}},qq($self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.unmpped.R1.pair.fastq);
						push @{$self->{"LIB"}{$lib}{'fq2'}},qq($self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.unmpped.R2.pair.fastq);
						push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.unmpped.R1.single.fastq);
						push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.unmpped.R2.single.fastq);
						push @{$self->{"LIB"}{$lib}{'fq'}},qq($self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.unmpped.R3.fastq);
					}
					$i++;
				}
			}
		}
		$bwa_cmd.="cd ..\n";
	}
	$bwa_cmd.="cd ..\n";
	if ($skip<@libraries) {
		return $bwa_cmd;
	} else {
		return "";
	}
}



sub runBowtie($$) {
	my ($self,$ref) = @_;
	if (!exists $self->{"software:bowtie"} || !-e $self->{"software:bowtie"} || !defined $self->{"software:bowtie"}) {
		return "";
	}
	$ref ||= 'ref';
	my $bowtie_cmd = qq(echo `date`; echo "run TopHat"\n);
	$bowtie_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	my $bowtie=checkPath($self->{'software:bowtie'});
	$bowtie_cmd.="export bowtie=$bowtie\n";
	my $samtools = checkPath($self->{"software:samtools"});
	$bowtie_cmd .= qq(export samtools="$samtools"\n);
	my $reference=checkPath($self->{"database:$ref"});
	$bowtie_cmd.="export reference=$reference\n";
	my $bowtiepara=$self->{'CustomSetting:bowtie'};
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
	my $workdir=checkPath($self->{"-workdir"});
	$bowtie_cmd.="export workdir=$workdir\n";
	my $alndir=checkPath($self->{"CustomSetting:aln_outdir"});
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
		my %fq=getlibSeq($self->{"LIB"}{$lib});
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
				my $MergeSamFiles="$1/MergeSamFiles.jar" if ($self->{"software:picard"}=~/(.*)\/[^\/\s]+$/);
				$MergeSamFiles=qq(java -Xmx\${heap} -Djava.io.tmpdir=./tmp_merge -jar $MergeSamFiles) if ($MergeSamFiles!~/^java/ && $MergeSamFiles !~ /\-jar/);
				my $MergeSamFilesPara=$self->{"CustomSetting:MergeSamFiles"};
				my $merge_bam=join " INPUT=",@bam;
				$bowtie_cmd .= "$MergeSamFiles INPUT=$merge_bam $MergeSamFilesPara OUTPUT=$lib.merge.bam\n";
				$bowtie_cmd .= "\${samtools} index $lib.merge.bam\n";
				$bowtie_cmd .= "\${samtools} rmdup $lib.merge.bam - |samtools rmdup -S - - | $samtools sort - $lib.merge.rmdup.sort\n";
				$bowtie_cmd .= "\${samtools} index $lib.merge.rmdup.sort.bam\n";
				$bowtie_cmd .= "rm -rf ./tmp_merge $lib.merge.bam\n" if (exists $self->{"CustomSetting:Clean"});
				@{$self->{$lib}{"$ref-bowtiebam"}}=();
				push @{$self->{$lib}{"$ref-bowtiebam"}},"$workdir/$alndir/$lib/$lib.merge.rmdup.sort.bam";
			} else {
				my $merge_bam=join " ",@bam;
				$bowtie_cmd .= qq(\${samtools} view -H $bam[0] |grep -v "^\@RG" | grep -v "^\@PG" >> $lib.inh.sam\n);
				foreach my $tophatbam(@bam)
				{
					$bowtie_cmd .=  qq(\${samtools} view -H $tophatbam |grep RG >> $lib.inh.sam\n);
				}
				$bowtie_cmd .=  qq(\${samtools} view -H $bam[0] |grep PG >> $lib.inh.sam\n);
				$bowtie_cmd .=  "\${samtools} merge -f -nr -h $lib.inh.sam $lib.merge.bam $merge_bam\n";
				$bowtie_cmd .= "\${samtools} rmdup $lib.merge.bam - |samtools rmdup -S - - | $samtools sort - $lib.merge.rmdup.sort\n";
				$bowtie_cmd .= "\${samtools} index $lib.merge.rmdup.sort.bam\n";
				$bowtie_cmd .= "rm -rf ./tmp_merge $lib.merge.bam\n" if (exists $self->{"CustomSetting:Clean"});
				@{$self->{$lib}{"$ref-bowtiebam"}}=();
				push @{$self->{$lib}{"$ref-bowtiebam"}},"$workdir/$alndir/$lib/$lib.merge.rmdup.sort.bam";
			}
		}
		else
		{
			@{$self->{$lib}{"$ref-bowtiebam"}}=@bam;
		}
	}
	$bowtie_cmd .= "cd ..\n";
	return ($bowtie_cmd);
}

sub runSOAP ($) {
	
}

sub runTopHat($) {
	my $self=shift;
	my $ref=shift;
	if (!exists $self->{"software:tophat"} || !-e $self->{"software:tophat"} || !defined $self->{"software:bowtie"}) {
		return "";
	}
	$ref ||= 'ref';
	my $tophat_cmd = qq(echo `date`; echo "run TopHat"\n);
	$tophat_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	my $tophat=checkPath($self->{'software:tophat'});
	$tophat_cmd.="export tophat=$tophat\n";
	my $bowtie=checkPath($self->{'software:bowtie'});
	$tophat_cmd.="export bowtie=$bowtie\n";
	my $samtools = checkPath($self->{"software:samtools"});
	$tophat_cmd .= qq(export samtools="$samtools"\n);
	my $reference=checkPath($self->{"database:$ref"});
	$tophat_cmd.="export reference=$reference\n";
	#my $refGene=checkPath($self->{'database:refGene'});
	#$tophat_cmd.="export refGene=$refGene\n";
	my $tophatpara=$self->{'CustomSetting:tophat'};
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
		$tophat_cmd .= "\${bowtie_build} \$reference $db_index\n" unless (checkIndex('bowtie2',$db_index)==1);
		$reference=$db_index;
	}
	$tophat_cmd .= qq(export REFERENCE=$reference\n);
	my $workdir=checkPath($self->{"-workdir"});
	$tophat_cmd.="export workdir=$workdir\n";
	$tophat_cmd.=qq(cd \$workdir\n);
	
	my @libraries=sort keys %{$self->{'LIB'}};
	my $tophat_outdir = (exists $self->{"CustomSetting:tophat_outdir"})? $self->{"CustomSetting:tophat_outdir"} : "tophat";
	$tophat_cmd .= "[[ -d $tophat_outdir ]] || mkdir $tophat_outdir\n" if (!-d qq($self->{"-workdir"}/$tophat_outdir));
	$tophat_cmd .= "cd $tophat_outdir\n";
	foreach my $lib(@libraries) {
		$tophat_cmd .= qq(echo `date`; echo "$lib"\n);
		#$tophat_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/tophat/$lib));
		#$tophat_cmd .= "cd $lib\n";
		my %fq=getlibSeq($self->{"LIB"}{$lib});
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
			my $rg="--rg-id $ID--rg-platform $PL --rg-library$LB --rg-sample $SM";
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
			my $rg="--rg-id $ID--rg-platform $PL --rg-library$LB --rg-sample $SM";
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
				my $MergeSamFilesPara=$self->{"CustomSetting:MergeSamFiles"};
				my $merge_bam=join " INPUT=",@bam;
				$tophat_cmd .= "$MergeSamFiles INPUT=$merge_bam $MergeSamFilesPara OUTPUT=$lib.tophat.merge.bam\n";
				$tophat_cmd .= "\${samtools} index $lib.merge.bam\n";
				@{$self->{$lib}{"tophatbam"}}=();
				$self->{$lib}{"tophatbam"}=$self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.tophat.merge.bam";
				if (exists $self->{$lib}{"$ref-bwabam"})
				{
					@{$self->{$lib}{"$ref-bwabam"}}=();
					push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.tophat.merge.bam";
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
				$self->{$lib}{"tophatbam"}=$self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.tophat.merge.bam";
				if (exists $self->{$lib}{"$ref-bwabam"})
				{
					@{$self->{$lib}{"$ref-bwabam"}}=();
					push @{$self->{$lib}{"$ref-bwabam"}},$self->{'-workdir'}."/".$self->{"CustomSetting:aln_outdir"}."/$lib/$lib.tophat.merge.bam";
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

sub runBLAT ($) {
	
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

sub runGATK ($$) {
	my $self=shift;
	my $ref=shift;
	$ref ||= 'ref';
	my $callVar_cmd = qq(echo `date`; echo "run Samtools-Picard-GATK"\n);
	$callVar_cmd .= qq(export PATH="$self->{"CustomSetting:PATH"}":\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	my $heap=$self->{"CustomSetting:heap"};
	$callVar_cmd .= qq(export heap="$heap"\n);
	my $samtools = checkPath($self->{"software:samtools"});
	$callVar_cmd .= qq(export samtools="$samtools"\n);
	my $samtools_path=$1 if ($samtools =~ /(\S+)\/samtools/);
	my ($picard,$picardpath,$MergeSamFiles,$MarkDuplicates,$bamtools);
	if (defined $samtools_path)
	{
		$callVar_cmd .= qq(export samtools_path="$samtools_path"\n);
		$callVar_cmd .= qq (export PATH="\${samtools}_path":"\${samtools}_path/bcftools":\$PATH\n);
	}
	if (exists $self->{"software:picard"})
	{
		$picard=$self->{"software:picard"};
		$picardpath=$1 if ($self->{"software:picard"}=~/(.*)\/[^\/\s]+$/);
		$callVar_cmd .= qq(export picard="$picard"\n);
		$picard = correctJavaCmd($picard,"\${heap}");
		$MergeSamFiles=(exists $self->{"software:MergeSamFiles"})?$self->{"software:MergeSamFiles"}:qq($picardpath/MergeSamFiles.jar);
		$callVar_cmd .= qq(export MergeSamFiles="$MergeSamFiles"\n);
		$MergeSamFiles=correctJavaCmd($MergeSamFiles,"\${heap}","./tmp_merge");
		$MarkDuplicates=(exists $self->{"software:MarkDuplicates"})?$self->{"software:MarkDuplicates"}:"$picardpath/MarkDuplicates.jar";
		$callVar_cmd .= qq(export MarkDuplicates="$MarkDuplicates"\n);
		$MarkDuplicates=correctJavaCmd($MarkDuplicates,"\${heap}","./tmp_rmdup");
	}
	if (defined $self->{"software:bamtools"})
	{
		$bamtools=checkPath($self->{"software:bamtools"});
		$callVar_cmd .= qq(export bamtools\="$bamtools"\n);
	}
	my $gatk=checkPath($self->{"software:gatk"}) if (exists $self->{"software:gatk"});
	$gatk=correctJavaCmd($gatk,$heap,"./tmp_gatk");
	$callVar_cmd .= qq(export gatk="$gatk"\n);
	my $reference=checkPath($self->{"database:$ref"});
	$callVar_cmd .= qq(export REFERENCE="$reference"\n);
	if (exists $self->{"database:dbSNP"})
	{
		$callVar_cmd .=qq(export dbSNP="$self->{"database:dbSNP"}"\n);
	}
	my $workdir=checkPath($self->{"-workdir"});
	my $vardir=checkPath($self->{"CustomSetting:var_outdir"});
	$callVar_cmd .= qq(export workdir="$workdir"\n);
	$callVar_cmd .= qq(export vardir="$vardir"\n);
	$callVar_cmd .= qq(cd \${workdir}\n);
	$callVar_cmd .= qq([[ -d \${vardir} ]] || mkdir -p \${vardir}\n) if (!-d qq($self->{"-workdir"}/$vardir));
	$callVar_cmd .= qq(cd \${vardir}\n);
	$callVar_cmd .= "\${samtools} index -a bwtsw \$REFERENCE\n" unless (checkIndex('samtools',$reference)==1);
	my @libraries=sort keys %{$self->{'LIB'}};
	my $multi=0;
	foreach my $lib(@libraries) {
		my $bam="";
		$callVar_cmd .= qq(echo `date`; echo "$lib"\n);
		if (!-e qq($self->{"-workdir"}/$self->{"CustomSetting:var_outdir"}/$lib))
		{
			$callVar_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) ;
		}
		$callVar_cmd .= "cd $lib\n";
		if (!exists $self->{$lib}{"$ref-bam"}) {
			if (exists $self->{$lib}{"$ref-bwabam"} && @{$self->{$lib}{"$ref-bwabam"}}>1) {
				my $merge_bam="";
				if (defined $MergeSamFiles && $MergeSamFiles ne "") {
## Merge BAM 
					#$callVar_cmd .= qq(mkdir -p ./tmp_merge);
					#$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
					#$callVar_cmd .= qq(export tmp_merge="./tmpmerge");
					#$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
					my $MergeSamFilesPara=(exists $self->{"CustomSetting:MergeSamFiles"})?$self->{"CustomSetting:MergeSamFiles"}:'USE_THREADING=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT';
					$merge_bam=join " INPUT\=",@{$self->{$lib}{"ref-bwabam"}};
					$callVar_cmd .= "$MergeSamFiles INPUT\=$merge_bam $MergeSamFilesPara OUTPUT\=$lib.merge.bam";
					$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
					if (exists $self->{"CustomSetting:Clean"} && ($self->{"CustomSetting:Clean"}=~/y/i || $self->{"CustomSetting:Clean"}=~/TRUE/i) )
					{
						$callVar_cmd .= qq(rm -rf ./tmp_merge);
						$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
					}
					$bam="$lib.merge.bam";
				} else {
					$merge_bam=join " ",@{$self->{$lib}{"ref-bwabam"}};
					$callVar_cmd .=  qq(\${samtools} view -H ${$self->{$lib}{"ref-bwabam"}}[0] |grep -v "^\@RG" | grep -v "^\@PG" >> $lib.inh.sam);
					$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
					foreach my $bwabam(@{$self->{$lib}{"ref-bwabam"}})
					{
						$callVar_cmd .=  qq(\${samtools} view -H $bwabam |grep "^\@RG" >> $lib.inh.sam);
						$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
					}
					$callVar_cmd .=  qq(\${samtools} view -H ${$self->{$lib}{"ref-bwabam"}}[0] |grep "^\@PG" >> $lib.inh.sam);
					$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
					$callVar_cmd .=  "\${samtools}merge -f -nr -h $lib.inh.sam $lib.merge.bam $merge_bam";
					$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
					$bam="$lib.merge.bam";
				}
				$self->{$lib}{"$ref-bam"}=qq($self->{"-workdir"}/$self->{"CustomSetting:var_outdir"}/$lib/$lib.merge.bam);
			}
			elsif (exists $self->{$lib}{"$ref-bwabam"}) {
				$bam=${$self->{$lib}{"$ref-bwabam"}}[0];
				#$callVar_cmd .= "\${samtools} mpileup -ugf \$REFERENCE $bam | bcftools view -bvcg - | bcftools view -cg - > $lib.var.vcf && vcfutils.pl varFilter -D100 $lib.var.vcf > $lib.var.flt.vcf";
				#$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
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
			$callVar_cmd .= qq(\${bamtools} filter -isMapped true -isPaired true -in $bam -out $fltbam);
			$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
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
		if (defined $MarkDuplicates && $MarkDuplicates ne "") {
			#$callVar_cmd .= qq(mkdir -p ./tmp_rmdup);
			#$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			#$callVar_cmd .= qq(export tmp_rmdup="./tmp_rmdup");
			#$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			my $MarkDuplicatesPara=(exists $self->{"CustomSetting:MarkDuplicates"})?$self->{"CustomSetting:MarkDuplicates"}:'VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true ASSUME_SORTED=true';
			my $rmdupbam="$1.rmdup.bam" if ($bam=~/([^\/\s]+)\.bam/);
			$callVar_cmd .= qq($MarkDuplicates INPUT=$bam OUTPUT=$rmdupbam M=$lib.duplicate_report.txt $MarkDuplicatesPara);
			$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			$callVar_cmd .= qq(\${samtools} index $rmdupbam);
			$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			$bam=$rmdupbam;
			if (exists $self->{"CustomSetting:Clean"} && ($self->{"CustomSetting:Clean"}=~/y/i || $self->{"CustomSetting:Clean"}=~/TRUE/i) )
			{
				$callVar_cmd .= qq(rm -rf ./tmp_rmdup);
				$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			}
		}else{
			my $rmdupbam="$1.rmdup.sort" if ($bam=~/([^\/\s]+)\.bam/);
			$callVar_cmd .= "\${samtools} rmdup $lib.merge.bam - | \${samtools} rmdup -S - - | \${samtools} sort -m 3000000000 - $rmdupbam";
			$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			$callVar_cmd .= "\${samtools} index $rmdupbam.bam";
			$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			$bam="$rmdupbam.bam";
		}
## Stat
#${bamtools} stats \
#  -insert \
#  -in $PWDS/${subjectID}.rmdup.bam \
#  > $PWDS/${subjectID}.rmdup.stats

		if (defined $self->{"software:bamtools"}) {
				my $stats="$1.stats" if ($bam=~/([^\/\s]+)\.bam/);
				$callVar_cmd .= qq(\${bamtools} stats -insert -in $bam > $stats\n);
				$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
		} 

		if (defined $gatk) {
##Local Realignment
#------------------
#Samtools calls short indels with local realignment, but it does not write a modified BAM file after the realignment.
#The GATK though provides such a tool that realigns reads in regions with suspected indel artifacts and generates a BAM with cleaned alignments.
			#$callVar_cmd .= "mkdir ./tmp_realign";
			#$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			$callVar_cmd .= qq($gatk -T RealignerTargetCreator -R \$REFERENCE -I $bam -o $lib.gatk.intervals);
			if (exists $self->{"database:dbSNP"})
			{
				$callVar_cmd .= qq( -know \$dbSNP);
			}
			$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? qq( -nt $self->{"CustomSetting:multithreads"} && ) : "\n";
			my $realignedbam="$1.realigned.bam" if ($bam=~/([^\/\s]+)\.bam/);
			$callVar_cmd .= qq($gatk -T IndelRealigner -R \$REFERENCE -I $bam -o $realignedbam -targetIntervals $lib.gatk.intervals -LOD 0.4 -compress 6 -l INFO);
			if (exists $self->{"database:dbSNP"})
			{
				$callVar_cmd .= qq( -know \$dbSNP);
			}
			$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? qq( && ) : "\n";
			if (defined $self->{"software:bamtools"})
			{
				my $stats="$1.stats" if ($realignedbam=~/([^\/\s]+)\.bam/);
				$callVar_cmd .= qq(\${bamtools} stats -insert -in $realignedbam  > $stats\n);
				$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
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
			#$callVar_cmd .= "mkdir ./tmp_covar";
			#$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			if (exists $self->{"database:dbSNP"})
			{
				$callVar_cmd .= qq($gatk -T BaseRecalibrator -R \$REFERENCE -knowSites \$dbSNP -l INFO -I $bam -o $lib.recalibration_report.grp -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate);
				$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? qq( && ) : "\n";

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
					$AnalyzeCovariates=checkPath("$1/resource/AnalyzeCovariates.jar") if ($self->{"software:gatk"}=~/(\S+)\/[^\/\s]+$/);
					if (defined $AnalyzeCovariates && -f $AnalyzeCovariates)
					{
						$AnalyzeCovariates=correctJavaCmd($AnalyzeCovariates,$heap,"./tmp_covar");
					}
				}
				if (defined $AnalyzeCovariates)
				{
					$callVar_cmd .= qq(mkdir analyzeCovar_v1);
					$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
					$callVar_cmd .= qq($AnalyzeCovariates -recalFile $lib.flt.recal_v1.csv -outputDrir analyzeCovar_v1 -ignoreQ 3);
					$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
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
			#$callVar_cmd .= qq($gatk -I $bam -R \$REFERENCE  -T BaseRecalibrator -l INFO \\\n);
			#$callVar_cmd .= qq(-recalFile recal_data.csv -cov ReadGroupCovariate -cov QualityScoreCovariate \\\n);
			#$callVar_cmd .= qq(-cov CycleCovariate -cov DinucCovariate -cov TileCovariate\n);
			#$callVar_cmd .= qq($gatk -I $bam -R \$REFERENCE -T TableRecalibration -l INFO \\\n);
			#$callVar_cmd .= qq(-recalFile recal_data.csv --output_bam $lib.merge.sort.alnRecal.bam\n);
			#$callVar_cmd .= "$samtools index $lib.merge.sort.alnRecal.bam\n";
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
				$callVar_cmd .= qq($gatk -T PrintReads -R \$REFERENCE -I $bam -BQSR $lib.recalibration_report.grp -o $recalbam && $samtools index $lib.recal.bam);
				$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? qq( && ) : "\n";
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
				$callVar_cmd .= qq($gatk -T BaseRecalibrator -R \$REFERENCE -knowSites \$dbSNP -I $bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -o $lib.flt.recal_v2.csv);
				$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? qq( -nt $self->{"CustomSetting:multithreads"} && ) : "\n";
			}
			#else
			#{
			#	$callVar_cmd .= qq($gatk -T BaseRecalibrator -I $bam -R \$REFERENCE -run_without_dbsnp_potentially_ruining_quality -o $lib.recalibration_report.grp --intermediate_csv_file $lib.recal.csv --plot_pdf_file $lib.comp.pdf);
			#	$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			#	my $recalbam = "$1.recal.bam" if ($bam=~/([^\/\s]+)\.bam$/);
			#	$callVar_cmd .= qq($gatk -T PrintReads -R \$REFERENCE -I $bam -BQSR $lib.recalibration_report.grp -o $recalbam && $samtools index $lib.recal.bam);
			#	$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? qq( && ) : "\n";
			#	$bam = $recalbam;
			#}

## BAQ calmd
			my $baqbam = "$1.baq.bam" if ($bam=~/([^\/\s]+)\.bam$/);
			$callVar_cmd .= "\${samtools} calmd -Abr $bam \$REFERENCE > $baqbam && \${samtools} index $baqbam";
			$bam = $baqbam;
			$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";

## Genotyping calling
			$callVar_cmd .= "$gatk -T UnifiedGenotyper -R \$REFERENCE -I $bam -baq CALCULATE_AS_NECESSARY -o $lib.gatk.var.vcf -U -S SILENT";
			$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? qq( -nt $self->{"CustomSetting:multithreads"} -nct $self->{"CustomSetting:multithreads"} && ) : "\n";
			if (exists $self->{"software:tabix"} || $self->{"software:bgzip"})
			{
				my $bgzip;
				if (exists $self->{"software:bgzip"})
				{
					$bgzip=checkPath($self->{"software:bgzip"});
				}
				else
				{
					my $tabix=checkPath($self->{"software:tabix"}) if (exists $self->{"software:tabix"});
					$bgzip=$tabix;
					$bgzip=~s/tabix$/bgzip/;
				}
				$callVar_cmd .= "$bgzip $lib.gatk.var.vcf";
				$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? qq( && ) : "\n";
			}

##coverage
#-------------------
#The GATK provides a tool to assess the depth of coverage over any number of intervals in a BAM file.
#The coverage will be displayed in the output for every single position plus as an average over each interval specified.
			if (exists $self->{"CustomSetting:TargetIntervalList"})
			{
				my $TL=checkPaht($self->{"CustomSetting:TargetIntervalList"});
				$callVar_cmd .= "$gatk -T DepthOfCoverage -I $bam -R \$REFERENCE -L $TL --minMappingQuality 20 -omitLocusTable $lib.Locus.coverage.txt";
				$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			}
			#$bam=$self->{'-workdir'}."/".$self->{"CustomSetting:var_outdir"}."/$lib/$lib.realigned.baq.bam";
			$self->{$lib}{"$ref-bam"}=$self->{'-workdir'}."/".$self->{"CustomSetting:var_outdir"}."/$lib/$bam";
			$self->{$lib}{"$ref-vcf"}=$self->{'-workdir'}."/".$self->{"CustomSetting:var_outdir"}."/$lib/$lib.gatk.var.vcf";
			if (exists $self->{"CustomSetting:Clean"} && ($self->{"CustomSetting:Clean"}=~/y/i || $self->{"CustomSetting:Clean"}=~/TRUE/i) )
			{
				$callVar_cmd .= "rm -rf $lib.merge.sort.bam $lib.realigned.bam tmp_realign tmp_gatk";
				$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			}
		} else {
			$callVar_cmd .= "\${samtools} mpileup -ugf \$REFERENCE $bam | bcftools view -bvcg - | bcftools view -cg - > $lib.var.vcf &&  vcfutils.pl varFilter -D100 $lib.var.vcf > $lib.var.flt.vcf";
			$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			$self->{$lib}{"$ref-bam"}=$self->{'-workdir'}."/".$self->{"CustomSetting:var_outdir"}."/$lib/$bam";
		}
		if (exists $self->{"software:sam2reads"}) {
			my $name=$1 if ($bam=~/([^\/\s]+).bam/);
			my $sam2reads=checkPath($self->{"software:sam2reads"});
			$callVar_cmd .= "\${samtools} view -F 4 $bam | $sam2reads -R1 $name\_mapped.R1.fastq -R2 $name\_mapped.R2.fastq -R3 $name\_mapped.R3.fastq -f fq";
			$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
			$callVar_cmd .= "\${samtools} view -f 4 $bam | $sam2reads -R1 $name\_unmapped.R1.fastq -R2 $name\_unmapped.R2.fastq -R3 $name\_unmapped.R3.fastq -f fq";
			$callVar_cmd .= (exists $self->{"CustomSetting:multimode"}) ? " && " : "\n";
		}
		if (exists $self->{"overlap"} && exists $self->{"sam2bed"} && exists $self->{"msort"}) {
			$callVar_cmd .= $self->stat_mappedreads("bam",$bam,"lib",$lib);
		}
		$multi++;
		$callVar_cmd  =~ s/\&*\s*$//;
		$callVar_cmd .= ((exists $self->{"CustomSetting:multithreads"}) && ($multi % $self->{"CustomSetting:multithreads"}!=0) && ($multi<@libraries)) ? " &\n" : "\n";
		$callVar_cmd .= "cd ..\n";
	}
	$callVar_cmd  =~ s/\&*\s*$/\n/;
	$callVar_cmd .= "cd ..\n";
	return $callVar_cmd;
}



#/usr/local/bin/dindel --analysis getCIGARindels --bamFile realigned.baq.bam --outputFile dindel_output --ref reference.fa
#python /usr/local/bin/makeWindows.py --inputVarFile dindel_output.variants.txt --windowFilePrefix realign_windows --numWindowsPerFile 1000
#perl -e 'my @f=glob("realign_windows.*.txt");foreach (@f){my $prefix="dindel_stage2_output_windows.$1" if ($_=~/windows\.(\d+)\./);system "/usr/local/bin/dindel --analysis indels --doDiploid --bamFile realigned.baq.bam --ref reference.fa --varFile $_ --libFile dindel_output.libraries.txt --outputFile $prefix";}'
#ls dindel_stage2_output_windows.*.glf.txt > dindel_stage2_outputfiles.txt
#python /usr/local/bin/mergeOutputDiploid.py --inputFiles dindel_stage2_outputfiles.txt --outputFile variantCalls.VCF --ref reference.fa 
sub runDindle ($) {
	my $self=shift;
	my $ref=shift;
	$ref ||= 'ref';
	if (!exists $self->{"software:dindle"}) {
		return "";
	}
	my $dindle=checkPath($self->{"software:dindle"});
	my $dindle_cmd = qq(echo `date`; echo "run Dindle"\n);
	$dindle_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$dindle_cmd .= qq(export dindle=$dindle\n);
	my $para = $self->{"CustomSetting:dindle"};
	$dindle_cmd .= qq(export para="$para"\n);
	my $multithreads = $self->{"CustomSetting:multithreads"};
	my $multirun = $self->{"CustomSetting:multithreads-run"};
	$dindle_cmd .= qq(export multirun=$multirun\n);
	my $reference = $self->{"database:$ref"};
	$dindle_cmd .= qq(export REFERENCE=$reference\n);
	my $dindle_output=(exists $self->{"CustomSetting:dindle_outdir"})?$self->{"CustomSetting:dindle_outdir"}:"dindle";
	$dindle_cmd .= "[[ -d $dindle_output ]] || mkdir $dindle_output\n";
	$dindle_cmd .= "cd $dindle_output\n";
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		if (exists $self->{$lib}{"$ref-bam"}) {
			$dindle_cmd .= "[[ -d $lib ]] || mkdir $lib\n";
			$dindle_cmd .= "cd $lib\n";
			my $bam=$self->{$lib}{"$ref-bam"};
			$dindle_cmd .= qq(\${dindle} --anlysis getCIGARindels --bamFile $bam --outputFile $lib.dindel_output --ref \${REFERENCE}\n);
			$dindle_cmd .= qq(makeWindows.py --inputVarFile $lib.dindel_output.variants.txt --windowFilePrefix $lib.realign_windows --numWindowsPerFile 1000\n);
			#/usr/local/bin/dindel --analysis indels --doDiploid --bamFile realigned.baq.bam --ref reference.fa --varFile $_ --libFile dindel_output.libraries.txt --outputFile $prefix
			$dindle_cmd .= qq(perl -e \'my \@f=glob\("$lib.realign_windows.*.txt"\););
			$dindle_cmd .= qq(for\(my \$i=0;\$i<\@f;\$i+=$multithreads\){my \@cmdary=();foreach my \$j\(\$i..\(\$i+$multithreads-1\)\){my \$prefix="$lib.dindel_stage2_output_windows.\$1" if \(\$f[\$i]=~/windows\\.\(\\d+\)\\./\););
			$dindle_cmd .= qq(push \@cmdary,qq\(\"$dindle --analysis indels --doDiploid --bamFile $bam --ref $reference --varFile \$f[\$i] --libFile $lib.dindel_output.libraries.txt --outputFile \$prefix\"\);}my \$cmd="$multirun ".join " ",\@cmdary;system \$cmd}\'\n);
			$dindle_cmd .= qq(ls $lib.dindel_stage2_output_windows.*.glf.txt > $lib.dindel_stage2_outputfiles.txt\n);
			$dindle_cmd .= qq(mergeOutputDiploid.py --inputFiles $lib.dindel_stage2_outputfiles.txt --outputFile $lib.variantCalls.VCF --ref \${REFERENCE}\n);
			$dindle_cmd .= qq(rm $lib.realign_windows.*.txt $lib.dindel_stage2_output_windows.*.glf.txt\n) if (exists $self->{"CustomSetting:Clean"} && ($self->{"CustomSetting:Clean"}=~/y/i || $self->{"CustomSetting:Clean"}=~/TRUE/i) );
			$dindle_cmd .= qq(cd ..\n);
		}
	}
	return $dindle_cmd;
}

#soapsnp -B aln.bam -d reference.fa -o cns -r 0.0005 -e 0.001 -u -L 150 -2 -Q J -s dbSNP -T region.out -m
#The dbSNP file consist of a lot of lines like this one:
#	chr1    201979756       1       1       0       0.161   0       0       0.839   rs568
#	The columns from left to right are: name of chromosome, coordinate on the chromosome, whether 
#	the SNP	has allele frequency information (1 is true, 0 is false), whether the SNP is validated 
#	by experiment (1 is true, 0 is false), whether the SNP is actually an indel (1 is true, 0 is false),
#	frequency of A, frequency of C, frequency of T, frequency of G, SNP id. For known SNP sites that do
#	not have allele frequency information, the frequency information can be arbitrarily determined as 
#	any positive values, which only imply what alleles have already been deposited in the database.
sub runSOAPsnp($) {
	
}

sub runSOAPsnv ($) {
	
}

sub runSOAPsv ($) {
	
}

sub runBreakDancer ($) {
	
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
	$cufflinks_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	my $cufflinks=checkPath($self->{'software:cufflinks'});
	$cufflinks_cmd .= qq(export cufflinks=$cufflinks\n);
	my $gffread=checkPath($self->{'software:gffread'}) if (exists $self->{'software:gffread'});
	$cufflinks_cmd .= qq(export gffread=$gffread\n) if (defined $gffread);
	my $bowtie=checkPath($self->{'software:bowtie'});
	$cufflinks_cmd .= qq(export bowtie=$bowtie\n);
	my $reference=checkPath($self->{'database:ref'});
	$cufflinks_cmd .= qq(export reference=$reference\n);
	my $refGene=checkPath($self->{'database:refGene'});
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
	my $para=$self->{'CustomSetting:cufflinks'};
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
	my $reference=checkPath($self->{"database:$ref"});
	my $para=$self->{'CustomSetting:cuffmerge'};
	my $cuffmerge="";
	if (exists $self->{'software:cuffmerge'}) {
		$cuffmerge=checkPath($self->{'software:cuffmerge'});
	} else {
		$cuffmerge=checkPath($self->{'software:cufflinks'});
		$cuffmerge=~s/cufflinks$/cuffmerge/;
	}
	my $cuffmerge_cmd = qq(echo `date`; echo "run Cuffmerge"\n);
	$cuffmerge_cmd .= qq(cd $self->{"-workdir"}\n);
	$cuffmerge_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$cuffmerge_cmd .= qq(export REFERENCE=$reference\n);
	$cuffmerge_cmd .= "[[ -d cuffmerge ]] || mkdir cuffmerge\n" if (!-d qq($self->{"-workdir"}/cuffmerge));
	$cuffmerge_cmd .= "cd cuffmerge\n";
	if (exists $self->{'database:$gene'}) {
		my $refGene=checkPath($self->{'database:$gene'});
		my $all_transcripts=join "\\n",@{$self->{'transcript-gtf'}};
		$cuffmerge_cmd.=qq(perl -e 'print "$all_transcripts\\n"' > GTF_list.txt\n);
		$cuffmerge_cmd.="$cuffmerge $para -g $refGene -s \$REFERENCE -p 10 GTF_list.txt\n";
	} else {
		my $all_transcripts=join "\\n",@{$self->{'transcript-gtf'}};
		$cuffmerge_cmd.=qq(perl -e 'print "$all_transcripts\\n"' > GTF_list.txt\n);
		$cuffmerge_cmd.=qq($cuffmerge $para -s \$REFERENCE -p $self->{"CustomSetting:multithreads"} GTF_list.txt\n);
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
				my $refGene=checkPath($self->{'database:$gene'});
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

#cuffcompare -i input_gtf_list -r ../reference/genes.gtf -o compare_out 2> cuffcompare.log &
#cuffcompare -s ../05.db/w14_v7.2.fasta -p 4 assembly_GTF_list.txt
sub runCuffCompare($) {
	my $self=shift;

	if (!exists $self->{"software:cuffcompare"} || !defined $self->{"software:cuffcompare"}) {
		return "";
	}
	my $cuffcompare_cmd = qq(echo `date`; echo "run Cuffcompare"\n);;
	$cuffcompare_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	my $reference=checkPath($self->{'database:ref'});
	$cuffcompare_cmd .= qq(export REFERENCE\="$reference"\n);
	if (exists $self->{'database:$gene'}) {
		my $refGene=$self->{'database:$gene'};
		$cuffcompare_cmd.=qq(export refGene="$refGene"\n);
	}
	my $para=$self->{'CustomSetting:cuffcompare'};
	$cuffcompare_cmd .= qq(export para="$para"\n);
	my $cuffcompare="";
	if (exists $self->{'software:cuffcompare'}) {
		$cuffcompare=checkPath($self->{'software:cuffcompare'});
	} else {
		$cuffcompare=checkPath($self->{'software:cufflinks'});
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


#Trinity.pl --seqType fq --JM 100G --left reads_1.fq  --right reads_2.fq --CPU 6
#TRINITY_RNASEQ_ROOT/util/alignReads.pl --left left.fq --right right.fq --seqType fq --target Trinity.fasta --aligner bowtie
#TRINITY_RNASEQ_ROOT/util/RSEM_util/run_RSEM.pl --transcripts Trinity.fasta --name_sorted_bam bowtie_out.nameSorted.sam.+.sam.PropMapPairsForRSEM.bam --paired
sub runTrinity($) {
	my $self=shift;
	my $trinity=checkPath($self->{"software:trinity"});
	my $trinity_cmd = qq(echo `date`; echo "run Trinity"\n);
	$trinity_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	if (exists $self->{'CustomSetting:TRINITY_RNASEQ_ROOT'}) {
		$trinity_cmd .= qq(export TRINITY_RNASEQ_ROOT="$self->{'CustomSetting:TRINITY_RNASEQ_ROOT'}"\n);
	} else {
		my $trinity_root=$1 if ($trinity=~/(\S+)\/[^\/]$/);
		$trinity_cmd .= qq(export TRINITY_RNASEQ_ROOT="$trinity_root"\n);
	}
	my $para=(exists $self->{"CustomSetting:trinity"})?$self->{"CustomSetting:trinity"}:qq(--seqType fq --JM 100G --CPU $self->{"CustomSetting:multithreads"});
	$trinity_cmd .= qq(cd $self->{"-workdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		$trinity_cmd .= qq(echo `date`; echo "$lib"\n);
		$trinity_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$lib));
		$trinity_cmd .= qq(cd $lib\n);
		$trinity_cmd .= qq([[ -d trinity_asm ]] || mkdir trinity_asm\n) unless (-d qq($self->{"-workdir"}/$lib/trinity_asm));
		$trinity_cmd .= qq(cd trinity_asm\n);
		my %read=getlibSeq($self->{"LIB"}{$lib});
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
	my $velveth=checkPath($self->{"software:velveth"});
	my $velvet_cmd = qq(echo `date`; echo "run Velvet-Oases"\n);
	$velvet_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$velvet_cmd .= qq(export velvet="$velveth"\n);
	my $velvetg="";
	if (exists $self->{"software:velvetg"})
	{
		$velvetg=checkPath($self->{"software:velvetg"});
		$velvet_cmd .= qq(export velveg="$velvetg"\n);
	}
	else
	{
		$velvetg=$self->{"software:velveth"};
		$velvetg=~s/velveth$/velvetg/;
		$velvet_cmd .= qq(export velveg="$velvetg"\n);
	}
	my $oases=checkPath($self->{"software:oascs"}) if (exists $self->{"software:oascs"});
	$velvet_cmd .= qq(export oases="$oases"\n) if (defined $oases);
	my $velvethpara = $self->{'CustomSetting:velveth'};
	$velvet_cmd .= qq(export velvethpara="$velvethpara"\n);
	my $velvetgpara = $self->{'CustomSetting:velvetg'};
	$velvet_cmd .= qq(export velvetgpara="$velvetgpara"\n);
	my $oasespara = $self->{'CustomSetting:oases'} if (exists $self->{'CustomSetting:oases'});
	$velvet_cmd .= qq(export oasespara="$oasespara"\n) if (defined $oasespara);
	my @libraries=sort keys %{$self->{'LIB'}};
	$velvet_cmd .= qq(cd $self->{"-workdir"}\n);
	$velvet_cmd .= qq([[ -d velvet ]] || mkdir velvet\n);
	foreach my $lib(@libraries) {
		$velvet_cmd .= qq(echo `date`; echo "$lib"\n);
		$velvet_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$lib));
		$velvet_cmd .= qq(cd $lib\n);
		my %read=getlibSeq($self->{"LIB"}{$lib});
		my $input="";
		if (exists $read{0}) {
			for (my $i=0;$i<@{$read{0}};$i++) {
				$input.="-short --".$self->check_fileformat(${$read{0}}[$i])." ${$read{0}}[$i]";
			}
		}
		if (exists $read{1} && exists $read{2} && @{$read{2}}==@{$read{1}}) {
			for (my $i=0;$i<@{$read{2}};$i++) {
				my $reads1=${$read{1}}[$i];
				my $reads2=${$read{2}}[$i];
				my $format=$self->check_fileformat($reads1);
				my $reads=shuffleSequences($reads1,$reads2,$format);
				$input.="-shortPaired --".$self->check_fileformat($reads)." $reads";
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

#SOAPdenovo-63mer all -s WGS_soapdenovo.config -K 25 -F -p 32 -a 92 -m 51 -M 1 -E -k 31 -F -L 100 -V -o WGS_25mer 1> 25mer.log 2> 25mer.log
sub runSOAPdenovo ($) {
	my $self=shift;
	if (!exists $self->{"software:soapdenovo"})
	{
		return "";
	}
	my $soapdenovo=checkPath($self->{"software:soapdenovo"}) if (exists $self->{"software:soapdenovo"});
	my $soapdenovo_cmd=qq(echo `date`; echo "run SOAPdenvo|SOAPdenovo-Trans"\n);
	$soapdenovo_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$soapdenovo_cmd .= qq(export soapdenvo="$soapdenovo"\n);
	my $config=checkPath($self->{"CustomSetting:soapdenovo_config"}) if (exists $self->{"CustomSetting:soapdenovo_config"});
	my $para=$self->{"CustomSetting:soapdenovo"};
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
		my %read=getlibSeq($self->{"LIB"}{$lib});
		if (exists $read{1} && exists $read{2}) {
			print OUT "[LIB]\n";
			for (my $i=0;$i<@{$read{2}};$i++) {
				my $type=(check_fileformat(${$read{1}}[$i])=~/fastq/i)?"fq":"fa";
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
				my $head=(check_fileformat(${$read{0}}[$i])=~/fastq/i)?"q":"f";
				print OUT qq($head=${$read{0}}[$i]\n);
			}
		}
		
		close OUT;
		return $config;
	}
}

#phrap seq.fas -new_ace -revise_greedy -shatter_greedy -forcelevel 0 -repeat_stringency 0.95 > phrap.out
sub runPhrap ($) {
	
}

sub runCAP3 ($) {
	
}

sub runPCAP ($) {
	
}


#/opt/454/bin/newAssembly BAC1
#/opt/454/bin/addRun BAC1 ~/01.data/W14/454/BAC/5_BAC_454/BAC_1.sff 
#/opt/454/bin/runProject -cpu 8 BAC1
sub runNewbler ($) {
	my $self = shift;
	if (!exists $self->{"CustomSetting::Newbler"} || !-d $self->{"CustomSetting::Newbler"}) {
		return "";
	}
	my $newbler_cmd = qq(echo `date`; echo "run Newbler"\n);
	$newbler_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$newbler_cmd .= qq(export NEWBLER="$self->{"CustomSetting::Newbler"}"\n);
	my $para = $self->{'CustomSetting:Newbler'} if (exists $self->{'CustomSetting:Newbler'});
	$newbler_cmd .= qq(export para=$para) if (defined $para);
	my $workdir .= checkPath($self->{"-workdir"});
	$newbler_cmd .= qq(export workdir=$workdir\n);
	$newbler_cmd .= qq(cd \${workdir}\n);
	$newbler_cmd .= "[[ -d newbler] || mkdir newbler\n" if (!-d qq($self->{"-workdir"}/newbler));
	$newbler_cmd .= qq(cd newbler\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries)
	{
		my %read=getlibSeq($self->{"LIB"}{$lib});
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
sub runMSR_CA($) {
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
	$cuffdiff_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	my $reference=checkPath($self->{"database:$ref"});
	$cuffdiff_cmd .= qq(export reference=$reference\n);
	my $refGene=(exists $self->{"database:$gene"})?checkPath($self->{"database:$gene"}):$self->{$gene};
	$cuffdiff_cmd .= qq(export refGene=$refGene\n);
	my $cuffdiff="";
	if (exists $self->{'software:cuffdiff'}) {
		$cuffdiff=checkPath($self->{'software:cuffdiff'});
	} else {
		$cuffdiff=checkPath($self->{'software:cufflinks'});
		$cuffdiff=~s/cufflinks$/cuffdiff/;
	}
	$cuffdiff_cmd .= qq(export cuffdiff=$cuffdiff\n);
	my $para = $self->{'CustomSetting:cuffdiff'};
	$cuffdiff_cmd .= qq(export para="$para"\n);
	my $workdir .= checkPath($self->{"-workdir"});
	$cuffdiff_cmd .= qq(export workdir=$workdir\n);
	$cuffdiff_cmd .= qq(cd \${workdir}\n);
	$cuffdiff_cmd .= "[[ -d cuffdiff] || mkdir cuffdiff\n" if (!-d qq($self->{"-workdir"}/cuffdiff));
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
	my $samtools = checkPath($self->{"software:samtools"});
	my $SamToFastq=checkPath($self->{"software:SamToFastq"});
	my $AlternativeSplicing=checkPath($self->{"software:AlternativeSplicing"});
	my $bwa=checkPath($self->{"software:bwa"});
	my $getAlnGene = checkPath($self->{"software:getAlnGene"});
	my $reference=checkPath($self->{"database:ref"});
	my $refGene=checkPath($self->{"database:refGene"});
	
	my @libraries=sort keys %{$self->{'LIB'}};
	my $runAS_cmd = qq(echo `date`; echo "run AlternativeSplicing"\n);
	$runAS_cmd .= qq(cd $self->{"-workdir"}\n);
	$runAS_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$runAS_cmd .= qq(export REFERENCE=$reference\n);
	$runAS_cmd .= qq(export refGene=$refGene\n);
	$runAS_cmd .= qq(mkdir $self->{"CustomSetting:as_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:as_outdir"}));
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
				$self->{$lib}{"ref-bam"}=$self->{'-workdir'}."/".$self->{"CustomSetting:as_outdir"}."/$lib/$lib.merge.sort.bam";
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
		$self->{$lib}{"AS_result-gff"}=$self->{'-workdir'}."/".$self->{"CustomSetting:as_outdir"}."/$lib/ASresult/AS_result.gff";
		$self->{$lib}{"Junctions-bed"}=$self->{'-workdir'}."/".$self->{"CustomSetting:as_outdir"}."/$lib/ASresult/Junctions.bed";
	}
	$runAS_cmd .= "cd ..\n";
	return ($runAS_cmd);
}

sub runASAP ($) {
	my $self = shift;
	my $samtools = checkPath($self->{"software:samtools"});
	my $asap = checkPath($self->{"software:asap"});
	my $para = (exists $self->{"CustomSetting:asap"}) ? $self->{"CustomSetting:asap"} : "";
	my $reference=checkPath($self->{"database:ref"});
	my $refGene=checkPath($self->{"database:refGene"});
	
	my @libraries=sort keys %{$self->{'LIB'}};
	my $runAS_cmd = qq(echo `date`; echo "run ASAP"\n);
	$runAS_cmd .= qq(cd $self->{"-workdir"}\n);
	$runAS_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$runAS_cmd .= qq(export REFERENCE=$reference\n);
	$runAS_cmd .= qq(export refGene=$refGene\n);
	$runAS_cmd .= qq(mkdir $self->{"CustomSetting:as_outdir"}\n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:as_outdir"}));
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
	my $GenePlot = checkPath($self->{"software:GenePlot"});
	my $plotas_cmd = qq(echo `date`; echo "run GenePlot"\n);
	$plotas_cmd .= qq(cd $self->{"-workdir"}\n);
	$plotas_cmd .= $self->runCuffcompare() if (!exists $self->{'compare-transcript-gtf'});
	#print qq(genelist=$self->{"database:list"}\n);
	my $genelist = checkPath($self->{"database:genelist"});
	my $refGene = checkPath($self->{"database:$gene"});
	
	$plotas_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
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
	my $macs=checkPath($self->{"software:macs"});
	my $para=(exists $self->{'CustomSetting:macs'}) ? $self->{'CustomSetting:macs'} : " --bw=180 --verbose=3 --diag  -B -S --nomodel --mfold 10,30";
	my $macs_cmd = qq(echo `date`; echo "run MACS"\n);
	$macs_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$macs_cmd .= qq(export macs=$macs\n);
	$macs_cmd .= qq(export para="$para"\n);
	$macs_cmd .= qq(cd $self->{"-workdir"}\n);
	my $macs_outdir = (exists $self->{"CustomSetting:macs_outdir"}) ? $self->{"CustomSetting:macs_outdir"} : "macs";
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

sub runRUM ($) {
	
}

sub runCisGenome ($) {
	
}

#########################################################
#                                                       #
#           Gene Prediction/Annotation                  #
#                                                       #
#########################################################

#gth -showintronmaxlen 6000 -gcmaxgapwidth 6000 -o MW_v2c.gth.xml -xmlout -maskpolyatails -paralogs \
#-genomic ../../00.db/MW_v2c.fasta -cdna ../../00.db/cassava_cdna.fa \
#-protein ../../00.db/Euphorbiaceae_protein_sequence.fasta > 1.log 2> 2.log
sub runGenomeThreader ($) {
	my $self=shift;
	if (!exists $self->{"software:gth"} || !-e $self->{"software:gth"})
	{
		return "";
	}
	my $gth=checkPath($self->{"software:gth"});
	my $gth_cmd = qq(echo `date`; echo "run GenomeThreader"\n);
	$gth_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
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
	my $gm=checkPath($self->{"software:gm"});
	my $gm_cmd = qq(echo `date`; echo "run GeneMark"\n);
	$gm_cmd = qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
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

sub runGeneScan ($) {
	
}

sub runPASA ($) {
	my $self=shift;
}

sub runEVM ($) {
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

sub runCRAC($) {
	my $self = shift;
	my $ref = shift;
	$ref ||= "ref";
	if (!exists $self->{"software:crac"} || !defined $self->{"software:crac"}){
		return "";
	}
	my $crac = $self->{"software:crac"};
	my $cracpara = $self->{"CustomSetting:crac"};
	$cracpara .= qq( -k 22 \n) if ($cracpara!~/\-k/);
	my $reference=checkPath($self->{"database:$ref"});
	my $crac_cmd = qq(echo `date`; echo "run CRAC"\n);;
	$crac_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$crac_cmd .= qq(export REFERENCE="$reference"\n);
	$crac_cmd .= qq(export crac=$crac\n);
	$crac_cmd .= qq(export seecerpara="$cracpara"\n);
	$crac_cmd .= qq(export qc_outdir=$self->{"CustomSetting:qc_outdir"}\n);
	$crac_cmd .= qq(cd $self->{"-workdir"}\n);
	$crac_cmd .= qq([[ -d \${qc_outdir} || mkdir \${qc_outdir}n) if (!-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}));
	$crac_cmd .= qq(cd $self->{"CustomSetting:qc_outdir"}\n);
	my @libraries=sort keys %{$self->{'LIB'}};
	foreach my $lib(@libraries) {
		$crac_cmd .= qq([[ -d $lib ]] || mkdir $lib\n) unless (-d qq($self->{"-workdir"}/$self->{"CustomSetting:qc_outdir"}/$lib));
		$crac_cmd .= "cd $lib\n";
		my %fq=getlibSeq($self->{"LIB"}{$lib});
		if (exists $fq{1} && exists $fq{2}){
			for (my $i=0;$i<@{$fq{2}};$i++) {
				my $fq1=${$fq{1}}[$i];
				my $fq2=${$fq{2}}[$i];
				my $bam=(@{$fq{2}}>1)?"$lib.pair.".($i+1).".sam":"$lib.pair.bam";
				my $chimera=(@{$fq{2}}>1)?"$lib.pair.".($i+1).".chimera":"$lib.pair.chimera";
				$crac_cmd .= qq(\${crac} \${cracpara} -i \$REFERENCE -r $fq1 $fq2 -o $bam );
				$crac_cmd .= qq( --paired-end-chimera $chimera ) if ($cracpara=~/chimera/);
				$crac_cmd .= qq(--nb-threads $self->{'CustomSetting:multithreads'}) if (exists $self->{'CustomSetting:multithreads'});
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

sub runCirCOS ($) {
	
}

sub runClustalW2($) {
	my $self=shift;
	my $which=`which clustalw2`;chomp $which;
	my $clustalw2=(-f $which && -e $which) ? $which : checkPath($self->{"software:clustalw2"});
	if (!defined $clustalw2)
	{
		return "";
	}
#clustalw2 -INFILE=chloroplast.ext -ALIGN -TREE -PIM -TYPE=DNA -OUTPUT=NEXUS -STATS=chloroplast_stats.log > chloroplast.log 2>&1 &
	my $clustalw2_cmd=qq(echo `date`; echo "run ClustalW2"\n);;
	$clustalw2_cmd.= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$clustalw2_cmd.=qq(export clustalw2=$clustalw2\n);
	my $para=(exists $self->{"CustomSetting:clustalw2"})?$self->{"CustomSetting:clustalw2"}:"-ALIGN -TREE -PIM -TYPE=DNA -OUTPUT=NEXUS -STATS=stats.log";
	$clustalw2_cmd.=qq(export clustalw2_para=$para\n);
	my $phylogen_outdir = (exists $self->{"CustomSetting:phylogen_outdir"}) ? $self->{"CustomSetting:phylogen_outdir"} : "phylogen";
	my $workdir=checkPath($self->{"-workdir"});
	$clustalw2_cmd.=qq(export workdir=$workdir\n);
	$clustalw2_cmd.=qq(export phylgen_outdir=$phylogen_outdir\n);
	$clustalw2_cmd.=qq(cd \${wordir}\n);
	$clustalw2_cmd.=qq(mkdir \${phylgen_outdir}\n) if (!-d qq($self->{"-workdir"}/$phylogen_outdir));
	$clustalw2_cmd.=qq(cd \${phylgen_outdir}\n);
	if (exists $self->{clustalw2_INFILE})
	{
		my $infile=checkPath($self->{clustalw2_INFILE});
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
	my $mb = (-f $which && -e $which) ? $which : checkPath($self->{"software:mb"});
	if (!defined $mb) {
		return "";
	}
	my $mb_cmd=qq(echo `date`; echo "run MrBayes"\n);
	$mb_cmd .= qq(export PATH=$self->{"CustomSetting:PATH"}:\$PATH\n) if (exists $self->{"CustomSetting:PATH"} && $self->{"CustomSetting:PATH"}!~/\/usr\/local\/bin/);
	$mb_cmd.=qq(export mb=$mb\n);
	my $conf=(exists $self->{"CustomSetting:mb"})?$self->{"CustomSetting:mb"}:"\tlset nst=6 rates=invgamma;\n\tmcmc ngen=20000000;\n\tsump relburnin=yes burninfrac=0.25;\n\tsumt relburnin=yes burninfrac=0.25;\n";
	$mb_cmd.=qq(export mb_conf=$conf\n);
	my $phylogen_outdir = (exists $self->{"CustomSetting:phylogen_outdir"}) ? $self->{"CustomSetting:phylogen_outdir"} : "phylogen";
	my $workdir=checkPath($self->{"-workdir"});
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


sub getlibSeq ($) {
	my $lib=shift;
	my %Seq=();
	for my $i(sort keys %$lib) {
		if (exists $lib->{$i}) {
			my ($A,$B)=($1,$2) if ($i=~/(\w+)([12])/);
			foreach my $reads(@{$lib->{$i}}) {
				if (defined $A && defined $B) {
					my $C=($B==1)?2:1;
					if (exists $lib->{"$A$C"}) {
						push @{$Seq{$B}},$reads;
					} else {
						push @{$Seq{0}},$reads;
					}
				} else {
					push @{$Seq{0}},$reads;
				}
			}
		}
	}
	return %Seq;
}

sub checkIndex($$) {
	my ($soft,$db)=@_;
	if ($soft eq 'bwa') {
		my @refidx=glob("$db.*");
		if (@refidx<5 || !-f "$db.bwt" || !-f $db) {
			print STDERR get_time()."\t\tno bwa index for $db\n";
			return 0;
		} else {
			return 1;
		}
	} elsif ($soft eq 'bowtie') {
		my @refidx=glob("$db.*.ebwt");
		if (@refidx<6 || !-f $db) {
			print STDERR get_time()."\t\tno bowtie index for $db\n";
			return 0;
		} else {
			return 1;
		}
	} elsif ($soft eq 'bowtie2') {
		$db=~s/\.fasta//;
		$db=~s/\.fa//;
		my @refidx=glob("$db.*.bt2");
		if (@refidx<6) {
			print STDERR get_time()."\t\tno bowtie2 index for $db\n";
			return 0;
		} else {
			return 1;
		}
	} elsif ($soft eq "samtools") {
		if ($db=~/fa/i && !-f "$db.fai") {
			print STDERR get_time()."\t\tno samtools index for $db\n";
			return 0;
		} elsif ($db =~ /bam$/ && !-f "$db.bai") {
			print STDERR get_time()."\t\tno bwa index for $db\n";
			return 0;
		} else {
			return 1;
		}
	}
}

sub correctJavaCmd
{
	my ($tools,$heap,$Djavaio)=@_;
	my $cmd="java";
	if (defined $tools && $tools !~ /^java/ && $tools !~ /\-jar/)
	{
		if (defined $heap)
		{
			$heap.="m" if ($heap=~/^\d+$/);
			$cmd.=" -Xmx$heap";
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
			$heap.="m" if ($heap=~/^\d+$/);
			$cmd=~s/java/java\ \-Xmx$heap\ /;
		}
		if (defined $Djavaio && $cmd !~ /\-Djava/)
		{
			$cmd=~s/\-jar/\-Djava\.io\.tmpdir\=$Djavaio\ \-jar/;
		}
	}
	return $cmd;
}

sub check_fileformat ($) {
	my $file=shift;
	if (-B $file) {
		if ($file=~/bam$/i)
		{
			my $index=checkIndex('samtools',$file);
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
		} elsif ($file=~/sff/) {
			return "sff";
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
	foreach my $process(@ary)
	{
		$out .= qq(user=`whoami`\np=\`ps -u \$user -f |grep $process |grep -v grep\`\nwhile [ "\$p" != "" ]\ndo\n\techo "$process is not finish yet! sleep 120s"\n\tsleep 120\n\tp=\`ps -u \$user -f |grep $process |grep -v grep\`\ndone\n);
	}
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

1;
