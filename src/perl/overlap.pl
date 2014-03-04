#!/usr/bin/perl -w

=head1 Announcement

Copyright (c) 2008-2013 BENM(Binxiao) Feng                            
All Rights Reserved.                                                   
Send all comments to BENM - BinxiaoFeng@gmail.com                     

This program is free software: you can redistribute it and/or modify  
it under the terms of the GNU General Public License as published by  
the Free Software Foundation, either version 3 of the License, or     
(at your option) any later version.                                   
This program is distributed in the hope that it will be useful,       
but WITHOUT ANY WARRANTY; without even the implied warranty of        
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
GNU General Public License for more details.                          
You should have received a copy of the GNU General Public License     
along with this program.  If not, see <http://www.gnu.org/licenses/>. 

overlap.pl  --  overlapping data crunching relations between one or two block sets

=head1 Description

This program is designed to find overlap relations of blocks, the input files
can be one or two tables,  just need to point out the column of ID, Stand
and End;This script is modified from fanw's findOverlap.pl

The algorithm is that:
(1)Order: sort the two block sets based on start positions, seperately;
(2)Traversal: walk along the chromosome, and find out overlap between the twoblock sets.
(3)Output:
    (i)if you had chosen  "--C" option, the output table would list the regions gone
        through combining the overlapping, and use the outstretched one instead of
        all the regions each other staying with overlapping; and output all the regions
        without redundance

    (ii)if you had chosen "--E" option, the output table would report who and who
        overlapped, as well as their regions(start and end postition), their own size
        and the overlapped size, it is sourced from fanw's program;

        The output is TAB delimited with each line consisting of
        Col 1:ID
        COl 2~5:FirstTable: RegionID Start End Length
        Col 6:OverlapNumber
        Col 7:OverlapSize
        Col 8:OverlapRate
        Col 9~:SecondTable: RegionID:Start,End,Length,OverlapSize
        Col ...:if exists more than one overlapped region the format is the same as Col 9

    (iii)if you had chosen "--F" option, the output table would list all the regions
        which has been filtered out overlapping with others in a certain degree, just left
        the relatively uniq regions.

=head1 Version

  Author: BENM, binxiaofeng@gmail.com
  Version: 4.1.8,  Date: 2008-2-25 UpDate: 2013-09-06

=head1 Usage

  --i1 <str>					input the first table file
  --f1 <int-int-int> or <int-int-int-int>	input the first table required column [id-start-end] [id1-id2-start-end]
  --i2 <str>					input the second table file
  --f2 <int-int-int> or <int-int-int-int>	input the second table required column [id-start-end] [id1-id2-start-end]
  --OL <num-str>	set the overlapped limited condition, see below "PS", default: 0
  --final	Limit the final overlap limited conditions
  --mN <int>	set the minimuim overlapping number: [0] (only could be used in existing "--E" option
                in "--E O" it means that the OverlapNumber should over this setting value(>)
                in "--E N" it means that the OverlapNumber should less than this setting value(<=)
  --C		combine the overlapped regions without redundance and list all the regions
  --F <[R|D>	filterout the regions with overlap, R: combine two blocks filterout the overlapped regions and renew the table;
                D: delete the overlapped blocks from the first table, which these blocks overlap with the second table.
  --E <[O|N]>	enumerate the overlapped blocks(O): overlap rate between two blocks bigger than the rate set by "--OL"
                or none(N) Redundance regions: overlap rate between two blocks smaller than the rate set by "--OL"
  --X <[O|N|G]>	extract the overlapped segment(O);
                extract the none redundance segment(N)(filterout redundance);
                extract the internal gap(G) regions without any segment covered.
  --M <str> || <1,str> || <2,str> || <1,str:2,str>	only list the region id matching situation if none or all list all:[all]
  --split <str>   set the interpolation for splitting data columns,
                such as [s+,s,t+,t] and some special symbol, default: t (i.e. @col=split/\t/)
  --expand <int-int:int-int>	expand regions size, "fst_left_expand-fst_right_expand:sec_left_expand-sec_right_expand", default: 0-0:0-0 or 0:0 (for both sides) or 0 (for both tables)
  --shrink <int-int:int-int>	shrink regions size, "fst_left_shrink-fst_right_shrink:sec_left_shrink-sec_right_shrink", default: 0-0:0-0 or 0:0 (for both sides) or 0 (for both tables)
  --verbose	output verbose information to screen
  --help	output help information to screen
  --Example	output Example to screen
  --version	output version information

=head1 Note
--OL    means if there are existing overlapping between two blocks and over this
        rate of overlapping in the "first" or "second" section of the tables
        (just could be used in "--E" option); or in the "small" size or in the
        "big" between two blocks, if the blocks were coincident with this rule,
        it will carry out the function work you need, else it will ignore to
        compare these two blocks but go on; if the limited condition is set by
        "num-bpsize", it will limited the overlapping size bigger than the "num"
        bp scale; if the limited condition is set by "numbp", it won't limited
        the overlap type, but the Total OverlapSize must bigger than or equal
        the "num" bp (>=num bp); Limit-type: "big", "small", "first", "second",
        "bpsize"; "--X O" only can be set as "--OL int-bp:int-covdepth"

=head1 Exmple

    check out "overlap.pl --Example"

=cut

use strict;
use Getopt::Long;
use Data::Dumper;

my ($input1,$inform1,$input2,$inform2,$OverlapLimited,$Final,$minNum,$Combine,$Filterout,
$Enumerate,$Extract,$matchingID,$Interpolation,$Expand,$Shrink,$SortRule,$Verbose,$Help,$Example,$Version);
my %opts;
GetOptions(
	\%opts,
	"i1:s"=>\$input1,
	"f1:s"=>\$inform1,
	"i2:s"=>\$input2,
	"f2:s"=>\$inform2,
	"OL:s"=>\$OverlapLimited,
	"final:s"=>\$Final,
	"mN:i"=>\$minNum,
	"C"=>\$Combine,
	"F:s"=>\$Filterout,
	"E:s"=>\$Enumerate,
	"X:s"=>\$Extract,
	"M:s"=>\$matchingID,
	"split:s"=>\$Interpolation,
	"expand:s"=>\$Expand,
	"shrink:s"=>\$Shrink,
	"sort:s"=>\$SortRule,
	"verbose"=>\$Verbose,
	"help"=>\$Help,
	"Example"=>\$Example,
	"version"=>\$Version
);
if (defined $Example)
{
	my $example = qq{
Example:
    1.if there are existing overlapping and the size of overlapping is over 50%
    more than the small section between two blocks, filter out these two blocks
    out of the table

    \$ perl overlap.pl  --i1 table1.txt --f1 0-2-3 --OL 0.5-small --F R >Overlap_result1.txt

    2.if there are existing overlapping and over 50% more than the small section
    between two regions, combine these two regions, meanhile preserve the unique
    regions,including the overlap not exceeding 0.5 of the small section,
    in order to get rid of the redundance

    \$ perl overlap.pl  --i1 table1.txt --f1 0-2-3 --i2 table2.txt --f2 0-2-3 --OL 0.5-small --C >Overlap_result2.txt

    3.list all the none overlapping or total length of overlapped segments less
    than the first table regions

    \$ perl overlap.pl  --i1 table1.txt --f1 0-1-2-3 --i2 table2.txt --f2 0-1-2-3 --OL 0.5-first --mN 2 --E N >Overlap_result3.txt

    4.list the two table with considerable overlapping(it can be controlled by
    the "--mR & mN"), and just list the accurately matching ID of "Variation"
    (set by "--M" option) in first table

    \$ perl overlap.pl --i1 table1.txt --f1 0-1-2-3 --i2 table2.txt --f2 0-1-2-3  --OL 0.1-small --mN 1 --E O --M 1,"^Variation$" >Overlap_result4.txt

    5. Treat file and split by tab and colon \"[\\t\\:]\"

    \$ perl overlap.pl --E O --i1 table1.txt --f1 0-1-2-3 --split \"[\\t\\:]\"
    
    6. Expand or Shrink region before overlap data crunching
    \$ perl overlap.pl --E O --i1 table1.txt --f1 0-1,2,3-2-3 --expand 10:20
    \$ perl overlap.pl --E O --i1 table1.txt --f1 0-1,2,3-2-3 --i2 table2.txt --f2 0-4,5-4-5 --shrink 25-50:30-40
    };
    die "$example\n";
}
die "  Author: BENM, binxiaofeng\@gmail.com\n  Version: 4.1,4  Date: 2012-10-26\n" if (defined $Version);
die `pod2text $0` if (($Help)||((!defined $Enumerate)&&(!defined $Filterout)&&(!$Combine)&&(!$Extract)));

die("
Combine Overlap\n
Usage:   overlap.pl --C\n
Options:
        --i1 <str>	input the first table file [in.table1]
        --f1 <int-int-int>  || <int-int-int-int>	input the first table format [id-start-end] or [id1-id2-start-end]
        --i2 <str> input the first table file [in.table2]
        --f2 <int-int-int>  || <int-int-int-int>	input the first table format [id-start-end] or [id1-id2-start-end]
        --OL <num-str>	set the overlapped limited condition [OverlapRate-LimitedType], LimitedType default:0-small
        --M <str> || <1,str> || <2,str> || <1,str:2,str>	only list the region id matching situation if none or all list all:[Matching ID]
        --C		combine the overlapped regions without redundance and list all the regions [null]\n
Example:1. perl overlap.pl --C --i1 table1.txt --f1 0-2-3 --OL 0 > Combine_Overlap_result1.txt
        2. perl overlap.pl --C --i1 table1.txt --f1 0-2-3 --i2 table2.txt --f2 0-2-3 --OL 0.1-big > Combine_Overlap_result2.txt\n
") if ( ($Combine) && ( (!defined $input1) || (!defined $inform1) ) );

die("
Filterout Overlap\n
Usage:   overlap.pl --F\n
Options:
        --i1 <str>	input the first table file [in.table1]
        --f1 <int-int-int>  || <int-int-int-int>	input the first table format [id-start-end] or [id1-id2-start-end]
        --i2 <str> input the first table file [in.table2]
        --f2 <int-int-int>  || <int-int-int-int>	input the first table format [id-start-end] or [id1-id2-start-end]
        --OL <num-str>	set the overlapped limited condition [OverlapRate-LimitedType], LimitedType default:0-small
        --M <str> || <1,str> || <2,str> || <1,str:2,str>	only list the region id matching situation if none or all list all:[Matching ID]
        --F <R|D>	filterout the regions with overlap
R: join two blocks filterout the overlapped regions and renew the table;
D: delete the overlapped blocks from the first table, which these blocks overlap with the second table.\n
Example:1. perl overlap.pl  --F R --i1 table1.txt --f1 0-2-3 --OL 0.5-small > Filter_Overlap_result1.txt
        2. perl overlap.pl  --F D --i1 table1.txt --f1 0-2-3 --i2 table2.txt --f2 0-2-3 --OL 100-bpsize > Filter_Overlap_result2.txt\n
") if  ( ($Filterout) && ( ( ( ($Filterout ne "R") && ($Filterout ne "D") ) && ( (!defined $input1) || (!defined $inform1) ) ) || ( ( ($Filterout eq "D") && (!defined $input2) ) ) ) );

die("
Enumerate Overlap\n
Usage:   overlap.pl --E <O|N>\n
Options:
        --i1 <str>	input the first table file [in.table1]
        --f1 <int-int-int>  || <int-int-int-int>	input the first table format [id-start-end] or [id1-id2-start-end]
        --i2 <str> input the first table file [in.table2]
        --f2 <int-int-int>  || <int-int-int-int>	input the first table format [id-start-end] or [id1-id2-start-end]
        --OL <num-str>	set the overlapped limited condition [OverlapRate-LimitedType] LimitedType can be null, default: 0
        --final		Limit the final overlap limited conditions
        --mN <int>	set the minimuim limited overlapping number: [1] for N option it is maxmuim, default: 1
        --M <str> || <1,str> || <2,str> || <1,str:2,str>	only list the region id matching situation if none or all list all:[Matching ID]
        --E <[O|N]>	enumerate the overlapped blocks
        --sort <str:str>	sort overlap list by: [id,start,end,size,overlap:big,small], default: overlap:big;
O: list the overlaping result
N: list the non-overlaping result
Example:1. perl overlap.pl --E O --i1 table1.txt --f1 0-1-2-3 --i2 table2.txt --f2 0-1-2-3  --OL 0.1-first --mN 1  --M Variation > List_Overlap_result.txt
        2. perl overlap.pl --E N --i1 table1.txt --f1 0-1-2-3 --i2 table2.txt --f2 0-1-2-3  --OL 100bp --mN 10 > None_Overlap_result.txt\n
") if ( ($Enumerate) && ( ( ($Enumerate ne "O") && ($Enumerate ne "N") ) || ( (!defined $input1) || (!defined $inform1) ) ) );

die("
eXtract Overlap\n
Usage:   overlap.pl --X <O|N>\n
Options:
        --i1 <str>	input the first table file [in.table1]
        --f1 <int-int-int>  || <int--int-int-int>	input the first table format [id-start-end] or [id1-id2-start-end]
        --i2 <str> input the first table file [in.table2]
        --f2 <int-int-int>  || <int--int-int-int>	input the first table format [id-start-end] or [id1-id2-start-end]
        --M <str> || <1,str> || <2,str> || <1,str:2,str>	only list the region id matching situation if none or all list all:[Matching ID]
        --X <[O|N|G]>	extract the overlapped segment(O);
			extract the none redundance segment(N)(filterout redundance);
			extract the gap(G) internal regions without any segment covered.
Example:1. perl overlap.pl --X O --i1 table1.txt --f1 0-2-3 --OL 100-bp:3-covdepth > Overlap_result1.txt
        2. perl overlap.pl --X N --i1 table1.txt --f1 0-2-3 > None_Overlap_result1.txt
        3. perl overlap.pl --X G --i1 table1.txt --f1 0-2-3 > Internal_gaps.txt
        4. perl overlap.pl --X O --i1 table1.txt --f1 0-1-2-3 --i2 table2.txt --f2 0-1-2-3 --M Variation > Overlap_result2.txt
        5. perl overlap.pl --X N --i1 table1.txt --f1 0-2-3 --i2 table2.txt --f2 0-2-3 > None_Overlap_result2.txt\n
") if ( ($Extract) && ( ($Extract ne "O") && ($Extract ne "N") && ($Extract ne "G") || ( (!defined $input1) || (!defined $inform1) ) ) );

## Global setting
$OverlapLimited = (defined $OverlapLimited) ? $OverlapLimited : 0;
$minNum = (defined $minNum) ? $minNum : 1;
$SortRule ||= "overlap:big";
my @orate=();
my ($minRate,$limitType,$minSize);
if (defined $OverlapLimited)
{
	my $ol=$1 if ($OverlapLimited=~/^([^\:]+)/);
	@orate=split /-/,$ol;
	$minRate = ($orate[0] =~ m/^(\d+)$/) ? $orate[0] : "";
	$minRate = 1/1.000000000000001 if ($minRate == 1);
	if ($orate[0] =~ m/(\d+)bp/)
	{
		$minSize = $1;
		$minRate = 0;
	}
	die "minRate:$minRate is not a available value!" if (($minRate !~ m/^\d+$/)&&($minRate !~ m/^\d+\.\d+$/)&&(($minRate !~ m/^\d+bp$/)));
	$limitType = ((defined $orate[1])&&($orate[1] =~ m/\w+/)) ? $orate[1] : "";
}

my @Shift=();
for (my $i=0;$i<2;$i++)
{
	$Shift[$i][0]=0;
	$Shift[$i][1]=0;
}
if (defined $Expand) {
	if ($Expand=~/^(\d+)\-(\d+)\:(\d+)\-(\d+)$/)
	{
		$Shift[0][0]=$1;
		$Shift[0][1]=$2;
		$Shift[1][0]=$3;
		$Shift[1][1]=$4;
	}
	elsif ($Expand=~/^(\d+)\:(\d+)$/)
	{
		$Shift[0][0]=$1;
		$Shift[0][1]=$1;
		$Shift[1][0]=$2;
		$Shift[1][1]=$2;
	}
	elsif ($Expand=~/^(\d+)$/)
	{
		$Shift[0][0]=$1;
		$Shift[0][1]=$1;
		$Shift[1][0]=$1;
		$Shift[1][1]=$1;
	}
}
if (defined $Shrink)
{
	if ($Shrink=~/^(\d+)\-(\d+)\:(\d+)\-(\d+)$/)
	{
		$Shift[0][0]-=$1;
		$Shift[0][1]-=$2;
		$Shift[1][0]-=$3;
		$Shift[1][1]-=$4;
	}
	elsif ($Shrink=~/^(\d+)\:(\d+)$/)
	{
		$Shift[0][0]-=$1;
		$Shift[0][1]-=$1;
		$Shift[1][0]-=$2;
		$Shift[1][1]-=$2;
	}
	elsif ($Shrink=~/^(\d+)$/)
	{
		$Shift[0][0]-=$1;
		$Shift[0][1]-=$1;
		$Shift[1][0]-=$1;
		$Shift[1][1]-=$1;
	}
}

my $bpsize = $minRate if ( (defined $limitType) && ($limitType =~ m/bp/i) );
$minNum ||= 0;
$matchingID ||= "all";
$Interpolation ||="t";

chomp $matchingID;
my $matchingT="";
my ($T1,$matchingT1,$T2,$matchingT2)=("","","","");
if ($matchingID ne "all")
{
	$matchingID=~s/\s+//g;
	if ($matchingID=~m/\:/g)
	{
		($matchingT1,$matchingT2)=split/\:/,$matchingID;
		($T1,$matchingT1)=split/\,/,$matchingT1;
		($T2,$matchingT2)=split/\,/,$matchingT2;
	}
	elsif (($matchingID=~m/\,/)&&($matchingID !~ m/\:/g))
	{
		($matchingT,$matchingID)=split/\,/,$matchingID;
		if ($matchingT==1)
		{
			$T1=1;
			$matchingT1=$matchingID;
		}
		if ($matchingT==2)
		{
			$T2=2;
			$matchingT2=$matchingID;
		}
	}
	elsif (($matchingID !~m/\,/)&&($matchingID !~ m/\:/g))
	{
		($T1,$matchingT1,$T2,$matchingT2)=(1,$matchingID,2,$matchingID);
	}
}

my ($id1,$rd1,$s1,$e1)=("","","","");
my ($id2,$rd2,$s2,$e2)=("","","","");
if (defined $inform1)
{
	if (table_col($inform1) == 4)
	{
		($id1,$rd1,$s1,$e1)=table_col($inform1);
	}
	else
	{
		($id1,$s1,$e1)=table_col($inform1);
	}
}

if (defined $inform2)
{
	if (table_col($inform2) == 4)
	{
		($id2,$rd2,$s2,$e2)=table_col($inform2);
	}
	else
	{
		($id2,$s2,$e2)=table_col($inform2);
	}
}

my %FirstTable=read_table($input1,$id1,$s1,$e1,$rd1,$T1,$matchingT1,1) if (defined $input1);
my %SecondTable=read_table($input2,$id2,$s2,$e2,$rd2,$T2,$matchingT2,2) if (defined $input2);

if ( (defined $input1) && (!defined $input2) )
{

	combine_overlap(\%FirstTable) if (defined $Combine);
    
	filterout_overlap(\%FirstTable) if ((defined $Filterout)&&($Filterout eq "R"));
    
	enumerate_overlap(\%FirstTable,\%FirstTable,$Enumerate) if (defined $Enumerate);
    
	extract_overlap(\%FirstTable) if (defined $Extract);

}

if (defined $input2)
{
	enumerate_overlap(\%FirstTable,\%SecondTable,$Enumerate) if (defined $Enumerate);
	filterout_overlap(\%FirstTable,\%SecondTable) if ((defined $Filterout)&&($Filterout eq "D"));
    
	if ((defined $Combine) ||  (defined $Filterout) || (defined $Extract))
	{
		my %CombineTable;
		foreach  my $id1 (sort keys %FirstTable)
		{
			my $table1=$FirstTable{$id1};
			push @{$CombineTable{$id1}},@$table1;
			delete $FirstTable{$id1};
		}
		delete @FirstTable {keys %FirstTable};
		foreach  my $id2 (sort keys %SecondTable)
		{
			my $table2=$SecondTable{$id2};
			push @{$CombineTable{$id2}},@$table2;
			delete $SecondTable{$id2};
		}
		delete @SecondTable {keys %SecondTable};
    
		combine_overlap(\%CombineTable) if (defined $Combine);
		filterout_overlap(\%CombineTable) if ((defined $Filterout)&&($Filterout eq "R"));
		extract_overlap(\%CombineTable) if (defined $Extract);
	}
}

print STDERR "\n Horray! Overlap work has completed!\n" if (defined $Verbose);

####################################################
################### Sub Routines ###################
####################################################

sub table_col
{
    my $string=shift;
    chomp $string;
    my @table=split /\-/,$string;
    return @table;
}

sub read_table
{
	my ($file,$id,$s,$e,$rd,$t,$matching_id,$input_num)=@_;
	my %Table=();
    
	open (IN,$file) || die "can't open $file for reading:$!";
	while (<IN>)
	{
		chomp;
		next if ($_ eq "");
		my @col;
		if ($Interpolation eq "t")
		{
			@col=split /\t/,$_;
		}
		elsif ($Interpolation eq "t+")
		{
			@col=split /\t+/,$_;
		}
		elsif ($Interpolation eq "s")
		{
			@col=split /\s/,$_;
		}
		elsif ($Interpolation eq "s+")
		{
			@col=split /\s+/,$_;
		}
		else
		{
			@col=split/$Interpolation/,$_;
		}
		my $ID=assign_id($id,\@col);
		next if ((defined $Enumerate)&&($input_num==2)&&(!exists $FirstTable{$ID}));
		if ( (defined $col[$s]) && (defined $col[$e]) && ($col[$s] =~ /\d+/) && ($col[$e] =~ /\d+/) )
		{
			my $RegionID="";
			if ((defined $rd)&&($rd ne ""))
			{
				$RegionID=assign_id($rd,\@col);
			}
			if (defined $Expand || defined $Shrink)
			{
				$col[$s]=$col[$s]-$Shift[$input_num-1][0];
				$col[$e]=$col[$e]+$Shift[$input_num-1][1];
			}
			my $Start = ($col[$s]>$col[$e]) ? $col[$e]: $col[$s];
			my $End = ($col[$s]>$col[$e]) ? $col[$s]: $col[$e];
			print STDERR "$_ this region is not standar format, and it has been corrected!\n" if (($col[$s]>$col[$e])&&($Verbose));
			if ((defined $rd)&&($rd ne ""))
			{
				push @{$Table{$ID}},[$Start,$End,$RegionID] if ((!defined $t)||($t eq "")||($RegionID =~ m/$matching_id/i)||($matchingID eq "all"));
			}
			else
			{
				push @{$Table{$ID}},[$Start,$End];
			}
		}
	}
	close IN;
    
	return %Table;
}

sub assign_id
{
	my ($id_col,$col)=@_;
	my $ID="";
	if (($id_col =~ /\,/)||($id_col =~ /\.\./))
	{
		my @Split=(split /\,/,$id_col);
		foreach my $sp(@Split)
		{
			next if ((!defined $sp)||($sp eq "")||($sp!~/\d+/));
			if ($sp=~/(\d+)\.\.(\d+)/)
			{
				for (my $i=$1;$i<=$2;$i++)
				{
					$ID.="$$col[$i]|";
				}
			}
			else
			{
				$ID.="$$col[$sp]|";
			}
		}
		$ID=~s/\|$//;
	}
	else
	{
		$ID=$$col[$id_col];
	}
	return $ID;
}

sub combine_overlap
{
	my ($hash_p,$printMark)=@_;
	$limitType ||= "small";

	if ((!defined $printMark)||($printMark ne "null"))
	{
		if ($rd1 ne "")
		{
			print "ID\tStart\tEnd\tLength\tCombinedNum\tCombinedTableID:Start,End,Length...\n";
		}
		else
		{
			print "ID\tStart\tEnd\tLength\n";
		}
	}

	foreach  my $table_id (keys %$hash_p)
	{
		my @array_order = sort {$a->[0]<=>$b->[0]} @{$hash_p->{$table_id}};
		my @combine=();
		push @combine,[@{$array_order[0]}];
		$combine[0][2]="$combine[0][2]:$combine[0][0],$combine[0][1],".($combine[0][1]-$combine[0][0]+1) if (@{$combine[0]}>2);
		print STDERR "combining overlap on $table_id\n" if($Verbose);
		for (my $i=1;$i<@array_order;$i++)
		{
			if ($array_order[$i][0]>$combine[-1][1])
			{
				push @combine,[@{$array_order[$i]}];
				$combine[-1][2]="$combine[-1][2]:$combine[-1][0],$combine[-1][1],".($combine[-1][1]-$combine[-1][0]+1) if (@{$combine[-1]}>2);
			}
			else
			{
				my ($S1,$E1,$regionID1)=@{$combine[-1]};
				my ($S2,$E2,$regionID2)=@{$array_order[$i]};
				if ( $OverlapLimited ne "" )
				{
					my @compare=compare_overlap ([$S1,$E1],[$S2,$E2]);
					if (scalar (@compare) ==2 )
					{
						($combine[-1][0],$combine[-1][1])=@compare;
						if ((@{$combine[-1]}>2)&&(defined $regionID2)&&($regionID2 ne ""))
						{
							$combine[-1][2].="\t$regionID2:$S2,$E2,".($E2-$S2+1);
						}
					}
					else
					{
						push @combine,[@{$array_order[$i]}];
						$combine[-1][2]="$combine[-1][2]:$combine[-1][0],$combine[-1][1],".($combine[-1][1]-$combine[-1][0]+1) if (@{$combine[-1]}>2);
					}
				}
			}
		}
		foreach my $table (@combine)
		{
			print ("$table_id\t$table->[0]\t$table->[1]\t",($table->[1]-$table->[0]+1));
			if (@$table>2)
			{
				my @CombineNumber=split /\t+/,$table->[2];
				print ("\t",scalar(@CombineNumber),"\t",$table->[2]);
			}
			print "\n";
		}
	}
}

sub enumerate_overlap
{
	my ($Fst_table,$Sec_table,$List_form)=@_;
	chomp $List_form;
	my $overlap_num=0;
	if  (($List_form eq "O") || ($List_form eq "N") )
	{
		if ( ($rd1 ne "") || (($rd2 ne "")&&(defined $input2)) )
		{
			print "ID\t";
			print "FirstTableID\t" if ($rd1 ne "");
			print "Start\tEnd\tLength\tOverlapNumber\tOverlapSize\tOverlapRate\t";
			print "SecondTableID:" if (($rd2 ne "")||((!defined $input2)&&($rd1 ne "")));
			print "Start,End,Length,OverlapSize...\n" ;
		}
		else
		{
			print "ID\tStart\tEnd\tLength\tOverlapNumber\tOverlapSize\tOverlapRate\tStart,End,Length,OverlapSize...\n" if ($List_form eq "N");
			print "ID\tStart\tEnd\tLength\tOverlapNumber\tOverlapSize\tOverlapRate\tStart,End,Length,OverlapSize...\n" if ($List_form eq "O");
		}
	}
	foreach  my $table_id (sort keys %$Fst_table)
	{
		my @fst_tab = (exists $Fst_table->{$table_id}) ? (sort {$a->[0] <=> $b->[0]} @{$Fst_table->{$table_id}}) : ();
		my @sec_tab = (!defined $input2) ? @fst_tab  : ((exists $Sec_table->{$table_id}) ? (sort {$a->[0] <=> $b->[0]} @{$Sec_table->{$table_id}}) : ());
	
		print STDERR "find overlap on $table_id\n" if($Verbose);
	
		my $sec_begin=0;
	
		for (my $i=0; $i<@fst_tab; $i++)
		{
			my $fst_size = $fst_tab[$i][1] - $fst_tab[$i][0] + 1;
			my @overlap=();
			my $total_overlap_size=0;
			for (my $j=$sec_begin; $j<@sec_tab; $j++)
			{
				next if ( (!defined $input2) && (!defined $inform2) &&
				(($rd1 ne "") && ($sec_tab[$j][2] eq $fst_tab[$i][2]) &&
				($fst_tab[$j][0]==$sec_tab[$j][0]) && ($fst_tab[$j][1]==$sec_tab[$j][1])) );
				if ($sec_tab[$j][1] < $fst_tab[$i][0])
				{
					next;
				}
				if ($sec_tab[$j][0] > $fst_tab[$i][1])
				{
					last;
				}
				$sec_begin=$j if ((scalar @overlap==0)&&(defined $input2));
				if ($limitType ne "")
				{
					my @compare=compare_overlap($fst_tab[$i],$sec_tab[$j]);
					next if ( ($List_form eq "O") && (@compare==4) );
					last if ( ($List_form eq "N") && (@compare==2) );
				}
		
				my $sec_size = $sec_tab[$j][1] - $sec_tab[$j][0] + 1;
				my $overlap_size = overlap_size($fst_tab[$i],$sec_tab[$j]);
				$total_overlap_size += $overlap_size;
		
				if ((defined $sec_tab[$j][2])&&($sec_tab[$j][2] ne ""))
				{
					my $sort_part=$overlap_size;
					if ($SortRule=~/id/i)
					{
						$sort_part=$sec_tab[$j][2];
					}
					elsif ($SortRule=~/start/i)
					{
						$sort_part=$sec_tab[$j][0];
					}
					elsif ($SortRule=~/end/i)
					{
						$sort_part=$sec_tab[$j][1];
					}
					elsif ($SortRule=~/size/i)
					{
						$sort_part=$sec_size;
					}
					push @overlap,["$sec_tab[$j][2]:$sec_tab[$j][0],$sec_tab[$j][1],$sec_size,$overlap_size",$sort_part];
				}
				else
				{
					my $sort_part=$overlap_size;
					if ($SortRule=~/start/i)
					{
						$sort_part=$sec_tab[$j][0];
					}
					elsif ($SortRule=~/end/i)
					{
						$sort_part=$sec_tab[$j][1];
					}
					elsif ($SortRule=~/size/i)
					{
						$sort_part=$sec_size;
					}
					push @overlap,["$sec_tab[$j][0],$sec_tab[$j][1],$sec_size,$overlap_size",$sort_part];
				}
			}
			if ($SortRule =~ /big/i)
			{
				@overlap=sort{$b->[1] cmp $a->[1]}@overlap;
			}
			else
			{
				@overlap=sort{$a->[1] cmp $b->[1]}@overlap;
			}
			my $overlap_out="";
			$overlap_num=scalar(@overlap);
			if ($overlap_num>0)
			{
				map{$overlap_out.="\t$$_[0]"}@overlap;
			}
			my $overlap_rate=$total_overlap_size/$fst_size;
			if (($rd1 ne "")&&($fst_tab[$i][2] ne ""))
			{
				if ( ($List_form eq "O") && ($overlap_num>=$minNum) )
				{
					if ((defined $Final)||($limitType eq ""))
					{
						next if ( (defined $minRate) && ($overlap_rate<$minRate) && (!defined $bpsize) );
					}
					next if ( (defined $minSize) && ($total_overlap_size<$minSize));
				}
				elsif ( ($List_form eq "N") && ($overlap_num<=$minNum) )
				{
					if ((defined $Final)||($limitType eq ""))
					{
						next unless ( ((!defined $bpsize)&&($total_overlap_size<=$fst_size*$minRate)) || ((!defined $minSize)&&(defined $bpsize)&&($total_overlap_size<=$bpsize)) );
					}
					next if ((defined $minSize) && ($total_overlap_size>$minSize));
				}
				else
				{
					next;
				}
				print ($table_id."\t".$fst_tab[$i][2]."\t".$fst_tab[$i][0]."\t".$fst_tab[$i][1]."\t".$fst_size."\t".$overlap_num."\t".$total_overlap_size."\t".$overlap_rate."$overlap_out\n");
			}
			else
			{
				if ( ($List_form eq "O")  && ($overlap_num>=$minNum) )
				{
					if ((defined $Final)||($limitType eq ""))
					{
						next if ((defined $minRate) && ($overlap_rate<$minRate) && (!defined $bpsize) );
					}
					next if ((defined $minSize) && ($total_overlap_size<$minSize));
				}
				elsif ( ($List_form eq "N") && ($overlap_num<=$minNum) )
				{
					if ((defined $Final)||($limitType eq ""))
					{
						next unless ( ((!defined $bpsize)&&($total_overlap_size<=$fst_size*$minRate)) || ((!defined $minSize)&&(defined $bpsize)&&($total_overlap_size<=$bpsize)) );
					}
					next if ((defined $minSize) && ($total_overlap_size>$minSize));
				}
				else
				{
					next;
				}
				print ($table_id."\t".$fst_tab[$i][0]."\t".$fst_tab[$i][1]."\t".$fst_size."\t".$overlap_num."\t".$total_overlap_size."\t".$overlap_rate."$overlap_out\n");
			}
		}
	}
}

sub filterout_overlap
{
	my ($table1_hash,$table2_hash)=@_;
	$limitType ||= "small";
	if (!defined $Extract)
	{
		if ($rd1 ne "")
		{
			print "ID\tFirstTableID\tStart\tEnd\n";
		}
		else
		{
			print "ID\tStart\tEnd\n";
		}
	}
	if ((!defined $table2_hash)||(keys %$table2_hash==0))
	{
		foreach  my $table_id (sort keys %$table1_hash)
		{
			my @table_order = sort {$a->[0]<=>$b->[0]} @{$table1_hash->{$table_id}};
			print STDERR "filterout overlap on $table_id\n" if ($Verbose);
			my @remain=();
			push @remain,[@{$table_order[0]}];
			my @last=();
			push @last,@remain;
			for (my $i=1;$i<@table_order;$i++)
			{
				if ((@remain>0)&&($table_order[$i][0]>$remain[-1][1]))
				{
					for (my $j=0;$j<@remain;$j++)
					{
						print "$table_id\t";
						print "$remain[$j][2]\t" if ((defined $remain[$j][2])&&($remain[$j][2] ne ""));
						print "$remain[$j][0]\t$remain[$j][1]\n";
					}
					@remain=();
					push @remain,[@{$table_order[$i]}];
				}
				else
				{
					my ($S1,$E1,$regionID1)=(@remain>0)?@{$remain[-1]}:@last;
					my ($S2,$E2,$regionID2)=@{$table_order[$i]};
					if ( $OverlapLimited ne "" )
					{
						my @compare=compare_overlap ([$S1,$E1],[$S2,$E2]);
						if (scalar (@compare) ==2 )
						{
							pop @remain if (@remain>0);
							@last=();
							push @last,@compare;
							$last[2] = $regionID2 if ((defined $regionID2)&&($regionID2 ne ""));
						}
						else
						{
							for (my $j=0;$j<@remain;$j++)
							{
								print "$table_id\t";
								print "$remain[$j][2]\t" if ((@{$remain[$j]}>2)&&($remain[$j][2] ne ""));
								print "$remain[$j][0]\t$remain[$j][1]\n";
							}
							push @remain,[@{$table_order[$i]}];
						}
					}
				}
			}
			if (@remain>0)
			{
				for (my $j=0;$j<@remain;$j++)
				{
					print "$table_id\t";
					print "$remain[$j][2]\t" if ((@{$remain[$j]}>2)&&($remain[$j][2] ne ""));
					print "$remain[$j][0]\t$remain[$j][1]\n";
				}
			}
		}
	}
	else
	{
		foreach  my $table_id (sort keys %$table1_hash)
		{
			print STDERR "filterout overlap on $table_id\n" if ($Verbose);
			my @table1_order = (exists $table1_hash->{$table_id})  ? (sort {$a->[0]<=>$b->[0]} @{$table1_hash->{$table_id}}) : ();
			my @table2_order = (exists $table2_hash->{$table_id})  ? (sort {$a->[0]<=>$b->[0]} @{$table2_hash->{$table_id}}) : ();
			if (@table2_order<1)
			{
				for (my $i=0;$i<@table1_order;$i++)
				{
					print "$table_id\t";
					print "$table1_order[$i][2]\t" unless ((defined $Extract)||(@{$table1_order[$i]}<3));
					print "$table1_order[$i][0]\t$table1_order[$i][1]\n";
				}
			}
			else
			{
				my $sec_begin=0;
				for (my $i=0;$i<@table1_order;$i++)
				{
					my @tmp_table=@{$table1_order[$i]};
					my ($S1,$E1,$regionID1)=@tmp_table;
					my $print=1;
					for (my $j=$sec_begin;$j<@table2_order;$j++)
					{
						my ($S2,$E2,$regionID2)=@{$table2_order[$j]};
						if ($E2<$S1)
						{
							$sec_begin=$j;
							next;
						}
						last if ($S2>$E1);
						my @compare=compare_overlap([$S1,$E1],[$S2,$E2]);
						if (@compare==2)
						{
							if ($S2>$S1) {
								print "$table_id\t";
								print "$regionID1\t" unless (@table1_order<3 || !defined $regionID1);
								print "$S1\t".($S2-1)."\n";
							}
							if ($E2<$E1) {
								$tmp_table[0]=$E2+1;
								($S1,$E1,$regionID1)=@tmp_table;
							} else {
								$print=0;
							}
							if ($S1==$S1 && $E2==$E1) {
								$print=0;
							}
						} else {
							last;
						}
					}
					if ($print==1)
					{
						print "$table_id\t";
						print "$regionID1\t" unless (@table1_order<3 || !defined $regionID1);
						print "$S1\t$E1\n";
					}
				}
			}
		}
	}
}

sub extract_overlap
{
	my $block_p=shift;
	my $limitedSize ||= 1;
	my $limitedDepth ||= 1;
	print STDERR "extracting overlap\n" if ($Verbose);
	if (defined $rd1)
	{
		if ($Extract eq "N")
		{
			print "ID\tStart\tEnd\tExtractedID\n";
		}
		elsif ($Extract eq "O")
		{
			if ( ($rd1 ne "") || (($rd2 ne "")&&(defined $input2)) )
			{
				print "ID\tStart\tEnd\tSize\tCoveredRegionsNum\tCoveredRegions(RegionID,Start,End,Length)...\n";
			}
			else
			{
				print "ID\tStart\tEnd\tSize\tCoveredRegionsNum\tCoveredRegions(Start,End,Length)...\n";
			}
			if ($OverlapLimited=~/([^\:\-\s]+)\-bp\:([^\:\-\s]+)\-covdepth/)
			{
				($limitedSize,$limitedDepth)=($1,$2);
			}
		}
	}
	else
	{
		print "ID\tStart\tEnd\n";
	}
	foreach  my $table_id (keys %$block_p)
	{
		print STDERR "extracting for $table_id segment\n" if (defined $Verbose);
		my @array_order = sort {$a->[0]<=>$b->[0]} @{$block_p->{$table_id}};
		if ($Extract eq "N")
		{
			my @nonOverlap=();
			for (my $i=0;$i<@array_order;$i++)
			{
				next if ((@nonOverlap>0)&&($array_order[$i][1]<$nonOverlap[-1][0]));
				if ((@nonOverlap>=$limitedDepth)&&($array_order[$i][1]<=$nonOverlap[-1][1]))
				{
					my $E=$nonOverlap[-1][1];
					if ($array_order[$i][0]<=$nonOverlap[-1][0])
					{
						pop @nonOverlap;
					}
					else
					{
						$nonOverlap[-1][1]=$array_order[$i][0]-1;
					}
					push @nonOverlap,[$array_order[$i][1]+1,$E];
					$nonOverlap[-1][2]=$array_order[$i][2] if (@{$array_order[$i]}>2);
				}
				elsif ((@nonOverlap>0)&&($array_order[$i][0]<=$nonOverlap[-1][1]))
				{
					my $E=$nonOverlap[-1][1];
					if ($array_order[$i][0]<=$nonOverlap[-1][0])
					{
						pop @nonOverlap;
					}
					else{
						$nonOverlap[-1][1]=$array_order[$i][0]-1;
					}
					push @nonOverlap,[$E+1,$array_order[$i][1]];
					$nonOverlap[-1][2]=$array_order[$i][2] if (@{$array_order[$i]}>2);
				}
				else
				{
					push @nonOverlap,[$array_order[$i][0],$array_order[$i][1]];
					$nonOverlap[-1][2]=$array_order[$i][2] if (@{$array_order[$i]}>2);
				}
			}
			foreach my $non(@nonOverlap)
			{
				print ("$table_id\t",(join "\t",@$non),"\n");
			}
		}
		elsif ($Extract eq "O")
		{
			#my %Overlap=();
			my @region=();
			my $fst_begin=0;
			my $sec_begin=1;
			my ($S,$E)=("","");
			for (my $i=$fst_begin;$i<@array_order;$i++)
			{
				next if ($i<$fst_begin);
				my ($S1,$E1,$regionID1)=@{$array_order[$i]};
				$S=$S1;
				$E=$E1;
				if ($limitedDepth<=1)
				{
					if (defined $regionID1)
					{
						push @region,[$S1,$E1,$regionID1];
					}
					else
					{
						push @region,[$S1,$E1];
					}
				}
				for (my $j=$sec_begin;$j<@array_order;$j++)
				{
					next if ($i==$j || $j<$sec_begin);
					my ($S2,$E2,$regionID2)=@{$array_order[$j]};
					if ($E2<$S)
					{
						$sec_begin=$j;
						next;
					}
					elsif ($S2>$E)
					{
						$fst_begin=$j;
						$sec_begin=$j;
						if (@region>=$limitedDepth)
						{
							call_coverdetph($limitedSize,$limitedDepth,$table_id,$S,$E,\@region);
						}
						($S,$E)=("","");
						@region=();
						last;
					}
					if (($S2>=$S)&&($S2<=$E))
					{
						$E=($E2>$E)?$E2:$E;
						#print STDERR "$S\t$E\t$S1\t$E1\t$S2\t$E2\n";
						if ((defined $regionID1)&&(defined $regionID2))
						{
							push @region,[$S1,$E1,$regionID1] if (@region==0);
							push @region,[$S2,$E2,$regionID2];
							$fst_begin=$j+1;
						}
						else
						{
							push @region,[$S1,$E1] if (@region==0);
							push @region,[$S2,$E2];
							$fst_begin=$j+1;
						}
					}
				}
			}
			if (@region>=$limitedDepth)
			{
				call_coverdetph($limitedSize,$limitedDepth,$table_id,$S,$E,\@region);
			}
		}
		elsif ($Extract eq "G")
		{
			my $pre=$array_order[0][1];
			my $pre_id=$array_order[0][2] if (@{$array_order[0]}>2);
			for (my $i=1;$i<@array_order;$i++)
			{
				if ($array_order[$i][0]>$pre+1)
				{
					print ("$table_id\t",($pre+1),"\t",($array_order[$i][0]-1));
					print "\t$pre_id\-$array_order[$i][2]" if (@{$array_order[$i]}>2);
					print "\n";
				}
				if ($pre<$array_order[$i][1])
				{
					$pre=$array_order[$i][1];
					$pre_id=$array_order[$i][2] if (@{$array_order[$i]}>2);
				}
			}
		}
	}
}

sub call_coverdetph
{
	my ($limitedSize,$limitedDepth,$table_id,$S,$E,$region)=@_;
	my %depth=();
	for (my $k=0;$k<@$region;$k++)
	{
		for(my $d=$$region[$k][0];$d<=$$region[$k][1];$d++)
		{
			$depth{$d}{$k}=1;
		}
	}
	my $pos="";
	my %covregion=();
	for (my $k=$S;$k<=$E;$k++)
	{
		my $covdepth=scalar(keys %{$depth{$k}});
		if ($covdepth<$limitedDepth || ($k==$E && $covdepth>=$limitedDepth))
		{
			if ($pos ne "" && $pos=~/\d+/)
			{
				map{$covregion{$_}=1;}keys %{$depth{$k}};
				if ($k-$pos+1>=$limitedSize)
				{
					my $num=scalar(keys %covregion);
					print ("$table_id\t$pos\t",$k,"\t",$k-$pos+1,"\t$num");
					my @covregion=sort{$a<=>$b}keys %covregion;
					my $out="";
					for (my $i=0;$i<@covregion;$i++)
					{
						my $j=$covregion[$i];
						$out.="\t";
						if (@{$$region[$j]}>2 && defined $$region[$j][2] && $$region[$j][2] ne "")
						{
							$out.="$$region[$j][2]:";
						}
						$out.="$$region[$j][0],$$region[$j][1],".($$region[$j][1]-$$region[$j][0]+1);
					}
					print "$out\n";
				}
			}
			$pos="";
			%covregion=();
			next;
		}
		elsif ($covdepth>=$limitedDepth)
		{
			$pos=$k if ($pos eq "");
			map{$covregion{$_}=1;}keys %{$depth{$k}};
		}
	}
}

sub merge_block
{
	my ($block1_p,$block2_p) = @_;
	my $S = ($block1_p->[0] > $block2_p->[0]) ?  $block1_p->[0] : $block2_p->[0];
	my $E = ($block1_p->[1] < $block2_p->[1]) ?  $block1_p->[1] : $block2_p->[1];
	return ($S,$E);
}

sub overlap_size
{
	my ($block1_p,$block2_p) = @_;

	#my $combine_start = ($block1_p->[0] < $block2_p->[0]) ?  $block1_p->[0] : $block2_p->[0];
	#my $combine_end   = ($block1_p->[1] > $block2_p->[1]) ?  $block1_p->[1] : $block2_p->[1];

	#my ($combine_start, $combine_end) = (sort {$a<=>$b} ($block1_p->[0],$block2_p->[0],$block1_p->[1],$block2_p->[1]))[0,-1];
	#my $overlap_size = ($block1_p->[1]-$block1_p->[0]+1) + ($block2_p->[1]-$block2_p->[0]+1) - ($combine_end-$combine_start+1);
	my ($overlap_start,$overlap_end) = (sort {$a<=>$b} ($block1_p->[0],$block2_p->[0],$block1_p->[1],$block2_p->[1]))[1,2];
	return ($overlap_end - $overlap_start + 1);
}

sub compare_overlap
{
	my $block1_p = shift;
	my $block2_p = shift;

	my $S1=$block1_p->[0];
	my $E1=$block1_p->[1];
	my $S2=$block2_p->[0];
	my $E2=$block2_p->[1];
	my $overlapSize=overlap_size($block1_p,$block2_p);
	my $Size1=($E1-$S1+1);
	my $Size2=($E2-$S2+1);
	my ($big,$small);
	my $S=($S1<$S2) ? $S1 : $S2;
	my $E=($E1>$E2) ? $E1 : $E2;
	if ($Size1>$Size2)
	{
		$big=$Size1;
		$small=$Size2;
	}
	else
	{
		$big=$Size2;
		$small=$Size1;
	}
	if ($limitType eq "small")
	{
		if ($overlapSize<=$small*$minRate)
		{
			return ($S1,$E1,$S2,$E2);
		}
		else
		{
			return ($S,$E);
		}
	}
	elsif ($limitType eq "big")
	{
		if ($overlapSize<=$big*$minRate)
		{
			return ($S1,$E1,$S2,$E2);
		}
		else
		{
			return ($S,$E);
		}
	}
	elsif ($limitType eq "first")
	{
		if ($overlapSize<=$Size1*$minRate)
		{
			return ($S1,$E1,$S2,$E2);
		}
		else
		{
			return ($S,$E);
		}
	}
	elsif($limitType eq "second")
	{
		if ($overlapSize<=$Size2*$minRate)
		{
			return ($S1,$E1,$S2,$E2);
		}
		else
		{
			return ($S,$E);
		}
	}
	elsif($limitType =~ /bp/)
	{
		if ($overlapSize<=$bpsize)
		{
			return ($S1,$E1,$S2,$E2);
		}
		else
		{
			return ($S,$E);
		}
	}
	else
	{
		last;
	}
}

__END__
