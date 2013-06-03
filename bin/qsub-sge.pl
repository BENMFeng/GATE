#!/usr/bin/perl

=head1 Name

qsub-sge.pl -- control processes running on linux SGE system

=head1 Description

This program throw the jobs and control them running on linux SGE system. It reads jobs 
from an input shell file. One line is the smallest unit of  a single job, however, you can also specify the 
number of lines to form a single job. For sequential commands, you'd better put them
onto a single line, seperated by semicolon. In anywhere, "&" will be removed 
automatically. The program will terminate when all its jobs are perfectly finished. 

If you have so many jobs, the efficency depends on how many CPUs you can get,
or which queque you have chosen by --queue option. You can use the --maxjob option to 
limit the number of throwing jobs, in order to leave some CPUs for other people. 
When each job consumes long time, you can use the --interval option to increase interval
time for qstat checking , in order to reduce the burden of the head node.

As SGE can only recognize absolute path, so you'd better use absolute path everywhere,
we have developed several ways to deal with path problems:
(1) We have added a function that converting local path to absolute
path automatically. If you like writting absolute path by yourself, then you'd better close this
function by setting "--convert no" option. 
(2) Note that for local path, you'd better write
"./me.txt" instead of only "me.txt", because "/" is the  key mark to distinguish path with
other parameters.  
(3) If an existed file "me.txt" is put in front of the redirect character ">", 
or an un-created file "out.txt" after the redirect character ">", 
the program will add a path "./" to the file automatically. This will avoid much
of the problems which caused by forgetting to write "./" before file name. 
However, I still advise you to write "./me.txt" instead of just "me.txt", this is a good habit.
(4) Please also note that for the re-direct character ">" and "2>", there must be space characters 
both at before and after, this is another good habit.

There are several mechanisms to make sure that all the jobs have been perfectly finished:
(1) We add an auto job completiton mark "This-Work-is-Completed!" to the end of the job, and check it after the job finished
(2) We check "GLIBCXX_3.4.9 not found" to make sure that the C/C++ libary on computing nodes are in good state
(3) We provide a "--secure" option to allow the users define their own job completition mark. You can print a mark
    (for example, "my job complete") to STDERR at the end of your program, and set --secure "my job complete" at 
	this program. You'd better do this when you are not sure about wheter there is bug in your program.
(4) We provide a "--reqsub" option, to throw the unfinished jobs automatically, until all the jobs are 
    really finished. By default, this option is closed, please set it forcely when needed. The maximum 
	reqsub cycle number allowed is 1000.
(5) Add a function to detect the died computing nodes automatically.
(6) Add checking "iprscan: failed" for iprscan
(7) Add a function to detect Eqw status.

Normally, The result of this program contains 3 parts: (Note that the number 24137 is the process Id of this program)
(1) work.sh.24137.globle,     store the shell scripts which has been converted to global path 
(2) work.sh.24137.qsub,       store the middle works, such as job script, job STOUT result, and job STDERR result
(3) work.sh.24137.error,      store the error job list, which has been throwed more than one times.

I advice you to always check the .error file after this program is finished. If it is empty, then everything
is good. Otherwise you should see what is the the problem and try to solve. 

For the resource requirement, by default, the --resource option is set to vf=1.9G, which means the total
memory restriction of one job is 1.9G. By this way, you can throw 8 jobs in one computing node, because the 
total memory restriction of one computing node is 15.5G. If your job exceeds the maximum memory allowed,
then it will be killed forcely. For large jobs, you must specify the --resource option manually, which 
has the same format with "qsub -l" option. If you have many small jobs, and want them to run faster, you
also need to specify a smaller memory requirement, then more jobs will be run at the same time. The key
point is that, you should always consider the memory usage of your program, in order to improve the efficency
of the whole cluster.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Author: Hu Yujie  huyj@genomics.org.cn
  Version: 8.2,  Date: 2008-10-30
  Modified: BENM,  binxiaofeng@gmail.com Date: 2009-6-13
  
=head1 Usage
  
  perl qsub-sge.pl <jobs.txt>
  --queue <str>		specify the queue to use, default all.q
  --interval <num>	set interval time of checking by qstat, default 10 seconds
  --lines <num>		set number of lines to form a job, default 1
  --maxjob <num>	set the maximum number of jobs to throw out, default 30000
  --convert <yes/no>	convert local path to absolute path, default yes  
  --secure <mark>	set the user defined job completition mark, can be any string
  --reqsub		reqsub the unfinished jobs untill they are finished.       
  --resource <str>	set the required resource used in qsub -l option, default vf=1.9G
  --over <str>		kill the jobs which exhausted too much memory, default: (minumium of MEMTOT-0.1G), such as 7.7G
  --prefix <str>	set the user-defined prefix name of the jobs, recommand not over 6 characters
  --verbose		output verbose information to screen   
  --help		output help information to screen  

=head1 Exmple
  
  1.work with default options (the most simplest way)
  perl qsub-sge.pl ./work.sh

  2.work with user specifed options: (to select queue, set checking interval time, set number of lines in each job, and set number of maxmimun running jobs)
  perl qsub-sge.pl --queue all.q -interval 1 -lines 3 -maxjob 10  ./work.sh

  3.do not convert path because it is already absolute path (Note that errors may happen when convert local path to absolute path automatically)
  perl qsub-sge.pl --convert no ./work.sh

  4.add user defined job completion mark (this can make sure that your program has executed to its last sentence)
  perl qsub-sge.pl -inter 1  -secure "my job finish" ./work.sh

  5.reqsub the unfinished jobs until all jobs are really completed (the maximum allowed reqsub cycle is 10000)
  perl qsub-sge.pl --reqsub ./work.sh

  6.work with user defined memory usage
  perl qsub-sge.pl --resource vf=1.9G ./work.sh

  7.recommend combination of usages for common applications (I think this will suit for 99% of all your work)
  perl qsub-sge.pl --queue all.q --resource vf=1.9G --maxjob 10 --reqsub --prefix MyWork ./work.sh

  8.set memory control
  perl qsub-sge.pl --queue all.q --resource vf=1.9G --maxjob 10 --over 15.6G --reqsub --prefix MyWork ./work.sh

=cut


use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

##get options from command line into variables and set default values
my ($Queue, $Interval, $Lines, $Maxjob, $Convert,$Secure,$Reqsub,$Resource,$Over,$Prefix,$Verbose,$Help);
GetOptions(
	"lines:i"=>\$Lines,
	"maxjob:i"=>\$Maxjob,
	"interval:i"=>\$Interval,
	"queue:s"=>\$Queue,
	"convert:s"=>\$Convert,
	"secure:s"=>\$Secure,
	"reqsub"=>\$Reqsub,
	"resource:s"=>\$Resource,
	"over:s"=>\$Over,
	"prefix:s"=>\$Prefix,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
#$Queue ||= "all.q";
$Interval ||= 10;
$Lines ||= 1;
$Maxjob ||= 30000;
$Convert ||= 'yes';
$Resource ||= "vf=1.9G";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $work_shell_file = shift;
my @name=split /\//,$work_shell_file;
my $work_prefix= (defined $Prefix) ? $Prefix : $name[-1];
$work_prefix =~ s/\.sh//g unless ($Prefix);
$Over =uc($Over);
if ($Over =~ /(\S+)G$/i)
{
	$Over = $1;
}
elsif ($Over =~ /(\S+)M$/i)
{
	$Over = $1/1024;
}
else
{
	$Over = `qhost | awk \'\$5~/G/{print \$5}\' |sed \'s/G//\' |awk \'BEGIN{min=100000}{min=(min<\$1)?min:\$1}END{print min}\'`;
	chomp $Over;
	$Over=$Over-0.1;
}
##global variables
my $work_shell_file_globle = $work_shell_file.".$$.globle";
my $work_shell_file_error = $work_shell_file.".$$.error";
my $Work_dir = $work_shell_file.".$$.qsub";
my $current_dir = `pwd`; chomp $current_dir;

if ($Convert =~ /y/i) {
	absolute_path($work_shell_file,$work_shell_file_globle);
}else{
	$work_shell_file_globle = $work_shell_file;
}

## read from input file, make the qsub shell files
my $line_mark = 0;
my $Job_mark="001";
mkdir($Work_dir);
my @Shell;  ## store the file names of qsub sell
open IN, $work_shell_file_globle || die "fail open $work_shell_file_globle";
while(<IN>){
	chomp;
	s/&/;/g;
	next unless($_);
	if ($line_mark % $Lines == 0) {
		open OUT,">$Work_dir/$work_prefix\_$Job_mark.sh" || die "failed creat $work_prefix\_$Job_mark.sh";
		push @Shell,"$work_prefix\_$Job_mark.sh";
		$Job_mark++;
	}
	s/;\s*$//;  ##delete the last character ";", because two ";;" characters will cause error in qsub
	s/;\s*;/;/g;
	print OUT $_."; echo This-Work-is-Completed!\n";

	if ($line_mark % $Lines == $Lines - 1) {
		close OUT;
	}
	
	$line_mark++;
}
close IN;
close OUT;


print STDERR "make the qsub shell files done\n" if($Verbose);


## run jobs by qsub, until all the jobs are really finished
my $qsub_cycle = 1;
while (@Shell) {

	## throw jobs by qsub
	##we think the jobs on died nodes are unfinished jobs
	my %Alljob; ## store all the job IDs of this cycle
	my %Runjob; ## store the real running job IDs of this cycle
	my %Error;  ## store the unfinished jobs of this cycle
	chdir($Work_dir); ##enter into the qsub working directoy
	my $job_cmd = "qsub -cwd -S /bin/sh ";  ## -l h_vmem=16G,s_core=8 
	$job_cmd .= "-q $Queue "  if(defined $Queue); ##set queue
	$job_cmd .= "-l $Resource " if(defined $Resource); ##set resource
	##warn $job_cmd;
	for (my $i=0; $i<@Shell; $i++) {
		while (1) {
			if ($i < $Maxjob || run_count(\%Alljob,\%Runjob) < $Maxjob) {
				my $jod_return = `$job_cmd $Shell[$i]`;
				my $job_id = $1 if($jod_return =~ /Your job (\d+)/);
				$Alljob{$job_id} = $Shell[$i];  ## job id => shell file name
				print STDERR "throw job $job_id in the $qsub_cycle cycle\n" if($Verbose);
				last;
			}else{
				print STDERR "wait for throwing next job in the $qsub_cycle cycle\n" if($Verbose);
				sleep $Interval;
			}
		}
	}
	chdir($current_dir); ##return into original directory 


	###waiting for all jobs fininshed
	while (1) {
		my $run_num = run_count(\%Alljob,\%Runjob);
		last if($run_num == 0);
		print STDERR "There left $run_num jobs runing in the $qsub_cycle cycle\n" if(defined $Verbose);
		sleep $Interval;
	}

	print STDERR "All jobs finished, in the firt cycle in the $qsub_cycle cycle\n" if($Verbose);


	##run the secure mechanism to make sure all the jobs are really completed
	open OUT, ">>$work_shell_file_error" || die "fail create $$work_shell_file_error";
	chdir($Work_dir); ##enter into the qsub working directoy
	foreach my $job_id (sort keys %Alljob) {
		my $shell_file = $Alljob{$job_id};
		
		##read the .o file
		my $content;
		if (-f "$shell_file.o$job_id") {
			open IN,"$shell_file.o$job_id" || warn "fail $shell_file.o$job_id";
			$content = join("",<IN>);
			close IN;
		}
		##check whether the job has been killed during running time
		if ($content !~ /This-Work-is-Completed!/) {
			$Error{$job_id} = $shell_file;
			print OUT "In qsub cycle $qsub_cycle, In $shell_file.o$job_id,  \"This-Work-is-Completed!\" is not found, so this work may be unfinished\n";
		}
		

		##read the .e file
		my $content;
		if (-f "$shell_file.e$job_id") {
			open IN,"$shell_file.e$job_id" || warn "fail $shell_file.e$job_id";
			$content = join("",<IN>);
			close IN;
		}
		##check whether the C/C++ libary is in good state
		if ($content =~ /GLIBCXX_3.4.9/ && $content =~ /not found/) {
			$Error{$job_id} = $shell_file;
			print OUT "In qsub cycle $qsub_cycle, In $shell_file.e$job_id,  GLIBCXX_3.4.9 not found, so this work may be unfinished\n";
		}

		##check whether iprscan is in good state
		if ($content =~ /iprscan: failed/) {
			$Error{$job_id} = $shell_file;
			print OUT "In qsub cycle $qsub_cycle, In $shell_file.e$job_id, iprscan: failed , so this work may be unfinished\n";
		}
		
		##check the user defined job completion mark
		if (defined $Secure && $content !~ /$Secure/) {
			$Error{$job_id} = $shell_file;
			print OUT "In qsub cycle $qsub_cycle, In $shell_file.o$job_id,  \"$Secure\" is not found, so this work may be unfinished\n";
		}


	}

	##make @shell for next cycle, which contains unfinished tasks
	@Shell = ();
	foreach my $job_id (sort keys %Error) {
		my $shell_file = $Error{$job_id};
		push @Shell,$shell_file;
	}
	
	$qsub_cycle++;
	if($qsub_cycle > 10000){
		print OUT "\n\nProgram stopped because the reqsub cycle number has reached 10000, the following jobs unfinished:\n";
		foreach my $job_id (sort keys %Error) {
			my $shell_file = $Error{$job_id};
			print OUT $shell_file."\n";
		}
		print OUT "Please check carefully for what errors happen, and redo the work, good luck!";
		die "\nProgram stopped because the reqsub cycle number has reached 10000\n";
	}

	chdir($current_dir); ##return into original directory 
	close OUT;
	print STDERR "The secure mechanism is performed in the $qsub_cycle cycle\n" if($Verbose);

	last unless(defined $Reqsub);
}

print STDERR "\nqsub-sge.pl finished\n" if($Verbose);


####################################################
################### Sub Routines ###################
####################################################

sub absolute_path{
	my($in_file,$out_file)=@_;
	my($current_path,$shell_absolute_path);

	#get the current path ;
	$current_path=`pwd`;   
	chomp $current_path;

	#get the absolute path of the input shell file;
	if ($in_file=~/([^\/]+)$/) {
		my $shell_local_path=$`;
		if ($in_file=~/^\//) {
			$shell_absolute_path = $shell_local_path;		
		}
		else{$shell_absolute_path="$current_path"."/"."$shell_local_path";}
	}	
	
	#change all the local path of programs in the input shell file;
	open (IN,"$in_file");
	open (OUT,">$out_file");
	while (<IN>) {
		chomp;
		##s/>/> /; ##convert ">out.txt" to "> out.txt"
		##s/2>/2> /; ##convert "2>out.txt" to "2> out.txt"
	my @words=split /\s+/, $_;
		
		##improve the command, add "./" automatically
		for (my $i=1; $i<@words; $i++) {
			if ($words[$i] !~ /\//) {
				if (-f $words[$i]) {
					$words[$i] = "./$words[$i]";
				}elsif($words[$i-1] eq ">" || $words[$i-1] eq "2>"){
					$words[$i] = "./$words[$i]";
				}
			}
			
		}
		for (my $i=0;$i<@words ;$i++) {
			if (($words[$i]!~/^\//) && ($words[$i]=~/\//)) {
				$words[$i]= "$shell_absolute_path"."$words[$i]";
				}
			}
	print OUT join("  ", @words), "\n";
	}
	close IN;
	close OUT;
}


##get the IDs and count the number of running jobs
##the All job list and user id are used to make sure that the job id belongs to this program 
##add a function to detect jobs on the died computing nodes.
sub run_count {
	my $all_p = shift;
	my $run_p = shift;
	my $run_num = 0;

	%$run_p = ();
	my $user = `whoami`; chomp $user;
	my @jobs = split /\n/,`qstat -u $user`;
	foreach my $job_line (@jobs) {
		$job_line =~s/^\s+//;
		my @job_field = split /\s+/,$job_line;
		next if($job_field[3] ne $user);
		if (exists $all_p->{$job_field[0]}){
			my %died;
			died_nodes(\%died); ##the compute node is down
			my $node_name = $1 if($job_field[7] =~ /(compute-\d+-\d+)/);
			#check_memory_use($job_field[0]);
			if (exists $died{$node_name} || ($job_field[4] ne "qw" && $job_field[4] ne "r" && $job_field[4] ne "t") ) { ##当有节点死掉的时候，使用qdel杀掉任务，并且将其放到%Error中
				`qdel $job_field[0]`;
				if ($died{$node_name}==1)
				{
					open OUT, ">>$work_shell_file_error" || die "fail create $$work_shell_file_error";
					print OUT "\n\nProgram of $job_field[2] killed because the node has been dead, the following jobs $job_field[2] running in $node_name killed:\n";
					close OUT;
				}
				elsif ($died{$node_name}==2)
				{
					open OUT, ">>$work_shell_file_error" || die "fail create $$work_shell_file_error";
					print OUT "\n\nProgram of $job_field[2] killed because it has exhausted too much memory, this job in $node_name killed:\n";
					close OUT;
				}
			}else{
				$run_p->{$job_field[0]} = $job_field[2]; ##job id => shell file name
				$run_num++;
			}
		}
	}

	return $run_num; ##qstat结果中的处于正常运行状态的任务，不包含那些在已死掉节点上的僵尸任务
}


##HOSTNAME                ARCH         NCPU  LOAD  MEMTOT  MEMUSE  SWAPTO  SWAPUS
##compute-0-24 lx26-amd64 8 - 15.6G - 996.2M -
sub died_nodes{
	my $died_p = shift;

	my @lines = split /\n/,`qhost`;
	shift @lines; shift @lines; shift @lines;  ##remove the first three title lines

	foreach  (@lines) {
		my @t = split /\s+/;
		my $node_name = $t[0];
		my $memory_tot =  $1 if ($t[4] =~ /(\S+)G$/i);
		if ($t[4] =~ /(\S+)M$/i)
		{
			$memory_tot =  $1/1024;
		}
		my $memory_use = $t[5];
		if($memory_use eq '-')
		{
			$died_p->{$node_name} = 1;
		}
		elsif ($memory_use =~ /(\S+)M$/i)
		{
			$memory_use = $1/1024;
		}
		elsif ($memory_use =~ /(\S+)G$/i)
		{
			$memory_use = $1;
		}
		if ($memory_tot>$Over)
		{
			if (($memory_use>=$Over)&&($memory_use/$memory_tot>=0.95)) #memory control
			{
				$died_p->{$node_name} = 2;
			}
		}
		else
		{
			if (($memory_tot>0)&&($memory_use/$memory_tot>=0.95))
			{
				$died_p->{$node_name} = 2;
			}
		}
	}
}

sub check_memory_use
{
	my $job_id = shift;
	`qstat -j $job_id | awk \'\$0~/resource_list\:/{print \$0}\'`;
	#hard resource_list:         virtual_free=14g,cpu=7
	`qstat -j $job_id | awk \'\$1~/usage/{print \$0}\'`;
	#usage    1:                 cpu=00:01:35, mem=20.45100 GBs, io=0.00000, vmem=426.355M, maxvmem=426.355M
}
