#!/usr/bin/perl -w

use strict;
use POSIX;
use Getopt::Long;

my $program_name=$1 if($0=~/([^\/]+)$/);
my $usage=<<USAGE;

Program: use this script to run processes in multiple threads
#Author: BENM, binxiaofeng\@gmail.com
#Date: v0.3, 2013-05-22

Usage: perl $program_name  <command_file|command_line>  
	-nt <int>	number of threads to use, default=3
	-q		run job batch queuing mode
	-m <str>	run in job break mark mode
	-verbose	output information of running progress
	-help		output help information to screen
Example:
perl $program_name -nt 2 -q -verbose "echo 'hello world';echo '1';sleep 5;ps aux|grep sleep|grep -v grep;echo '2';echo 'sleep';sleep 4;ps aux|grep sleep|grep -v grep;echo '3'" 
perl $program_name -nt 2 -verbose "echo 'hello world';echo '1';sleep 5;ps aux|grep sleep|grep -v grep;echo '2';echo 'sleep';sleep 4;ps aux|grep sleep|grep -v grep;echo '3'"
perl $program_name work.sh
perl $program_name -nt 4 -m break -verbose "echo 'hello world';echo '1';sleep 5;ps aux|grep sleep|grep -v grep;break;echo '2';echo 'sleep';break;sleep 4;ps aux|grep sleep|grep -v grep;break;echo '3'" 

USAGE

my ($nt,$queue,$mark,$Debug,$Verbose,$Help);
my %opts;
GetOptions(\%opts, "nt:i"=>\$nt,"q!"=>\$queue,"m:s"=>\$mark,"debug!"=>\$Debug,"verbose!"=>\$Verbose,"help!"=>\$Help);
die $usage if ( @ARGV==0 ||$Help );

$nt ||= 3;

my $user=`whoami`;

my %pids=();
my @cmd = ();
for (my $i=0;$i<@ARGV;$i++) {
	if (-f $ARGV[$i]) {
		open(IN,$ARGV[$i]) or die "Can't open $ARGV[$i]\n";
		while (<IN>) {
			s/\s*$//;
			push @cmd,$_;
		}
		close IN;
	} else {
		my $all_cmd=$ARGV[$i];
		$all_cmd=~s/\&\&/\;/g;
		push @cmd,split /\;/,$all_cmd;
	}
}


my $total=scalar(@cmd);
my $num=0;
print STDERR "\n\tcmd num:  $total\n\tmulti-threads num:  $nt\n\n" if $Verbose;
my $new_add=$nt;
my $job_mark=0;
while(@cmd>0) {
	my @processes=();
	print STDERR ("\tresidue jobs num: ",scalar(@cmd),"\n") if $Debug;
	for (my $i=0;$i<$new_add;$i++) {
		if ($job_mark==1) {
			checkpid();
			last if (keys %pids>0);
		}
		if (defined $mark && @cmd>0 && $cmd[0] eq $mark){
			shift @cmd;
			$job_mark=1;
			print STDERR "\tbreaking\n" if $Verbose;
			last;
		} else {
			$job_mark=0;
		}
		push @processes,shift(@cmd) if (@cmd>0);
	}
	if (@processes>0) {
		$num+=scalar(@processes);
		my $progress=int($num*100/$total+.4999);
		print STDERR "\tthrow out jobs progress: $progress%\n" if $Verbose;
		$new_add=runmulti(\@processes);
		print STDERR "\tadd jobs num: $new_add\n" if $Debug;
	}

}
while (keys %pids>0) {
	checkpid();
}
print STDERR "\tAll processes were accomplishment\n" if $Verbose;

sub runmulti {
	my $cmd_p=shift;
	my $zombies = 0;
	my $kid_proc_num = 0;
	$SIG{CHLD} = sub { $zombies++ };
	for(my $i=0; $i<@$cmd_p; $i++) {
		my $pid = fork();
		if( !defined($pid) ) { exit 1; }
		unless($pid) {
			exec $$cmd_p[$i];
			exit 0;
		}
		print STDERR "\tpid: $pid\n" if $Debug;;
		$pids{$pid}=1 if (defined $pid && $pid!=0);
		$kid_proc_num++;
		print STDERR "\tkid processes num: $kid_proc_num\n" if $Debug;
	}
	while (1) { 
		if($zombies > 0) {
			$zombies = 0;
			my $collect;
			checkpid(\$kid_proc_num);
			if ($kid_proc_num==0) { return $nt }
			elsif($kid_proc_num<$nt) { if ($queue) {next;} else {return $nt-$kid_proc_num} }
			else { next; }
		}
	}
}

sub checkpid {
	my $kid_proc=shift;
	my $collect;
	while(($collect = waitpid(-1, WNOHANG)) > 0) {
		print STDERR "\tfinished pid: $collect\n" if $Debug;
		delete $pids{$collect};
		$$kid_proc-- if (defined $kid_proc);
	} 
}
