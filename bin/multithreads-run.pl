#!/usr/bin/perl -w
#Author: BENM, binxiaofeng@gmail.com
#Date: 2013-05-22
use strict;
use POSIX;
use Getopt::Long;

my $program_name=$1 if($0=~/([^\/]+)$/);
my $usage=<<USAGE;

Program: use this script to run processes in multiple threads

Usage: perl $program_name  <command_file|command_line>  
	-nt <int>	number of threads to use, default=3
	-verbose	output information of running progress
	-help		output help information to screen
Example:
perl multithreads-run.pl "echo 'hello world';echo '1';sleep 5;ps aux|grep sleep|grep -v grep;echo '2';echo "sleep";sleep 4;ps aux|grep sleep|grep -v grep;echo '3'" -nt 2 --verbose

USAGE

my ($nt,$Verbose,$Help);
my %opts;
GetOptions(\%opts, "nt:i"=>\$nt,"verbose!"=>\$Verbose,"help!"=>\$Help);
die $usage if ( @ARGV==0 ||$Help );
#die qq(perl $0 <processes>\nExample: perl $0 "echo 'hello world';sleep 10;echo '1'" "echo 'hallo kitty';sleep 5; echo '2';"\n) if (@ARGV==0);

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
while(@cmd>0) {
	my @processes=();
	#print STDERR ("\tresidue jobs: ",scalar(@cmd),"\n");
	for (my $i=0;$i<$new_add;$i++) {
		push @processes,shift(@cmd) if (@cmd>0);
	}
	if (@processes>0) {
		$num+=scalar(@processes);
		my $progress=int($num*100/$total+.4999);
		print STDERR "\tthrow out progress: $progress%\n" if $Verbose;
		$new_add=runmulti(\@processes);
		#print STDERR "\tnew add jobs: $new_add\n" if $Verbose;
	}

}
while (keys %pids>0) {
	my $collect;
	while(($collect = waitpid(-1, WNOHANG)) > 0) {
		delete $pids{$collect};
	}
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
		#print STDERR "\tpid: $pid\n";
		$pids{$pid}=1 if (defined $pid && $pid!=0);
		$kid_proc_num++;
		#print STDERR "\tkid processes num: $kid_proc_num\n" if $Verbose;
	}
	while (1) { 
		if($zombies > 0) {
			$zombies = 0;
			my $collect;
			while(($collect = waitpid(-1, WNOHANG)) > 0) {
				delete $pids{$collect};
				$kid_proc_num--;
			} 
			if ($kid_proc_num==0) { return $nt }
			elsif($kid_proc_num<$nt) { return $nt-$kid_proc_num }
			else { next; }
		}
	}
}
