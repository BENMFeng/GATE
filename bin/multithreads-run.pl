#!/usr/bin/perl -w
use strict;
use POSIX;
die qq(perl $0 <processes>\nExample: perl $0 "echo 'hello world';sleep 10;echo '1'" "echo 'hallo kitty';sleep 5; echo '2';"\n) if (@ARGV==0);
my @cmd = @ARGV;
my $zombies = 0;
my $kid_proc_num = 0;

$SIG{CHLD} = sub { $zombies++ };

for(my $i=0; $i<@cmd; $i++) {
	my $pid = fork();
	if( !defined($pid) ) { exit 1; }
	unless($pid) {
		system "$cmd[$i]";
		exit 0;
	}
	$kid_proc_num++;
} 

while (1) { 
	if($zombies > 0) {
		$zombies = 0;
		my $collect;
		while(($collect = waitpid(-1, WNOHANG)) > 0) {
			$kid_proc_num--;
		} 
		if($kid_proc_num==0) { last; }
		else { next; }
	}
} 
