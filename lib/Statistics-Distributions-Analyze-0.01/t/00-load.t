#!perl -T

use Test::More tests => 1;

BEGIN {
	use_ok( 'Statistics::Distributions::Analyze' );
}

diag( "Testing Statistics::Distributions::Analyze $Statistics::Distributions::Analyze::VERSION, Perl $], $^X" );
