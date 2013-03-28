# $Id: Analyze.pm,v 1.01 2008/10/05 Swansun Exp $ 

package Statistics::Distributions::Analyze;

use strict;
no warnings;
use Statistics::Descriptive;
use Statistics::Distributions;
use vars qw($VERSION);
use constant SIGNIFICANT => 4; # number of float digits to be returned

$VERSION = '0.01';


1;

# new {{{
sub new {
    my $this = shift;

    warn "[new analyze]\n" if $ENV{DEBUG} >= 2;

    $this = bless {}, $this;

    return $this;
}
# }}}

sub init {
	my $this = shift;

	$this->{array_all} = [];

	## set sub array & total array
	my $i = 0;
	foreach my $array_ref(@{$this->{v}}) {
		push @{$this->{array_all}}, @$array_ref;

		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data($array_ref);

		$this->{'v'.$i} = $stat; 

		$i ++;
	}

	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data($this->{array_all});
	$this->{v_all} = $stat;
	
}


sub Anova {
	my $this = shift;
	my $vector = shift;

	if( ref($vector) eq "ARRAY" ) {
        # good arguments
    } else {
        #croak "argument to variance() too strange";
		die "argument to variance() too strange";
	}

	$this->{v} = $vector;
	$this->init;

	# get value C 
	$this->get_C($this->{v_all});

	# get value SS total
	$this->get_SS($this->{v_all});

	## 取 SS组间
	$this->get_SS_between_group([0..((scalar @$vector) - 1)]); 

	# 取 SS组内
	$this->{value_SS_in_group} = $this->{value_SS} - $this->{value_SS_between_group} ; 

	# get 样本组数 k
	$this->{value_k} = scalar @$vector ;
	
	# v 组间
	$this->{value_v_between_group} = $this->{value_k} - 1;	

	# MS组间
	$this->{value_MS_between_group} = $this->{value_SS_between_group} / $this->{value_v_between_group};  

	# v组内
	$this->get_v_in_group([0..((scalar @$vector) - 1)]); 
	
	## MS组内
	$this->{value_MS_in_group} = $this->{value_SS_in_group} / $this->{value_v_in_group};	

	## 取 F 值
	my $F = $this->{value_MS_between_group} / $this->{value_MS_in_group};

	## 取P值
	my $P =Statistics::Distributions::fprob ($this->{value_v_between_group}, $this->{value_v_in_group}, $F);	# 取P值

	return (&float_precision($F), &float_precision($P));

}

sub T_test{
	my $this = shift;
	my $vector = shift;

	$this->{v} = $vector;
	$this->init;

	my $count_1 = $this->{v0}->count;
	my $mean_1 = $this->{v0}->mean;
	my $std_1 = $this->get_std_deviation($this->{v0});

	my $count_2 = $this->{v1}->count;
	my $mean_2 = $this->{v1}->mean;
	my $std_2 = $this->get_std_deviation($this->{v1});


	my $Sc2 = ( ($std_1 ** 2) * ($count_1 - 1) + ($std_2 ** 2) * ($count_2 - 1) ) /( ($count_1 - 1) +  ($count_2 - 1) );

	my $S_x1_x2 = ( $Sc2 * ( ($count_1+$count_2)/($count_1*$count_2) ) ) ** (1/2);

	my $T = ($mean_1 - $mean_2) / ($S_x1_x2);

	my $v = ($count_1 - 1) +  ($count_2 - 1) ;

	my $P = Statistics::Distributions::tprob ($v, $T); 

	$P = $P * 2;

	return (&float_precision($T), &float_precision($P));
}


sub get_std_deviation {
	my $this = shift;
	my $obj = shift;

	my $variance = $obj->variance();

	my $std_deviation = $variance ** (1/2);

	return $std_deviation;

}

sub get_mean_of_each_array(){
	my $this = shift;

	my $number_of_array = scalar(@{$this->{v}});

	my @averages;

	foreach my $i (0..$number_of_array-1) {
		my $key = 'v' . $i;
		
		my $sum = $this->{$key}->mean();

		push @averages, $this->{$key}->mean();

	}

	$this->{mean_of_each_array} = \@averages;
}


sub get_C(){
	my $this = shift;
	my $obj = shift;

	my $sum = $obj->sum();

	my $count = $obj->count();

	$this->{value_C} = ($sum ** 2) / $count;
	
}

sub get_SS(){
	my $this = shift;
	my $obj = shift;

	my $sum = 0;
	foreach my $num (@{$this->{array_all}}) {
		$sum += $num ** 2;
	}

	my $C = $this->get_C($obj);

	$this->{value_SS} = $sum - $C;
	
}

sub get_SS_between_group(){
	my $this = shift;
	my $number_of_array = shift;

	my $ss_tmp = 0;

	foreach my $i (0..scalar(@$number_of_array)-1) {
		my $key = 'v' . $i;
		
		my $sum = $this->{$key}->sum();
		my $count = $this->{$key}->count();
		$ss_tmp += ($sum ** 2)/$count;

	}

	$this->{value_SS_between_group} = $ss_tmp - $this->{value_C};

}


sub get_v_in_group(){
	my $this = shift;
	my $number_of_array = shift;

	my $v_tmp = 0;

	foreach my $i (0..scalar(@$number_of_array)-1) {
		my $key = 'v' . $i;
		
		my $count = $this->{$key}->count();
		$v_tmp += $count;

	}

	$this->{value_v_in_group} = $v_tmp - $this->{value_k};

}

sub float_precision(){
	my $number = shift;

	my $str = '%.' . SIGNIFICANT . 'f';

	my $num = sprintf($str, $number);

	return $num;
}

__END__

=head1 NAME

    Statistics::Distributions::Analyze - A module to caculate Anova, T-test, etc.

=head1 SYNOPSIS

	use Statistics::Distributions::Analyze;

	my $stats = Statistics::Distributions::Analyze->new();

	my @datas = ([1,2,3], [4,5,6], [7,8,9]);

	# for Anova
	my ($F, $P) = $stats->Anova(\@datas);
	print "F is $F, P is $P\n";

	# for T-test
	my @datas = ([1,2,3], [4,5,6]);
	my ($F, $P) = $stats->T_test(\@datas);
	print "F is $F, P is $P\n";


=head1 OTHERS

More algorithms will be out soon

=head1 AUTHOR

Please contact me with ANY suggestions.

Swansun Huang <swansun95 at gmail.com>

=head1 SEE ALSO

L<Statistics::Descriptive>, L<Statistics::Distributions>

=head1 COPYRIGHT

Copyright (c) 2008 by the Swansun Huang.  All rights reserved.

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.


=cut
