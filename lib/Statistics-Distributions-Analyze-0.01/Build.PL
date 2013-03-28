use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Statistics::Distributions::Analyze',
    license             => 'perl',
    dist_author         => 'Swansun Huang <swansun95@gmail.com>',
    dist_version_from   => 'lib/Statistics/Distributions/Analyze.pm',
    dist_abstract       => 'A module to caculate Anova, T-test, etc',
    build_requires => {
        'Test::More' => 0,
	'Statistics::Descriptive' => 2.6,
	'Statistics::Distributions' => 1.02,
    },
    add_to_cleanup      => [ 'Statistics-Distributions-Analyze-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
