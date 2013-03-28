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