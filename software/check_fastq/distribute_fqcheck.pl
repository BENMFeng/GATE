#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $usage=<<"USAGE";

	Script   : draw distribution of base and quality along reads from .fqcheck file
	Usage    : perl $0 <read1_fqcheck_file> [read2_fqcheck_file] -o <out_file_prefix>
	Exmple   : perl $0 s_1_1.fqcheck s_1_2.fqcheck -o s_1.fqcheck

USAGE
my ($out);
GetOptions("o:s"=>\$out);
die $usage unless ($out && $ARGV[0]);
my @fqcheck_files = @ARGV;

my $gnuplot = find_gnuplot ();
my $set_plot = set_plot ($gnuplot);

my @last_cycle = (0);
open OUT1, ">$out.base";
open OUT2, ">$out.qual";
for (my $i=0; $i<@fqcheck_files; $i++) {
	open IN, "<$fqcheck_files[$i]" or die "Error: cannot open $fqcheck_files[$i]";
	$last_cycle[$i+1]=$last_cycle[$i];
	while (my $line = <IN>) {
		next unless ($line =~ /^base/);
		my @value = split /\s+/, $line;
		$value[1] += $last_cycle[$i];
		$last_cycle[$i+1] = $value[1] if ($value[1]>$last_cycle[$i+1]);
		print OUT1 join "\t", @value[1 ... 6];
		print OUT1 "\n";
		for (my $i=7; $i<@value; $i++) {
			print OUT2 "$value[1]\t",($i-7),"\t$value[$i]\n" if ($value[$i]>=5);
		}
	}
	close IN;
}
close OUT1;
close OUT2;

my $plot;
$plot = $set_plot;
$plot .= "set output '$out.base.ps'\n";
$plot .= "set xrange [0:$last_cycle[-1]+1]\n";
$plot .= "set yrange [-2:52]\n";
$plot .= "set title 'Base percentage composition along reads'\n";
$plot .= "set xlabel 'Position along reads'\n";
$plot .= "set ylabel 'Percent'\n";
$plot .= "plot '$out.base' u 1:2 w lp t 'A',";
$plot .= " '' u 1:3 w lp t 'C',";
$plot .= " '' u 1:4 w lp t 'G',";
$plot .= " '' u 1:5 w lp t 'T',";
$plot .= " '' u 1:6 w l  t 'N'";
$plot .= ", '' u ($last_cycle[1]+0.5):(50) w impulses lt 3 notitle", if ($last_cycle[2]);
$plot .= "\n";
open TMP, ">$out.base.tmp";
print TMP $plot;
close TMP;
system "$gnuplot $out.base.tmp";
system "convert ps:$out.base.ps $out.base.png";
unlink "$out.base", "$out.base.tmp", "$out.base.ps";

$plot = $set_plot;
$plot .= "set output '$out.qual.ps'\n";
$plot .= "set xrange [0:$last_cycle[-1]+1]\n";
$plot .= "set yrange [-7:42]\n";
$plot .= "set title 'Distribution of qualities'\n";
$plot .= "set xlabel 'Position along reads'\n";
$plot .= "set ylabel 'Quality'\n";
$plot .= "set key off\n";
$plot .= "set style fill solid noborder\n";
$plot .= "plot '$out.qual' u (\$1):(\$2):(\$3*0.0005):(0.3) w boxxy lt 3;\n";
open TMP, ">$out.qual.tmp";
print TMP $plot;
close TMP;
system "$gnuplot $out.qual.tmp";
system "convert ps:$out.qual.ps $out.qual.png";
unlink "$out.qual", "$out.qual.tmp", "$out.qual.ps";

sub find_gnuplot {
	my $gnuplot = "/usr/local/bin/gnuplot";
	$gnuplot = "gnuplot" unless (-e $gnuplot);
	return $gnuplot;
}

sub set_plot {
	my $gnuplot = shift;
	my $version = `$gnuplot -V`;
	$version or die "Error: cannot execute $gnuplot";

	my $set_plot;
	if ($version =~ /gnuplot 4\.0/) {
		$set_plot = "set terminal postscript portrait color\n";
		$set_plot .= "set size 0.914, 0.48\n";
	}
	else {
		$set_plot = "set terminal postscript portrait color size 6.4, 4.8\n";
	}
	$set_plot .= "set grid\n";
	return $set_plot;
}

__END__
