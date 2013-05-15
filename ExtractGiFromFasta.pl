# ExtractGiFromFasta.pl
# perl ./ExtractGiFromFasta.pl /path/to/input/dir/ /path/to/output/dir/

use lib '/nv/hp10/jstern7/perl5reinstall/lib';
use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';

my $inputDir = $ARGV[0];
print "inputDir is $inputDir\n";
my @files = <$inputDir/*>;
my $outputDir = $ARGV[1];
print "outputDir is $outputDir\n";
my $outputFile;

print "files " . join(" ", @files) . "\n";

foreach my $file (@files) {
	if ($file =~ m/(\d+)\.fasta/) {
		$outputFile = $outputDir . $1 . "_" . $1;
		print "outputting to $outputFile\n";
		open IN, $file;
		my @gis=();
		while (<IN>) {
			my $line = $_;
			chomp($line);
			if ($line =~ m/gi\|(\d+)/) {
				push @gis, $1;
			}
		}
		close IN;
		
		open OUT, ">$outputFile";
		foreach $gi (@gis) {
			print OUT "$gi\t$gi\t100\n";
		}
		close OUT;
	}
}