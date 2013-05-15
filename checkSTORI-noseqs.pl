#command: perl /path/to/checkSTORI.pl /path/to/file-A /path/to/file-B
#this is only for checking if a run is complete using the jobfiles

use lib '/nv/hp10/jstern7/perl5reinstall/lib';
use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';
use Data::Dumper;
use List::MoreUtils qw / uniq /;

my $inputA=shift(@ARGV);
my $inputB=shift(@ARGV);


my $finishedCount=0;
my @temp = @{GetFams($inputA)}; #read the dataset into memory
my $runtime_a = $temp[1];

my @temp = @{GetFams($inputB)}; #read the dataset into memory
my $runtime_b = $temp[1];



print "runtime_a: $runtime_a\n";
print "runtime_b: $runtime_b\n";


sub GetFams {
	my $inputFile = shift(@_);
	open (input, $inputFile);
	my @fileArr=<input>;
	my $startLine=-1;
	my $runtime=-1;

	my $finishedFlag=0;
	for (my $r=$#fileArr; $r>=($#fileArr-40); $r--) {
		#print "check file: $fileArr[$r]";
		if ($fileArr[$r] =~ m/Epilogue/) {
			#print "match\n";
			$finishedFlag=1; 
			$finishedCount++;
		}
		if ($fileArr[$r] =~ m/vmem\=\d+kb\,walltime\=(.+)\n/) {
			$runtime=$1;
		}
	}
	if ($finishedFlag==1) {
		print "$inputFile has finished.\n";
	}
	
	my %fam = ();

	my @temp = (\%fam, $runtime);
	return \@temp;
}