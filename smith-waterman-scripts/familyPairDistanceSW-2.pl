#usage: perl familyPairDistanceSW-2.pl <A.fasta> <B.fasta>

use strict; use warnings;
use Statistics::Descriptive;
my @identity_m=();
my @bitscore_m=();
my @length_m=();

my $STORMdir = "/nv/hp10/jstern7/STORI3";

for (my $x=0; $x<3; $x++) {

	my $Afile = $ARGV[0];
	my $Bfile = $ARGV[1];
	
	my $cmd = $STORMdir . "/smith-waterman-scripts/ssearch36 -p -q -a -w 80 -m 8 -z 11 -f -11 -g -1 -s BL62 -E 100 " . $Afile . " " .  $Bfile . " > $STORMdir/smith-waterman-scripts/temp.ssearch36.1";
	system($cmd);
	$cmd = $STORMdir . "/smith-waterman-scripts/ssearch36 -p -q -a -w 80 -m 8 -z 11 -f -11 -g -1 -s BL62 -E 100 " . $Bfile . " " .  $Afile . " > $STORMdir/smith-waterman-scripts/temp.ssearch36.2";
	system($cmd);
	
	$cmd = "perl $STORMdir/smith-waterman-scripts/summarize2_ssearch36.pl $STORMdir/smith-waterman-scripts/temp.ssearch36.1 $STORMdir/smith-waterman-scripts/temp.ssearch36.2";
	my $output = qx($cmd);
	
	#print $output . "\n";
	
	my @identity=(); 
	my @bitscore=();
	my @length=();
	
	if ($output =~ m/Identity\s(.+)\n/) {
		@identity = split / /, $1;
	}
	if ($output =~ m/Bitscore\s(.+)\n/) {
		@bitscore = split / /, $1;
	}
	if ($output =~ m/len\s(.+)\n/) {
		#print "aln lengths $1\n";
		@length = split / /, $1;
		
	}
	push @identity_m, [@identity];
	push @bitscore_m, [@bitscore];
	push @length_m, [@length];
}

my @id_lo = ($identity_m[0][0], $identity_m[1][0], $identity_m[2][0]);
my @id_me = ($identity_m[0][1], $identity_m[1][1], $identity_m[2][1]);
my @id_hi = ($identity_m[0][2], $identity_m[1][2], $identity_m[2][2]);
my $identity_lo = Min(\@id_lo);
my $identity_mean = Mean(\@id_me);
my $identity_max = Max(\@id_hi);

my @bit_lo = ($bitscore_m[0][0], $bitscore_m[1][0], $bitscore_m[2][0]);
my @bit_me = ($bitscore_m[0][1], $bitscore_m[1][1], $bitscore_m[2][1]);
my @bit_hi = ($bitscore_m[0][2], $bitscore_m[1][2], $bitscore_m[2][2]);
my $bitscore_lo = Min(\@bit_lo);
my $bitscore_mean = Mean(\@bit_me);
my $bitscore_max = Max(\@bit_hi);

my @len_lo = ($length_m[0][0], $length_m[1][0], $length_m[2][0]);
my @len_me = ($length_m[0][1], $length_m[1][1], $length_m[2][1]);
my @len_hi = ($length_m[0][2], $length_m[1][2], $length_m[2][2]);
my $length_lo = Min(\@len_lo);
my $length_mean = Mean(\@len_me);
my $length_max = Max(\@len_hi);

print "Attribute -2stderr Mean +2stderr\n";
print "identity $identity_lo $identity_mean $identity_max\n";
print "bitscore $bitscore_lo $bitscore_mean $bitscore_max\n";
print "length $length_lo $length_mean $length_max\n";

sub Mean { #exclude -1
	my($ref_arr) = @_;
	my @tempArr = @{$ref_arr};
	my @array;
	foreach my $elt (@tempArr) {
		if ($elt != -1) {
			push @array, $elt;
		}
	}
	
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@array);
	my $answer = $stat->mean();
	return $answer;
}

sub Max {
	my($ref_arr) = @_;
	my @array = @{$ref_arr};
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@array);
	my $answer = $stat->max();
	return $answer;
}

sub Min {
	my($ref_arr) = @_;
	my @array = @{$ref_arr};
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@array);
	my $answer = $stat->min();
	return $answer;
}