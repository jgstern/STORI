use Statistics::Descriptive;
use Data::Dumper;

my $inputFile1 = $ARGV[0];
my $inputFile2 = $ARGV[1];

my @identityArr =();
my @bitscoreArr =();
my @lengthArr =();

my %scoresHash=();
my %scoresHash2=();
my %scoresHash_p=();
my %scoresHash2_p=();

open input1, $inputFile1;
while (<input1>) {
	my $line = $_;
	chomp($line);
	#               qid  sid    %ident      len misma gap qstar qend ssta send evalue  bit
	if ($line =~ m/(.+)\t(.+)\t(\d+\.\d+)\t(\d+)\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t.+\t(.+)/) {
		my $qid = $1;
		my $sid = $2;
		my $ident = $3;
		
		my $len = $4;
		
		my $bit = $5;
		my $qid_s = substr($qid,0,25);
		my $sid_s = substr($sid,0,25);
		
		#print "$ident $bit $len\n";
		
		$scoresHash{$qid_s}{$sid_s}{"ident"} = $ident;
		$scoresHash{$qid_s}{$sid_s}{"bit"} = $bit;
		$scoresHash{$qid_s}{$sid_s}{"len"} = $len;
		
	}
}
close input1;

open input2, $inputFile2;
while (<input2>) {
	my $line = $_;
	chomp($line);
	#               qid  sid    %ident       len misma gap qstar qend ssta send evalue  bit
	if ($line =~ m/(.+)\t(.+)\t(\d+\.\d+)\t(\d+)\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t.+\t(.+)/) {
		my $qid = $1;
		my $sid = $2;
		my $ident = $3;
		my $len = $4;
		my $bit = $5;
		my $qid_s = substr($qid,0,25);
		my $sid_s = substr($sid,0,25);
		#print "$ident $bit $len\n";
		$scoresHash2{$qid_s}{$sid_s}{"ident"} = $ident;
		$scoresHash2{$qid_s}{$sid_s}{"bit"} = $bit;
		$scoresHash2{$qid_s}{$sid_s}{"len"} = $len;
	}
}
close input2;


foreach my $query (keys %scoresHash) {
	my %temp = %{$scoresHash{$query}};
	foreach my $hit (keys %temp) {
		if (exists $scoresHash2{$hit}{$query}) {
			$scoresHash_p{$query}{$hit} = \%{$scoresHash{$query}{$hit}};
			$scoresHash2_p{$query}{$hit} = \%{$scoresHash2{$hit}{$query}};
		}
	}
}

#print Dumper \%scoresHash2_p;


my $n=0;
foreach my $query (keys %scoresHash_p) {
	my %temp = %{$scoresHash{$query}};
	foreach my $hit (keys %temp) {
		my $identity = $scoresHash{$query}{$hit}{"ident"};
		my $bitscore = $scoresHash{$query}{$hit}{"bit"};
		my $aln_len = $scoresHash{$query}{$hit}{"len"};
		push @identityArr, $identity;
		push @bitscoreArr, $bitscore;
		push @lengthArr, $aln_len;
		$n++;
	}
}

print "n=$n\n";
my $sqrtn = sqrt($n);
my $min = -1;
my $loBar = -1;
my $mean = -1;
my $hiBar = -1;
my $max = -1;

print $inputFile1 . "\n";
print "Attribute -2stderr Mean +2stderr\n";
DoStats(\@identityArr);
print "\%Identity $loBar $mean $hiBar\n";
DoStats(\@bitscoreArr);
print "Bitscore $loBar $mean $hiBar\n";
DoStats(\@lengthArr);
print "len $loBar $mean $hiBar\n";



sub DoStats {
	my $arr_ref = shift(@_);
	$min = Round(Min($arr_ref));
	$max = Round(Max($arr_ref));
	$mean = Mean($arr_ref);
	my $sd = Sd($arr_ref);
	my $stderr = ($sd/$sqrtn);
	$loBar = Round(($mean - (2*$stderr)));
	$hiBar = Round(($mean + (2*$stderr)));
	$mean = Round($mean);
}


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

sub Round {
	my $num=shift(@_);
	if ($num =~ m/\./) {
		$num = sprintf("%.2f", $num);
	}
	return $num;
}

sub Sd { #exclude -1
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
	my $answer = $stat->standard_deviation();
	return $answer;
}
