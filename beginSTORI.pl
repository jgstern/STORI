#command: perl /path/to/beginSTORI.pl runNumber /path/to/dir/containing/source/dir taxaFile windowSize finalMaxFams
use lib '/nv/hp10/jstern7/perl5reinstall/lib';
use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';
use Data::Dumper;

my $blastdb_dir = "/nv/hp10/jstern7/scratch/universal120312/blast";
my $blastdbcmdPath = "/nv/hp10/jstern7/STORI3/blastdbcmd";
my $STORIpbsPath = "/nv/hp10/jstern7/STORI3/STORI-pbs_t";
my $home = "/nv/hp10/jstern7/STORI3";
my $runNumber=shift(@ARGV);
my $parentDirPath=shift(@ARGV);
my $searchOrderFile=shift(@ARGV);
my $windowSize=shift(@ARGV);
my $finalMaxFams=shift(@ARGV);


if ((!(defined $runNumber)) || (!(defined $parentDirPath)) || (!(defined $searchOrderFile)) || (!(defined $windowSize)) || (!(defined $finalMaxFams))) {
	print "usage: perl /path/to/beginSTORI.pl runNumber /path/to/dir/containing/source/dir taxaFile windowSize finalMaxFams\n";
	exit;
}

my $sourceDirPathA = $parentDirPath . "/" . $runNumber . "a";
my $sourceDirPathB = $parentDirPath . "/" . $runNumber . "b";


my $queryGiFileA = "query-gis\[" . $runNumber . "a-iter0\].txt";
my $queryGiFileB = "query-gis\[" . $runNumber . "b-iter0\].txt";

my $runA_id=-1;
my $runB_id=-1;

my $runIDfile = $parentDirPath . "/tempRunIDfile.txt";
my $cmd = "touch $runIDfile";
system($cmd);

open beginTemp, ">$runIDfile";
print beginTemp "runA $runA_id\;\nrunB $runB_id\;\n";
close beginTemp;


my $go="yes";
if ((-d $sourceDirPathA) || (-d $sourceDirPathB)) {
	$go="no";
	print "a sourceDirPath already exists. enter yes if you want to continue & overwrite.\n";
	$go=<>;
	
	if ($go =~ m/yes/) {
		my $cmd = "rm -r $sourceDirPathA";
		print "cmd is: " . $cmd . "\n";
		system($cmd);
		my $cmd = "rm -r $sourceDirPathB";
		print "cmd is: " . $cmd . "\n";
		system($cmd);
	}
}
if ($go =~ m/yes/) {
	my $cmd = "mkdir " . $sourceDirPathA;
	print "cmd is: " . $cmd . "\n";
	system($cmd);
	my $cmd = "mkdir " . $sourceDirPathB;
	print "cmd is: " . $cmd . "\n";
	system($cmd);
	
	my $orderAfile = $sourceDirPathA . "/taxa-master[". $runNumber . "a].txt";
	my $orderBfile = $sourceDirPathB . "/taxa-master[". $runNumber . "b].txt";
	
	my $cmd = "cp $home/taxa-master[" . $searchOrderFile . "].txt " .  $orderAfile;
	print "cmd is: " . $cmd . "\n";
	system($cmd);
	my $cmd = "cp $home/taxa-master[" . $searchOrderFile . "].txt " .  $orderBfile;
	print "cmd is: " . $cmd . "\n";
	system($cmd);
	
	my $queryGiPathA = $sourceDirPathA . "/". $queryGiFileA;
	my $queryGiPathB = $sourceDirPathB . "/". $queryGiFileB;
	
	my %annot_hash = ();
	#my @first8 = @{MakeRandomSearchOrder($orderAfile)};
	#my @second8 = @{MakeRandomSearchOrder($orderBfile)};
	my @taxaArr = @{MakeTaxaArr($orderAfile)};
	
	my $satisFlag="no"; my $annotQuery="derp";
	while (!($satisFlag =~ m/yes/)) {
		my $cmd = "rm " . $queryGiPathA;
		print "cmd is: " . $cmd . "\n";
		system($cmd);
		my $cmd = "rm " . $queryGiPathB;
		print "cmd is: " . $cmd . "\n";
		system($cmd);
		$annotQuery="derp";
		print "\nPlease enter an expression to match with protein names: ";
		$annotQuery=<>;
		print "you entered " . $annotQuery . "\n";
		
		print "what offset factor\? (usually 3)\n";
		my $offset_factor = <>;
		
		my $seed1 = MakeQueryGis_new(\@taxaArr, $queryGiPathB, $annotQuery, $offset_factor);
		my $seed2 = MakeQueryGis_new(\@taxaArr, $queryGiPathA, $annotQuery, $offset_factor);
		
		my $r = $seed1 / $seed2;
		
		if ( ($r > 3) || ($r < 0.33) ) {
			print "warning- larger than expected discrepency between number of seed gis. consider repeating the search.\n\n";
		}
		if (($seed1 < ($finalMaxFams/2)) || ($seed2 < ($finalMaxFams/2))) {
			print "warning- one of the seed files is smaller than expected given the number of max families. consider reducing the number of families you expect\n\n";
		}
		
		print "satisfied?\n";
		$satisFlag=<>;
	}
	
	print "\nI have created suggestions for the query-gi files. Please adjust them now, and type \"blastoff\" when you have finished and are ready for the STORI to begin.\nThe files are located here:\n$sourceDirPathA\n$sourceDirPathB\n";
	my $input;
	while (!($input eq "blastoff")) {
		print "3 2 1\>";
		$input=<>;
		chomp($input);
	}
	
	if ($input eq "blastoff") {
		my $command= "msub --timeout=900 -v runNumber=" . $runNumber . "a,sourceDirPath=" . $sourceDirPathA . ",windowSize=" . $windowSize . ",finalMaxFams=" . $finalMaxFams . ",sc=0 " . $STORIpbsPath;
		print "cmd is: " . $command . "\n";
		$runA_id=qx($command);
		sleep 30;
		my $command= "msub --timeout=900 -v runNumber=" . $runNumber . "b,sourceDirPath=" . $sourceDirPathB . ",windowSize=" . $windowSize . ",finalMaxFams=" . $finalMaxFams . ",sc=0 " . $STORIpbsPath;
		print "cmd is: " . $command . "\n";
		$runB_id=qx($command);
		print "liftoff!\n";
	}
}

my $runIDfile = $parentDirPath . "/tempRunIDfile.txt";
my $cmd = "touch $runIDfile";
system($cmd);

open beginTemp, ">$runIDfile";
chomp($runA_id); chomp($runB_id);
$runA_id =~ s/\n*(\d+)/$1/;
$runB_id =~ s/\n*(\d+)/$1/;
print beginTemp "runA $runA_id\;\nrunB $runB_id\;\n";
close beginTemp;

sub MakeRandomSearchOrder {
	my $taxaFile = shift(@_);
	open taxa, $taxaFile;
	my @taxaArr=<taxa>;
	close taxa;
	
	my $count=0;
	foreach my $taxon (@taxaArr)
	{
		chomp($taxon);  #removes any newline characters at the end of the string
		$taxaArr[$count] = $taxon;
		$count++;
	}
	
	fisher_yates_shuffle(\@taxaArr);
	
	open taxa, ">$taxaFile";
	foreach my $taxon (@taxaArr)
	{
		print taxa $taxon . "\n";
	}
	close taxa;
	
	return \@taxaArr;
}

sub MakeTaxaArr {
	my $taxaFile = shift(@_);
	open taxa, $taxaFile;
	my @taxaArr=<taxa>;
	close taxa;
	
	my $count=0;
	foreach my $taxon (@taxaArr)
	{
		chomp($taxon);  #removes any newline characters at the end of the string
		$taxaArr[$count] = $taxon;
		$count++;
	}
	fisher_yates_shuffle(\@taxaArr);
	return \@taxaArr;
}

sub fisher_yates_shuffle {
		my $array = shift;
		my $i;
		for ($i = @$array; --$i; ) {
			my $j = int rand ($i+1);
			next if $i == $j;
			@$array[$i,$j] = @$array[$j,$i];
		}
	}


sub MakeQueryGis_new {
	%annot_hash = ();
	my $arr_ref = shift(@_);
	my $queryGiPath = shift(@_);
	my $queryText = shift(@_);
	my $offset_factor = shift(@_);
	
	chomp($queryText);
	my @smallTaxaArr=@{$arr_ref};
	print "searching taxID ";
	foreach my $taxon (@smallTaxaArr) {			
		print $taxon . " ";
		my $cmd="$blastdbcmdPath -entry all -db $blastdb_dir\/$taxon -outfmt \"\%g \%t\"";
		my $annot = qx($cmd);
		my @annotArr = split /\n/, $annot;
		
		while (@annotArr) {
			my $subject = shift(@annotArr);
			if ($subject =~ m/($queryText)/) {
				$subject =~ m/(\d+)\s(.+)/;
				my $gi = $1;
				my $annotation = $2;
				$annot_hash{$taxon}{$gi} = $annotation;
			}
		}
		
	}
	
	@smallTaxaArr = keys %annot_hash;  #we only want to choose from taxa that contain hits
	my $taxaSize = ($#smallTaxaArr + 1);
	
	#my $n = 1;
	#if (int(0.05 * $taxaSize) > $n) {
	#	$n = int(0.05 * $taxaSize);
	#}
	print "\nsmallTaxaArr\: $#smallTaxaArr $taxaSize\n";
	fisher_yates_shuffle(\@smallTaxaArr);
	#$n--;
	#my @taxaSample = @smallTaxaArr[0..$n];
	my @taxaSample = ();
	
	my @seedGis=();
	while (($#seedGis < $finalMaxFams) && ($#smallTaxaArr > -1)) {
		my $randomTaxon = shift(@smallTaxaArr);
		push @taxaSample, $randomTaxon;
		my %temp = %{$annot_hash{$randomTaxon}};
		foreach my $gi (keys %temp) {
			push @seedGis, $gi;
		}
	}
	print "taxa sampled\: " . join(" ", @taxaSample) . "\n";
	
	#my @seedGis=();
	#foreach my $taxon (@taxaSample) {
#		my %temp = %{$annot_hash{$taxon}};
	#	foreach my $gi (keys %temp) {
	#		push @seedGis, $gi;
	#	}
	#}
	
	open giFile, ">>$queryGiPath";
	
	my $x=0;
	foreach my $gi (@seedGis) {
		my $randScore = GenerateRandomScore($offset_factor);
		print giFile $gi . " " . $randScore . "\n";
		print giFile "h fam" . $x . "\n\n";
		$x++;
	}
	
	close giFile;
	
	my $ans = $#seedGis;
	print "Number of seed GIs\: $ans\n\n";
	return $ans;
}



sub GenerateRandomScore {
	my $ans;
	
	my $offset_factor = shift(@_);
	my $offset_raw = $offset_factor * ($windowSize * $windowSize);
	
	my $min = (($windowSize * $windowSize) + $offset_raw);
	my $max = (1.1 * $min);
	my $range = ($max - $min);
	my $ans = (int rand ($range)) + $min;
	
	return $ans;
}