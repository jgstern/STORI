#command: perl /path/to/continueSTORI.pl /path/to/file-A /path/to/file-B runNumberA runNumberB /path/to/dir/containing/source/dir taxaFile windowSize finalMaxFams

use lib '/nv/hp10/jstern7/perl5reinstall/lib';
use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';
use Data::Dumper;

my $STORIdir = "/nv/hp10/jstern7/STORI3/";
my $STORIpbsPath = $STORIdir . "STORI_48hr-pbs";
my $inputA=shift(@ARGV);
my $inputB=shift(@ARGV);
my $runNumberA=shift(@ARGV);
my $runNumberB=shift(@ARGV);
my $parentDirPath=shift(@ARGV);
if ($parentDirPath =~ m/\/$/) {
	$parentDirPath =~ s/^(.+)\/$/$1/;
}
my $taxaFile=shift(@ARGV);
my $windowSize=shift(@ARGV);
my $finalMaxFams=shift(@ARGV);
my $sc = shift(@ARGV);


my %famsA = %{GetFams($inputA)};		#read the two datasets into memory
#print Dumper \%famsA;
my %famsB = %{GetFams($inputB)};
#print Dumper \%famsB;

#Find families between the two datasets that have common sequences
my @pairScores=();
foreach my $famFromA (keys %famsA) {
	my %consensusFam = ();
	my %chosenFam = ();
	my $chosenFamName=-1;
	my $consensusSize = -1;
	foreach my $famFromB (keys %famsB) {
		my $agreementScore = GetAgreementScore($famsA{$famFromA},$famsB{$famFromB});
		my @temp = ($famFromA, $famFromB, $agreementScore);
		push @pairScores, [@temp];
	}
}
my @pairScores_sorted = sort { $b->[2] <=> $a->[2] } @pairScores;  #sort descending
#print Dumper \@pairScores_sorted;


my $sourceDirPathA = $parentDirPath . "/" . $runNumberA;
my $sourceDirPathB = $parentDirPath . "/" . $runNumberB;

my $orderAfile = $sourceDirPathA . "/taxa-master[". $runNumberA . "].txt";
my $orderBfile = $sourceDirPathB . "/taxa-master[". $runNumberB . "].txt";


my $runA_id=0;
my $runB_id=0;

my $go="yes";
if ($go =~ m/yes/) {
	my $cmd = "mkdir " . $sourceDirPathA;
	print "cmd is: " . $cmd . "\n";
	system($cmd);
	my $cmd = "mkdir " . $sourceDirPathB;
	print "cmd is: " . $cmd . "\n";
	system($cmd);

	if (-e $orderAfile) {
		my $cmd = "cp $orderAfile $orderAfile" . "o";
		system($cmd);
	}
	if (-e $orderBfile) {
		my $cmd = "cp $orderBfile $orderBfile" . "o";
		system($cmd);
	}
	my $cmd = "cp " . $STORIdir . "taxa-master[" . $taxaFile . "].txt " .  $orderAfile;
	print "cmd is: " . $cmd . "\n";
	system($cmd);
	my $cmd = "cp " . $STORIdir . "taxa-master[" . $taxaFile . "].txt " .  $orderBfile;
	print "cmd is: " . $cmd . "\n";
	system($cmd);

	MakeRandomSearchOrder($orderAfile);
	MakeRandomSearchOrder($orderBfile);
	
	my $queryGiPathA = $parentDirPath . "/" .  $runNumberA . "/query-gis[" . $runNumberA . "-iter0].txt"; 
	my $queryGiPathB = $parentDirPath . "/" .  $runNumberB . "/query-gis[" . $runNumberB . "-iter0].txt"; 
	my $cmd = "touch $queryGiPathA";
	system($cmd);
	my $cmd = "touch $queryGiPathB";
	system($cmd);
	open(giFileA, ">$queryGiPathA");
	foreach my $pair_ref (@pairScores_sorted) {
		if ($pair_ref->[2] > 0) {
			my $famA = $pair_ref->[0];
			my $famB = $pair_ref->[1];
			if (exists $famsA{$famA}) {
				if (exists $famsB{$famB}) {
					my %consensus = %{GetConsensusFam($famsA{$famA},$famsB{$famB})};
					my $consensusName = $famA . $famB;
					delete $famsA{$famA};
					delete $famsB{$famB};
					print $consensusName . "\n";
					#print Dumper \%consensus;
					foreach my $gi (values %consensus) {
						print giFileA $gi . " 200\n";
					}
					print giFileA "h " . $consensusName . "\n\n";
				}
			}
		}
	}
	close giFileA;
	my $cmd = "cp $queryGiPathA $queryGiPathB";
	system($cmd);
	
	sleep 30;
	my $command= "msub --timeout=900 -v runNumber=" . $runNumberA . ",sourceDirPath=" . $sourceDirPathA . ",windowSize=" . $windowSize . ",finalMaxFams=" . $finalMaxFams . ",sc=" . $sc . " " . $STORIpbsPath;
	print "cmd is: " . $command . "\n";
	$runA_id=qx($command);
	
	sleep 30;
	my $command= "msub --timeout=900 -v runNumber=" . $runNumberB . ",sourceDirPath=" . $sourceDirPathB . ",windowSize=" . $windowSize . ",finalMaxFams=" . $finalMaxFams . ",sc=" . $sc . " " . $STORIpbsPath;
	print "cmd is: " . $command . "\n";
	
	$runB_id=qx($command);
	
	my $runIDfile = $parentDirPath . "/tempRunIDfile.txt";
	my $cmd = "touch $runIDfile";
	system($cmd);
	
	open beginTemp, ">$runIDfile";
	chomp($runA_id); chomp($runB_id);
	$runA_id =~ s/\n*(\d+)/$1/;
	$runB_id =~ s/\n*(\d+)/$1/;
	print beginTemp "runA $runA_id\;\nrunB $runB_id\;\n";
	close beginTemp;
}





sub GetAgreementScore {
	my %fam1 = %{shift(@_)};
	my %fam2 = %{shift(@_)};
	my $score=0;
	foreach my $taxon (keys %fam1) {
		if (exists $fam2{$taxon}) {
			$fam1{$taxon} =~ m/\{*(\d+)\}*/;
			my $gi1 = $1;
			$fam2{$taxon} =~ m/\{*(\d+)\}*/;
			my $gi2 = $1;
			if ($gi1 == $gi2) {
				$score++;
			}
		}
	}
	return $score;
}

sub GetConsensusFam {
	my %fam1 = %{shift(@_)};
	my %fam2 = %{shift(@_)};
	my %consensus=();
	foreach my $taxon (keys %fam1) {
		if (exists $fam2{$taxon}) {
			$fam1{$taxon} =~ m/\{*(\d+)\}*/;
			my $gi1 = $1;
			$fam2{$taxon} =~ m/\{*(\d+)\}*/;
			my $gi2 = $1;
			if ($gi1 == $gi2) {
				$consensus{$taxon} = $gi1;
			}
		}
	}
	return \%consensus;
}

sub GetFams {
	my $inputFile = shift(@_);
	open (input, $inputFile);
	my @fileArr=<input>;
	my $startLine=-1;
	for (my $r=$#fileArr; $r>=0; $r--) {
		if ($fileArr[$r] =~ m/final\stable\smain/) {
			$startLine = ($r + 2);
			$r=-1;
		}
	}
	
	my %fam = ();
	
	my @famsArr=();
	if ($startLine>0) {
		for (my $r=$startLine; $r<=$#fileArr; $r++) {
			my $line=$fileArr[$r];
			chomp($line);
			if ($line=~ m/^taxon/) {
				$line =~ m/taxon\s(.+)/;
				my $famsTemp = $1;
				chomp($famsTemp);
				@famsArr = split / /, $famsTemp;
				my $i=0;
				foreach my $fam (@famsArr) {   #prevent the family names from getting too long
					my @chars = split //, $fam;
					my $start; my $end;
					my $new = $fam;
					if ($#chars > 200) {
						my $max = ($#chars - 200);
						$start = int( rand($max) );
						$end = ($start + 200);
						my @subset = @chars[$start..$end];
						$new = join("", @subset);
					}
					$famsArr[$i] = $new;
					$i++;
				}
			}
			elsif ($line =~ m/^\d+/) {
				$line =~ m/(\d+)\s\s(.+)/;
				my $gisTemp = $2;
				chomp($gisTemp);
				my @giArr = split / /, $gisTemp;
				my $taxon = $1;
				foreach my $family (@famsArr) {
					if ($family =~ m/.+/) {
						my $gi = shift(@giArr);
						if ($gi != -1) {
							$fam{$family}{$taxon} = $gi;
							#print "fam{$family}{$taxon}=$gi\n";
						}
					}
				}
			}
			else {
				$r=($#fileArr + 1);
			}
		}
	}
	return \%fam;
}

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