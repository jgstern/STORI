=head1 NAME

STORIcontrol.pl

=head1 SYNOPSIS

	perl ./STORIcontrol.pl
	STORI>start <run-name> <scratch/dir> <taxa file> <windowSize> <finalMaxFams>
	STORI>start hemoglobin_eumetazoa_1x_STORI /nv/hp10/jstern7/scratch/STORM3_runfiles eumetazoa 4 20
	STORI>start hemoglobin_euk_8x_STORI /nv/hp10/jstern7/scratch/STORM3_runfiles eukaryota 4 20
	STORI>stop hemoglobin_euk_8x_STORI
	STORI>exit

=head1 DESCRIPTION

This script controls the rest of STORI, except for
STORIstats, which is independent of STORIcontrol.

STORIcontrol allows the user to start and stop protein family retrievals.
Each retrieval ("run") consists of two PBS jobs, meaning that in its present
version, STORI requires a Torque/Moab server and at least 2 compute CPUs.

Most of the computation occurs in STORI.pl, because this is the script that
executes BLASTP searches repeatedly. Each run comprises two independent
instances of STORI.pl, one in the first PBS job, and the other in the second PBS job.

Here is an
overview of the different scripts, to help explain the purpose of STORIcontrol.pl:

=head2 STORIcontrol.pl accepts user input repeatedly

=over 3

=item * If the user commands "start <arguments>"
	
	Then STORIcontrol.pl saves a record of the new run, and begins the run
	by invoking beginSTORI.pl
	
	beginSTORI.pl guides the user through choosing seed proteins
	and passes this information to STORI-pbs_t, which runs STORI.pl
	on one processor for 24 hrs. (2x)
			
=item * If the user commands "stop <run name>"
	
	Then STORIcontrol.pl deletes the record of this run and gives a qdel
	command. NB funny things may occur if the user stops runs that failed
	to start.
		
=item * If the user does nothing
	
	Then after 10 minutes STORIcontrol will cycle through its record of
	existing runs, and update the status of each run using checkSTORI.pl.
	The run parameters are stored in job_data_STORI.txt, @runArr, and %run_params.
	
	If both jobs in a run have completed
	
		Then STORIcontrol will examine the last three agreement scores for the
		run to determine whether or not the run has converged. (Each job could
		be considered a Markov chain, and when the chains are sufficiently similar,
		then the run has converged.)
		
		If the run converged, then it is labelled "converged" and nothing more happens
		with it.
		
		If the run has not converged, then STORIcontrol will initiate either a 24 hour 
		continuation, or a 48 hour continuation. These are simply new instantiations of
		STORI.pl as above, using the result of the jobs that just finished. The 48-hr
		continuations occur if the agreement score of the last jobs was at least 0.9, and
		there are at least 3 complete job-pairs (ie $STORIcount ≥ 2).
		Otherwise the continuation has a 24-hr wall clock limit.
		
	Otherwise
	
		STORIcontrol simply updates the agreement score between the sequence groupings,
		which is an output of checkSTORI.pl. If checkSTORI is checking a run that began
		recently, it is possible that one or both of the runs have not yet generated any
		output, in which case the agreement score = -1. STORIstats will only be able
		to summarize runs whose agreement score is > 0.

=back

=head1 This Document is Supplementary Material to the following study:

Resolving the Tree of Life with STORI: Selectable Taxon Ortholog Retrieval Iteratively

Journal of Molecular Evolution, 2013

Joshua G. Stern and Eric A. Gaucher

School of Biology, Georgia Institute of Technology

Atlanta, GA 30332 USA

eric.gaucher@biology.gatech.edu

=cut

use strict; use warnings;
use lib '/nv/hp10/jstern7/perl5reinstall/lib';
use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';

use Fcntl ':flock';
use Data::Dumper;
use List::MoreUtils qw / uniq /;
use Statistics::Descriptive;
use Getopt::Long;

my $home = "/nv/hp10/jstern7";
my $arrayFile = "/nv/hp10/jstern7/STORI3/job_data_STORI.txt";
my $checkSTORIpath = "/nv/hp10/jstern7/STORI3/checkSTORI.pl";
my $checkSTORIpath2 = "/nv/hp10/jstern7/STORI3/checkSTORI-noseqs.pl";
my $continueSTORIpath_t = "/nv/hp10/jstern7/STORI3/continueSTORIfast_t.pl";
my $continueSTORIpath_48hr = "/nv/hp10/jstern7/STORI3/continueSTORI_48hr.pl";
my $beginSTORIpath_t = "/nv/hp10/jstern7/STORI3/beginSTORI.pl";
my $beginSTORIpath_l = "/nv/hp10/jstern7/STORI3/beginSTORI_l.pl";
my $GetMissingSeqspath = "/nv/hp10/jstern7/STORI3/GetMissingSeqs.pl";
#my $jobStatsFile = "/nv/hp10/jstern7/STORI3/STORIcontrol_job_statistics.txt";
my $help;
GetOptions(
	   'h'             => \$help,
	   );

if( $help ) {
    exec('perldoc',$0);
    exit(0);
}

my $go="yes";
#if (-e $jobStatsFile) {
#	$go="no";
#	print "jobStatsFile already exists. enter yes if you already backed it up and/or wish to overwrite.\n";
#	$go=<>;
#}

#if ($go =~ m/yes/) {
#	open(jobStats, ">$jobStatsFile"); }

my $cmd = "touch $arrayFile";
system($cmd);

my @runArr = @{LoadBackup($arrayFile)};
print "runArr size: $#runArr \n";
my %runHash = ();
my %converged = ();
my %not_converged = ();
my %paused = ();
my %run_params = %{LoadParams($arrayFile)};


print "\n\t\tWelcome to STORI Control!
	   \n\t\tSTORI ©2012 JG Stern & EA Gaucher
	   \n\t\t\tCommands:
	   \t\t\tstart <run-name> <scratch/dir> <taxa file> <windowSize> <finalMaxFams>
	   \t\t\tsetDefaultScratchDir <path>
	   \t\t\tstop <run-name>
	   \t\t\tGetMissingSeqs <1 or more run_names>
	   \t\t\texit\n";


my $repeat=1;
while ($repeat==1) {
	my $input="";
	$SIG{"ALRM"} = sub {die "TimeOut"};  #keep looping even without input
	eval {alarm 600; print "STORI> "; $input = <>; alarm 0};
	if ($@) {print "\n";}
	else {print "You entered $input"}
	
	ParseAndUpdate();
	
	if ($input =~ m/show\s(.+)/) {
		if ($1 =~ m/groups/) {
			print "showing the groups\n"; }
		else {
			print "showing the runs\n"; 
			print Dumper \%runHash;
			print Dumper \%run_params;
			}
	}
	elsif ($input =~ m/start\s(.+)\s(.+)\s(.+)\s(.+)\s(.+)/) {
		my $runName=$1; my $parentDir = $2; my $taxaFile = $3; my $windowSize=$4; my $finalMaxFams=$5;
			print "starting runName=$runName parentDir=$parentDir taxaFile=$taxaFile windowSize=$windowSize finalMaxFams=$finalMaxFams\n";
			StartRetrieval_t($runName, $parentDir, $taxaFile, $windowSize, $finalMaxFams, 0);
	}
	elsif ($input =~ m/startl\s(.+)\s(.+)\s(.+)\s(.+)\s(.+)/) {
		my $runName=$1; my $parentDir = $2; my $taxaFile = $3; my $windowSize=$4; my $finalMaxFams=$5;
			print "starting runName=$runName parentDir=$parentDir taxaFile=$taxaFile windowSize=$windowSize finalMaxFams=$finalMaxFams\n";
			StartRetrieval_l($runName, $parentDir, $taxaFile, $windowSize, $finalMaxFams, 0);
	}
	elsif ($input =~ m/start/) {
		print "usage: start <run-name> <parent/dir> <taxa/file> <windowSize> <finalMaxFams>\n";
	}
	elsif ($input =~ m/domainJump\s(.+)\s(.+)\s(.+)\s(.+)\s(.+)\s(.+)/) {
		my $fromRun=$1; my $newRun=$2; my $taxaFile=$3; my $scratchDir=$4; my $windowSize=$5; my $finalMaxFams=$6;
			print "starting domainJump fromRun=$fromRun newRun=$newRun taxaFile=$taxaFile scratchDir=$scratchDir windowSize=$windowSize finalMaxFams=$finalMaxFams\n";
			DomainJump($fromRun, $newRun, $taxaFile, $scratchDir, $windowSize, $finalMaxFams, 0);
	}
	
	elsif ($input =~ m/stop\s(.+)/) {
		my $runName=$1;
		print "stopping runName=$runName\n";
		StopRun($runName);
	}
	
	elsif ($input =~ m/pause\s(.+)/) {
		my $runName=$1;
		print "pausing runName=$runName\n";
		PauseRun($runName);
	}
	
	elsif ($input =~ m/GetMissingSeqs\s(.+)/) {
		my $runNames=$1;
		my @runsArr = split / /, $runNames;

		foreach my $run (@runsArr) {
			GetMissingSeqs($run, $run_params{$run}{"parentDir"},$run_params{$run}{"taxaFile"});
		}
	}
	elsif ($input =~ m/exit/) {
		$repeat=0;
	}
}

#close jobStats;
print "bye!\n";


sub StopRun {
	my $runName = shift(@_);
	chomp($runName);
	my @STORIruns = (keys %{$runHash{$runName}});
	my $current=Max(\@STORIruns);
	my $idA=$runHash{$runName}{$current}{"ids"}->[0];
	my $idB=$runHash{$runName}{$current}{"ids"}->[1];
	#my $cmd = "qdel $idA";
	my $cmd = "canceljob $idA";
	print "cmd is: $cmd\n";
	system($cmd);
	#$cmd = "qdel $idB";
	$cmd = "canceljob $idB";
	print "cmd is: $cmd\n";
	system($cmd);
	delete $runHash{$runName};
	#print Dumper \%runHash;
	delete $run_params{$runName};
	SaveBackup($arrayFile);
	@runArr = @{LoadBackup($arrayFile)};
	
	#print "(StopRun runHash)\n";
	#print Dumper \%runHash;
}

sub PauseRun {
	my $runName = shift(@_);
	chomp($runName);
	my @STORIruns = (keys %{$runHash{$runName}});
	my $current=Max(\@STORIruns);
	my $idA=$runHash{$runName}{$current}{"ids"}->[0];
	my $idB=$runHash{$runName}{$current}{"ids"}->[1];
	#my $cmd = "qdel $idA";
	my $cmd = "canceljob $idA";
	print "cmd is: $cmd\n";
	system($cmd);
	#$cmd = "qdel $idB";
	$cmd = "canceljob $idB";
	print "cmd is: $cmd\n";
	system($cmd);
	
	#delete $runHash{$runName};
	
	#print Dumper \%runHash;
	$paused{$runName} = 1;
	SaveBackup($arrayFile);
	@runArr = @{LoadBackup($arrayFile)};
	
	#print "(PauseRun runHash)\n";
	#print Dumper \%runHash;
}


sub ContinueRetrieval_48hr {
	my $runName = shift(@_);
	my $parentDir = shift(@_);
	my $taxaFile = shift(@_);
	my $windowSize = shift(@_);
	my $finalMaxFams = shift(@_);
	my $STORIcount = shift(@_);
	
	my @temp=(keys %{$runHash{$runName}});
	my $current = Max(\@temp);

	my $runNameA = $runName . "a";
	my $runNameB = $runName . "b";	
	
	my $sourceFilesDirA = $run_params{$runName}{"parentDir"} . "/" . $runNameA;
	my $sourceFilesDirB = $run_params{$runName}{"parentDir"} . "/" . $runNameB;
	
	my $runA_id=$runHash{$runName}{$current}{"ids"}->[0];
	my $runB_id=$runHash{$runName}{$current}{"ids"}->[1];
	
	my $fileA = $sourceFilesDirA . "/STORI_out_sc" . ($STORIcount-1) . "_" . $runNameA . ".txt";
	my $fileB = $sourceFilesDirB . "/STORI_out_sc" . ($STORIcount-1) . "_" . $runNameB . ".txt";
	
	my $cmd = "perl $continueSTORIpath_48hr $fileA $fileB $runNameA $runNameB $parentDir $taxaFile $windowSize $finalMaxFams $STORIcount";
	print "cmd is: $cmd\n";
	system($cmd);
	my $runIDfile = $parentDir . "/tempRunIDfile.txt";
	#print "attempting open $runIDfile yo wassup\n\n\n";
	print "\n";
	open (beginTemp2, "$runIDfile");
	#sleep 5;
	while (my $output = <beginTemp2>) {
		if ($output =~ m/runA\s(\d+)\;/) {
			$runA_id = $1; }
		if ($output =~ m/runB\s(\d+)\;/) {
			$runB_id = $1; }
		}
	close beginTemp2;
	
	
	my $temp = "ids"; my @tempArr = ($runA_id, $runB_id);
	$runHash{$runName}{$STORIcount}{$temp} = [@tempArr];
	print "(end of ContinueRetrieval) runHash{$runName}{$STORIcount}{ids} = $runA_id $runB_id \n";
	@tempArr = (-1);
	$runHash{$runName}{$STORIcount}{"scores"} = [@tempArr];
	$runHash{$runName}{$STORIcount}{"runtime"} = "-1";
	$run_params{$runName}{"taxaFile"}=$taxaFile;
	$run_params{$runName}{"parentDir"}=$parentDir;
	$run_params{$runName}{"windowSize"}=$windowSize;
	$run_params{$runName}{"finalMaxFams"}=$finalMaxFams;
	
	#print "end of ContinueRetrieval runHash \n";
	#print Dumper \%runHash;
	
	SaveBackup($arrayFile);
	@runArr = @{LoadBackup($arrayFile)};
	sleep 5;
}

sub ContinueRetrieval_t {
	my $runName = shift(@_);
	my $parentDir = shift(@_);
	my $taxaFile = shift(@_);
	my $windowSize = shift(@_);
	my $finalMaxFams = shift(@_);
	my $STORIcount = shift(@_);
	
	my @temp=(keys %{$runHash{$runName}});
	my $current = Max(\@temp);

	my $runNameA = $runName . "a";
	my $runNameB = $runName . "b";	
	
	my $sourceFilesDirA = $run_params{$runName}{"parentDir"} . "/" . $runNameA;
	my $sourceFilesDirB = $run_params{$runName}{"parentDir"} . "/" . $runNameB;
	
	my $runA_id=$runHash{$runName}{$current}{"ids"}->[0];
	my $runB_id=$runHash{$runName}{$current}{"ids"}->[1];
	
	my $fileA = $sourceFilesDirA . "/STORI_out_sc" . ($STORIcount-1) . "_" . $runNameA . ".txt";
	my $fileB = $sourceFilesDirB . "/STORI_out_sc" . ($STORIcount-1) . "_" . $runNameB . ".txt";
	
	my $cmd = "perl $continueSTORIpath_t $fileA $fileB $runNameA $runNameB $parentDir $taxaFile $windowSize $finalMaxFams $STORIcount";
	print "cmd is: $cmd\n";
	system($cmd);
	my $runIDfile = $parentDir . "/tempRunIDfile.txt";
	#print "attempting open $runIDfile yo wassup\n\n\n";
	print "\n";
	open (beginTemp2, "$runIDfile");
	#sleep 5;
	while (my $output = <beginTemp2>) {
		if ($output =~ m/runA\s(\d+)\;/) {
			$runA_id = $1; }
		if ($output =~ m/runB\s(\d+)\;/) {
			$runB_id = $1; }
		}
	close beginTemp2;
	
	
	my $temp = "ids"; my @tempArr = ($runA_id, $runB_id);
	$runHash{$runName}{$STORIcount}{$temp} = [@tempArr];
	print "(end of ContinueRetrieval) runHash{$runName}{$STORIcount}{ids} = $runA_id $runB_id \n";
	@tempArr = (-1);
	$runHash{$runName}{$STORIcount}{"scores"} = [@tempArr];
	$runHash{$runName}{$STORIcount}{"runtime"} = "-1";
	$run_params{$runName}{"taxaFile"}=$taxaFile;
	$run_params{$runName}{"parentDir"}=$parentDir;
	$run_params{$runName}{"windowSize"}=$windowSize;
	$run_params{$runName}{"finalMaxFams"}=$finalMaxFams;
	
	#print "end of ContinueRetrieval runHash \n";
	#print Dumper \%runHash;
	
	SaveBackup($arrayFile);
	@runArr = @{LoadBackup($arrayFile)};
	sleep 5;
}





#GetMissingSeqs($run, $run_params{$run}{"parentDir"}, $run_params{$run}{"taxaFile"});
sub GetMissingSeqs {
	my $runName = shift(@_);
	my $parentDir = shift(@_);
	my $taxaFile = shift(@_);
	
	my @temp=(keys %{$runHash{$runName}});
	my $current = Max(\@temp);

	my $scores_ref=$runHash{$runName}{$current}{"scores"};
	
	my $runNameA = $runName . "a";
	my $runNameB = $runName . "b";	
	
	my $sourceFilesDirA = $run_params{$runName}{"parentDir"} . "/" . $runNameA;
	my $sourceFilesDirB = $run_params{$runName}{"parentDir"} . "/" . $runNameB;
	
	my $fileA = $sourceFilesDirA . "/STORI_out_sc" . $current . "_" . $runNameA . ".txt";
	my $fileB = $sourceFilesDirB . "/STORI_out_sc" . $current . "_" . $runNameB . ".txt";
	
	$runNameA = $runName . "_fin_a";
	$runNameB = $runName . "_fin_b";	
	
	my $cmd = "perl $GetMissingSeqspath $fileA $fileB $runNameA $runNameB $parentDir $taxaFile";
	print "cmd is: $cmd\n";
	system($cmd);
	
	$runName = $runName . "_fin_";
	
	my @tempArr = (-2, -2);
	$runHash{$runName}{"0"}{"ids"} = [@tempArr];
	@tempArr = @{$scores_ref};	#the scores for the finalized run are just the same scores as the input to the finalizer. Not really accurate, but close enough.
	$runHash{$runName}{"0"}{"scores"} = [@tempArr];
	$runHash{$runName}{"0"}{"runtime"} = "1";
	$run_params{$runName}{"taxaFile"}=$taxaFile;
	$run_params{$runName}{"parentDir"}=$parentDir;
	$run_params{$runName}{"windowSize"}=-1;
	$run_params{$runName}{"finalMaxFams"}=-1;
		
	$converged{$runName} = 1;
	
	SaveBackup($arrayFile);
	@runArr = @{LoadBackup($arrayFile)};
	sleep 5;
}



sub StartRetrieval_t {
	#print "(StartRetrieval runHash)\n";
	#print Dumper \%runHash;

	my $runName = shift(@_);
	my $parentDir = shift(@_);
	my $taxaFile = shift(@_);
	my $windowSize = shift(@_);
	my $finalMaxFams = shift(@_);
	my $STORIcount = shift(@_);
	
	#my $output = qx($cmd);
	my $runA_id=-1;
	my $runB_id=-1;
	
	my $cmd = "perl $beginSTORIpath_t $runName $parentDir $taxaFile $windowSize $finalMaxFams";
	print "cmd is: $cmd\n";
	system($cmd);
	my $runIDfile = $parentDir . "/tempRunIDfile.txt";
	#print "attempting open $runIDfile yo wassup\n\n\n";
	print "\n";
	open (beginTemp2, "$runIDfile");
	#sleep 5;
	while (my $output = <beginTemp2>) {
		if ($output =~ m/runA\s(\d+)\;/) {
			$runA_id = $1; }
		if ($output =~ m/runB\s(\d+)\;/) {
			$runB_id = $1; }
		}
	close beginTemp2;
	
	if ( ($runA_id > 0) && ($runB_id > 0) ) {
		my $temp = "ids"; my @tempArr = ($runA_id, $runB_id);
		$runHash{$runName}{$STORIcount}{$temp} = [@tempArr];
		@tempArr = (-1);
		$runHash{$runName}{$STORIcount}{"scores"} = [@tempArr];
		$runHash{$runName}{$STORIcount}{"runtime"} = "-1";
		$run_params{$runName}{"taxaFile"}=$taxaFile;
		$run_params{$runName}{"parentDir"}=$parentDir;
		$run_params{$runName}{"windowSize"}=$windowSize;
		$run_params{$runName}{"finalMaxFams"}=$finalMaxFams;
		SaveBackup($arrayFile);
		@runArr = @{LoadBackup($arrayFile)};
	}
}

sub StartRetrieval_l {
	#print "(StartRetrieval runHash)\n";
	#print Dumper \%runHash;

	my $runName = shift(@_);
	my $parentDir = shift(@_);
	my $taxaFile = shift(@_);
	my $windowSize = shift(@_);
	my $finalMaxFams = shift(@_);
	my $STORIcount = shift(@_);
	
	#my $output = qx($cmd);
	my $runA_id=-1;
	my $runB_id=-1;
	
	my $cmd = "perl $beginSTORIpath_l $runName $parentDir $taxaFile $windowSize $finalMaxFams";
	print "cmd is: $cmd\n";
	system($cmd);
	my $runIDfile = $parentDir . "/tempRunIDfile.txt";
	#print "attempting open $runIDfile yo wassup\n\n\n";
	print "\n";
	open (beginTemp2, "$runIDfile");
	#sleep 5;
	while (my $output = <beginTemp2>) {
		if ($output =~ m/runA\s(\d+)\;/) {
			$runA_id = $1; }
		if ($output =~ m/runB\s(\d+)\;/) {
			$runB_id = $1; }
		}
	close beginTemp2;
	
	if ( ($runA_id > 0) && ($runB_id > 0) ) {
		my $temp = "ids"; my @tempArr = ($runA_id, $runB_id);
		$runHash{$runName}{$STORIcount}{$temp} = [@tempArr];
		@tempArr = (-1);
		$runHash{$runName}{$STORIcount}{"scores"} = [@tempArr];
		$runHash{$runName}{$STORIcount}{"runtime"} = "-1";
		$run_params{$runName}{"taxaFile"}=$taxaFile;
		$run_params{$runName}{"parentDir"}=$parentDir;
		$run_params{$runName}{"windowSize"}=$windowSize;
		$run_params{$runName}{"finalMaxFams"}=$finalMaxFams;
		SaveBackup($arrayFile);
		@runArr = @{LoadBackup($arrayFile)};
	}
}




sub ParseAndUpdate {
	#print "(begin ParseAndUpdate runHash)\n";
	#print Dumper \%runHash;
	
	%paused=();
	%converged=();
	%not_converged=();
	
	if ($#runArr >= 0) {
	foreach my $row_ref (@runArr) {			#load the runArr, which is from the file, into a hash
		my @tempRow = @{$row_ref};
		my $runName = shift(@tempRow);
		my $runid_A = shift(@tempRow);
		my $runid_B = shift(@tempRow);
		my $score=""; 
		my $STORIcount=0;
		
		#print "tempRow= " . join(" ", @tempRow) . "\n";
		
		while (@tempRow) {
			my @scores=();
			while (!($score =~ m/r/)) {
				if ($score =~ m/\d+/) {
					push @scores, $score; }
				$score = shift(@tempRow);
			}
			my $temp="scores";
			$runHash{$runName}{$STORIcount}{$temp} = [@scores];
			
			chomp($score);
			$score =~ m/r(.+)/;
			my $time=$1;
			$temp="runtime";
			$runHash{$runName}{$STORIcount}{$temp} = $time;
			$score = shift(@tempRow);
			
			#populate %converged and %not_converged and %paused, so that we don't confuse ourselves trying to continue runs that have already converged (would only a problem if the converged jobfile has been moved to a different dir)
			if (defined $score) {
				if ($score =~ m/converged/) {
					$converged{$runName} = 1;
					if (exists $not_converged{$runName}) {
						delete $not_converged{$runName};
					}
				}
				elsif ($score =~ m/paused/) {
					$paused{$runName} = 1;
					if (exists $not_converged{$runName}) {
						delete $not_converged{$runName};
					}
				}
			}
			else {
				$not_converged{$runName} = 1;
			}
			
			$STORIcount++;		
		}	
		$STORIcount--;
		my $temp = "ids"; my @tempArr = ($runid_A, $runid_B);
		$runHash{$runName}{$STORIcount}{$temp} = [@tempArr];
	}
	}
	
	
	#print Dumper \%runHash;
	my @tempArr = (keys %not_converged);
	if ($#tempArr >= 0) {
	foreach my $runName (keys %not_converged) {		#cycle through the not-converged runs in runHash, updating scores and runtimes, and start a new run if prev are finished and not converged
		my @STORIcounts = (keys %{$runHash{$runName}});
		my @STORIcounts_sorted = sort { $a <=> $b } @STORIcounts;  #sort ascending
		
		my $STORIcount = $STORIcounts_sorted[$#STORIcounts_sorted];
		my $id_a = $runHash{$runName}{$STORIcount}{"ids"}->[0];
		my $id_b = $runHash{$runName}{$STORIcount}{"ids"}->[1];
		
		my $runNameA = $runName . "a";
		my $runNameB = $runName . "b";
		
		my $sourceFilesDirA = $run_params{$runName}{"parentDir"} . "/" . $runNameA;
		my $sourceFilesDirB = $run_params{$runName}{"parentDir"} . "/" . $runNameB;
		
		my $fileA = $sourceFilesDirA . "/STORI_out_sc" . $STORIcount . "_" . $runNameA . ".txt";
		my $fileB = $sourceFilesDirB . "/STORI_out_sc" . $STORIcount . "_" . $runNameB . ".txt";
		
		
		#first check the jobfile (because if the run reaches its wall clock limit, there will be no indication in the out file)
		my $hr_frac;
		my $sec_frac;
		
		my $runtime_a_f = -1;
		my $runtime_b_f = -1;
		
		my $cmd = "perl $checkSTORIpath2 $home/STORI-pbs.o" . $id_a . " $home/STORI-pbs.o" . $id_b;
		print "cmd is: $cmd\n";
		my $output = qx($cmd);
		
		$output =~ m/runtime_a:\s(.+)\n/;
		my $runtime_a = $1;	
		if (!($runtime_a eq -1)) {
			$runtime_a =~ m/(\d+):(\d+):(\d+)/;
			$hr_frac = ($2 / 60);
			$sec_frac = ($3 / 3600);
			$runtime_a_f = ($1 + $hr_frac + $sec_frac);
		}
		
		$output =~ m/runtime_b:\s(.+)\n/;
		my $runtime_b = $1;
		if (!($runtime_b eq -1)) {
			$runtime_b =~ m/(\d+):(\d+):(\d+)/;
			$hr_frac = ($2 / 60);
			$sec_frac = ($3 / 3600);
			$runtime_b_f = ($1 + $hr_frac + $sec_frac);
		}
		
		#now check the out file
		$cmd = "perl $checkSTORIpath $fileA $fileB -q";
		print "cmd is: $cmd\n";
		$output = qx($cmd);
		
		$output =~ m/runtime_a:\s(.+)\n/;
		$runtime_a = $1;	
		if (!($runtime_a eq -1)) {
			$runtime_a =~ m/(\d+):(\d+):(\d+)/;
			$hr_frac = ($2 / 60);
			$sec_frac = ($3 / 3600);
			$runtime_a_f = ($1 + $hr_frac + $sec_frac);
		}
		
		$output =~ m/runtime_b:\s(.+)\n/;
		$runtime_b = $1;
		if (!($runtime_b eq -1)) {
			$runtime_b =~ m/(\d+):(\d+):(\d+)/;
			$hr_frac = ($2 / 60);
			$sec_frac = ($3 / 3600);
			$runtime_b_f = ($1 + $hr_frac + $sec_frac);
		}
		
		
		#runtimes have been established. now time to update scores
		#print "(end of ParseAndUpdate) output = $output \n";
		$output =~ m/score\s(.+)\s$/;
		my @scores = split / /, $1;
		$runHash{$runName}{$STORIcount}{"scores"} = [@scores];
		#print "runHash{$runName}{$STORIcount}{scores} = [" . join(" ",@scores) . "]\n";
		
		
		if (($runtime_a_f != -1) && ($runtime_b_f != -1)) {
			my $cpuTime = ($runtime_a_f + $runtime_b_f);
			#if ($runtime_a_f > $runtime_b_f) {
			#	$runHash{$runName}{$STORIcount}{"runtime"} = $runtime_a_f;
			#	print "runHash{$runName}{$STORIcount}{runtime}= $runtime_a_f \n";
			#}
			#else {
			#	$runHash{$runName}{$STORIcount}{"runtime"} = $runtime_b_f;
			#	print "runHash{$runName}{$STORIcount}{runtime}= $runtime_b_f \n";
			#}
			
			$runHash{$runName}{$STORIcount}{"runtime"} = $cpuTime;
			print "runHash{$runName}{$STORIcount}{runtime}= $cpuTime \n";
			
			
			if ($STORIcount < 2) {
				my $tempScore = Mean($runHash{$runName}{$STORIcount}{"scores"});			#when both sub-runs have completed each of their latest iterations, output the stats to a file
				my $t_runA_id = $runHash{$runName}{$STORIcount}{ids}->[0];
				my $t_runB_id = $runHash{$runName}{$STORIcount}{ids}->[1];
				#print jobStats "$runName\: $tempScore $t_runA_id $runtime_a_f $t_runB_id $runtime_b_f\n";
				ContinueRetrieval_t($runName, $run_params{$runName}{"parentDir"}, $run_params{$runName}{"taxaFile"}, $run_params{$runName}{"windowSize"}, $run_params{$runName}{"finalMaxFams"}, ($STORIcount+1));
				$not_converged{$runName} = 1;
			}
			elsif ($STORIcount >= 2) {
				my $a = $STORIcount;
				my $b = ($a - 1);
				my $c = ($b - 1);
				print "checking convergence... ";
				if (!(CheckConvergence2($runHash{$runName}{$a}{"scores"},$runHash{$runName}{$b}{"scores"},$runHash{$runName}{$c}{"scores"}))) {
					my $tempScore = Mean($runHash{$runName}{$STORIcount}{"scores"});			#when both sub-runs have completed each of their latest iterations, output the stats to a file
					my $t_runA_id = $runHash{$runName}{$STORIcount}{ids}->[0];
					my $t_runB_id = $runHash{$runName}{$STORIcount}{ids}->[1];
					#print jobStats "$runName\: $tempScore $t_runA_id $runtime_a_f $t_runB_id $runtime_b_f\n";
					if ($tempScore < 0.9) {
						ContinueRetrieval_t($runName, $run_params{$runName}{"parentDir"}, $run_params{$runName}{"taxaFile"}, $run_params{$runName}{"windowSize"}, $run_params{$runName}{"finalMaxFams"}, ($STORIcount+1));
					}
					elsif ($tempScore >= 0.9) {
						ContinueRetrieval_48hr($runName, $run_params{$runName}{"parentDir"}, $run_params{$runName}{"taxaFile"}, $run_params{$runName}{"windowSize"}, $run_params{$runName}{"finalMaxFams"}, ($STORIcount+1));
					}
					$not_converged{$runName} = 1;
				}
				else {
					if (exists $not_converged{$runName}) {
						delete $not_converged{$runName};
					}
					$converged{$runName} = 1;
				}
			}
		}
		else {
			$runHash{$runName}{$STORIcount}{"runtime"} = -1;
		}
		
	}
	}
	
	#print "(end ParseAndUpdate runHash)\n";
	#print Dumper \%runHash;
	
	SaveBackup($arrayFile);
	@runArr = @{LoadBackup($arrayFile)};

}


sub CheckConvergence2 {
	my $convergenceFlag=0;
	my $ref_a = shift(@_);
	my $ref_b = shift(@_);
	my $ref_c = shift(@_);
	my @arr1 = @{$ref_a};
	my @arr2 = @{$ref_b};
	my @arr3 = @{$ref_c};
	
	my $avg1 = Mean(\@arr1);
	my $avg2 = Mean(\@arr2);
	my $avg3 = Mean(\@arr3);
	
	$avg1 = sprintf("%.2f", $avg1);
	$avg2 = sprintf("%.2f", $avg2);
	$avg3 = sprintf("%.2f", $avg3);
	
	if (($avg1 > $avg2) && ($avg2 > $avg3)) {
		$convergenceFlag=0;
	}
	elsif (($avg1 > 0.9) && ($avg2 > 0.9) && ($avg3 > 0.9)) {
		if (abs($avg1 - $avg2) < 0.04) {
			$convergenceFlag=1;
		}
	}
	
	#print jobStats "Check Convergence: $avg1 $avg2 $avg3 ";
	
	#if ($convergenceFlag==1) {
	#	print jobStats "convergence achieved.\n";
	#}
	#else {
	#	print jobStats "nope\n";
	#}
	
	return $convergenceFlag;
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

sub LoadBackup {		#loads a specified CSV file into an array and returns the array's memory reference
	#print "(begin LoadBackup runHash)\n";
	#print Dumper \%runHash;
	
	my $file = shift(@_);
	
	if (-e $file) {
		open (arrFile, $file);
		unless (flock(arrFile, 1)) {
			warn "File $file already locked; waiting...\n";
			flock(arrFile, 1) or die;
		}
		
		my @tempArr = <arrFile>;
		close arrFile;
		my @parsedArr=();
		
		SET_LOOP: foreach my $line (@tempArr) {
			#print "parsing line: $line \n";
			if ($line =~ m/BEGIN\sPARAMETERS/) {
				last SET_LOOP;
			}
			chomp($line);
			my @row = split /\,/, $line;
			push @parsedArr, [@row];
		}
		
		#print "(end LoadBackup runHash)\n";
		#print Dumper \%runHash;
		
		return \@parsedArr;
	}
}

sub LoadParams {

	#print "(LoadParams runHash)\n";
	#print Dumper \%runHash;
	
	my %run_params=();
	my $file = shift(@_);
	my $readFlag=0;
	if (-e $file) {
		open (arrFile, $file);
		unless (flock(arrFile, 1)) {
			warn "File $file already locked; waiting...\n";
			flock(arrFile, 1) or die;
		}
		
		my @tempArr = <arrFile>;
		my @tempParamsArr=();
		close arrFile;
		my $r =0;
		foreach my $line (@tempArr) {
			chomp($line);
			if ($line =~ m/BEGIN\sPARAMETERS/) {
				$readFlag=1;
			}
			if ($readFlag==1) {
				my @row = split /\,/, $line;
				push @tempParamsArr, [@row];
				$r++;
			}
		}
		
		foreach my $ref (@tempParamsArr) {
			my @tempRow = @{$ref};
			my $runName = shift(@tempRow);
			foreach my $elt (@tempRow) {
				$run_params{$runName}{"parentDir"} = shift(@tempRow);
				$run_params{$runName}{"taxaFile"} = shift(@tempRow);
				$run_params{$runName}{"windowSize"} = shift(@tempRow);
				$run_params{$runName}{"finalMaxFams"} = shift(@tempRow);
			}
		}
	}
	return \%run_params;
}

sub SaveBackup {		#output runHash and run_params to a CSV file at specified loc, overwriting any previous contents
	my $file = shift(@_);
	
	open (backupFile, ">$file");		#should I use flock here? I don't know. I cant make sense of its documentation. As long as you don't run 2 instances of this program in the same folder, everything will be fine.
	
	#print "(begin SaveBackup runHash)\n";
	#print Dumper \%runHash;
	
	foreach my $runName (keys %runHash) {
		my @temp=(keys %{$runHash{$runName}});
		my $current = Max(\@temp);
		my @ids = @{$runHash{$runName}{$current}{"ids"}};
		my @STORIS_sorted = sort { $a <=> $b } @temp;  #sort ascending
		print backupFile $runName . "\," . join("\,", @ids) . "\,";
		foreach my $STORInum (@STORIS_sorted) {
			my @scores = @{$runHash{$runName}{$STORInum}{"scores"}};
			my $runtime = $runHash{$runName}{$STORInum}{"runtime"};
			print backupFile join("\,", @scores) . "\,r" . $runtime . ",";
		}
		if (exists $converged{$runName}) {
			print backupFile "converged,";
		}
		if (exists $paused{$runName}) {
			print backupFile "paused,";
		}
		print backupFile "\n";
	}
	print backupFile "BEGIN PARAMETERS\n";
	foreach my $runName (keys %run_params) {
		print backupFile $runName . "\," . $run_params{$runName}{"parentDir"} . "\," . $run_params{$runName}{"taxaFile"} . "\," . $run_params{$runName}{"windowSize"} . "\," . $run_params{$runName}{"finalMaxFams"} . "\n";
	}

	close backupFile;
	
	#print "(end SaveBackup runHash)\n";
	#print Dumper \%runHash;
}


sub Max {
	my($ref_arr) = @_;
	my @array = @{$ref_arr};
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@array);
	my $answer = $stat->max();
	return $answer;
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