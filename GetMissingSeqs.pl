#example call for testing
#perl /nv/hp10/jstern7/STORI/GetMissingSeqs.pl /nv/hp10/jstern7/scratch/STORI_test/L20_bacta/STORI_out_sc2_L20_bacta.txt /nv/hp10/jstern7/scratch/STORI_test/L20_bactb/STORI_out_sc2_L20_bactb.txt L20_bact_fina L20_bact_finb /nv/hp10/jstern7/scratch/STORI_test bacteria
#perl /nv/hp10/jstern7/STORI/GetMissingSeqs.pl /nv/hp10/jstern7/scratch/STORI_test/L8_euka/STORI_out_sc2_L8_euka.txt /nv/hp10/jstern7/scratch/STORI_test/L8_eukb/STORI_out_sc2_L8_eukb.txt L8_euk_fin_a L8_euk_fin_b /nv/hp10/jstern7/scratch/STORI_test eukaryota

#GetMissingSeqs

#input: two run output files, run names, taxa file
#output: two new output files with missing seqs filled in

#algorithm: 
# foreach phylum
# 	foreach family
# 		get the gis from the taxa that have them
#		if #full taxa > 2
# 			get the empty taxa & put in a short term list. also add to a longer-term hash of all taxa that are empty in this phylum.
# 			foreach empty taxon
# 				blast the full gis against it
#				if there is a best hit that is consistent >50% of the time, save it and the avg bitscore
# 	foreach taxon in the emptyTaxa_phylum hash
# 		determine whether there are conflicts in family assignments of missing seqs.
# 		if so, assign the missing seq to the family in which it has the highest avg bitscore
# make a final assignment table, combining the original input seqs with the missing seqs
# output to file in same format as the input


# my $cmd = "perl $GetMissingSeqspath $fileA $fileB $runNameA $runNameB $parentDir $taxaFile";

use lib '/nv/hp10/jstern7/perl5reinstall/lib';
use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';
use Data::Dumper;

my $blastdbDir="/nv/hp10/jstern7/scratch/universal120312/blast/";
my $STORIdir = "/nv/hp10/jstern7/STORI3";
my $blastdbcmdPath = "/nv/hp10/jstern7/STORI3/blastdbcmd";
my $blastpPath = "/nv/hp10/jstern7/STORI3/blastp"; 

my $inputA=shift(@ARGV);
my $inputB=shift(@ARGV);
my $runNumberA=shift(@ARGV);
my $runNumberB=shift(@ARGV);
my $parentDirPath=shift(@ARGV);
if ($parentDirPath =~ m/\/$/) {
	$parentDirPath =~ s/^(.+)\/$/$1/;
}
my $taxaFile=shift(@ARGV);



my @taxaArr=();
my $sc = 0;


my %fams_nophyla=();
my %fams_phyla=();
my %finalAssignments=();

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


my $cmd = "mkdir " . $sourceDirPathA;
#print "cmd is: " . $cmd . "\n";
system($cmd);
my $cmd = "mkdir " . $sourceDirPathB;
#print "cmd is: " . $cmd . "\n";
system($cmd);

my $cmd = "cp " . $STORIdir . "/taxa-master[" . $taxaFile . "].txt " .  $orderAfile;
#print "cmd is: " . $cmd . "\n";
system($cmd);
my $cmd = "cp " . $STORIdir . "/taxa-master[" . $taxaFile . "].txt " .  $orderBfile;
#print "cmd is: " . $cmd . "\n";
system($cmd);

my $outputFileA="$sourceDirPathA/STORI_out_sc" . $sc . "_" . $runNumberA . ".txt";
open (outA, ">$outputFileA");
my $outputFileB="$sourceDirPathB/STORI_out_sc" . $sc . "_" . $runNumberB . ".txt";
open (outB, ">$outputFileB");

foreach my $pair_ref (@pairScores_sorted) {
	if ($pair_ref->[2] > 0) {
		my $famA = $pair_ref->[0];
		my $famB = $pair_ref->[1];
		if (exists $famsA{$famA}) {
			if (exists $famsB{$famB}) {
				my %consensus = %{GetConsensusFam($famsA{$famA},$famsB{$famB})};
				
				my $consensusName = $famA . $famB;
				$fams_nophyla{$consensusName} = {%consensus};
				
				delete $famsA{$famA};
				delete $famsB{$famB};
			}
		}
	}
}



#read in the phyla file
my %phylaHash=%{FillPhylaHash($taxaFile)};

#print "here is fams_nophyla\n";
#print Dumper \%fams_nophyla;

#OutputTable(\%fams_nophyla);

#print "here is phylaHash\n";
#print Dumper \%phylaHash;

#sort the fams_nophyla into fams_phyla
foreach my $fam (keys %fams_nophyla) {
	foreach my $phylum (keys %phylaHash) {
		my %tempTaxa = %{$phylaHash{$phylum}};
		foreach my $tax (keys %tempTaxa) {
			$fams_phyla{$fam}{$phylum}{$tax} = $fams_nophyla{$fam}{$tax} if (exists $fams_nophyla{$fam}{$tax});
		}
	}
}

#print "here is fams_phyla\n";
#print Dumper \%fams_phyla;

foreach my $phylum (keys %phylaHash) {
	my %emptyTaxa_phylum=();
	my %missingSeqs_phylum=();
	foreach my $family (keys %fams_phyla) {
		my %family_missingSeqs=();
		my @fullGis = @{GetFullGis($fams_phyla{$family}{$phylum})};                         #get the gis from the taxa that have them
		if ($#fullGis >= 1) {
			my @emptyTaxa = @{GetEmptyTaxa($fams_phyla{$family}{$phylum}, $phylaHash{$phylum})};					#get the empty taxa & put in a short term list, just for this family.
			#print "emptyTaxa\= " . join(" ", @emptyTaxa) . "\n";
			#foreach my $taxon (keys %emptyTaxa_phylum) {
			foreach my $taxon (@emptyTaxa) {
				$emptyTaxa_phylum{$taxon} = 1;												#also add the empty taxa to a longer-term hash of all taxa that are empty, for any family, in this phylum.
				$family_missingSeqs{$taxon} = [@{GetConsensusBestHit($fams_phyla{$family}{$phylum}, $taxon)}];	#blast the full gis against the empty taxon. if there is a best hit that is consistent >50% of the time, save it and the avg bitscore
				
				#print "family_missingSeqs\{$taxon\}\= " . join(" ", @{$family_missingSeqs{$taxon}}) . "\n";
			}
		}
		#print "family_missingSeqs\n";
		#print Dumper \%family_missingSeqs;
		$missingSeqs_phylum{$family} = {%family_missingSeqs};
		#print "missingSeqs_phylum\n";
		#print Dumper \%missingSeqs_phylum;
	}
	my %duplicateAssignments=();
	my %temp_taxgi=();
	foreach my $taxon (keys %emptyTaxa_phylum) {											#determine whether there are conflicts in family assignments of missing seqs.
		foreach my $family (keys %fams_phyla) {
			if (exists $duplicateAssignments{$missingSeqs_phylum{$family}{$taxon}}) {
				my @temp = @{$duplicateAssignments{$missingSeqs_phylum{$family}{$taxon}->[0]}};
				push @temp, $family;
				$duplicateAssignments{$missingSeqs_phylum{$family}{$taxon}->[0]} = [@temp];
			}
			else {
				my @temp = ($family);
				my $gi = $missingSeqs_phylum{$family}{$taxon}->[0];
				$duplicateAssignments{$gi} = [@temp] if ($gi > 0);
				$temp_taxgi{$gi} = $taxon if ($gi > 0);				#keep track of the parent taxa
			}
		}
	}
	
	foreach my $gi (keys %duplicateAssignments) {												#determine whether there are conflicts in family assignments of missing seqs. if so, assign the missing seq to the family in which it has the highest avg bitscore
		my @famAssignments = @{$duplicateAssignments{$gi}};
		my $famChoice="";
		my $taxon = $temp_taxgi{$gi};
		if ($#famAssignments > 0) {
			my @scoredFams=();
			foreach my $fam (@famAssignments) {
				my @famAndScore = ($fam, $missingSeqs_phylum{$fam}{$taxon}->[1]);
				push @scoredFams, [@famAndScore];
			}
			my @scoredFams_sorted = sort { $b->[1] <=> $a->[1] } @scoredFams;  #sort descending
			$famChoice = $scoredFams[0][0];
		}
		else {
			$famChoice = $famAssignments[0];
		}
		$finalAssignments{$famChoice}{$taxon} = $gi if (NotPresentInOriginal($taxon, $gi));											#make a final assignment table, combining the original input seqs with the missing seqs
	}
	

}

	#print "here is finalAssignments just missing\n";
	#print Dumper \%finalAssignments;

foreach my $fam (keys %fams_nophyla) {															#make a final assignment table, combining the original input seqs with the missing seqs
	my %temp= %{$fams_nophyla{$fam}};
	foreach my $tax (keys %temp) {
		$finalAssignments{$fam}{$tax} = $fams_nophyla{$fam}{$tax};
	}
}

#print "here is finalAssignments done\n";
#print Dumper \%finalAssignments;

OutputTable(\%finalAssignments);

close outA;
close outB;



sub NotPresentInOriginal {
	my $taxon = shift(@_);
	my $gi = shift(@_);
	my $ans = 1;
	
	foreach my $fam (keys %fams_nophyla) {
		if ($fams_nophyla{$fam}{$taxon} == $gi) {
			$ans = 0;
		}
	}
	
	return $ans;
}

sub OutputTable {
	my $table_ref = shift(@_);
	my %table = %{$table_ref};
	print outA "\n\n\n final table main\n\n";
	print outB "\n\n\n final table main\n\n";
	my @familyNames = keys %table;
	print outA "taxon ";
	print outB "taxon ";
	foreach my $family (@familyNames) {
		print outA $family . " ";
		print outB $family . " "; }
	print outA "\n";
	print outB "\n";
	foreach my $taxon (@taxaArr) {
		print outA $taxon . "  ";
		print outB $taxon . "  ";
		foreach my $family (@familyNames) {
			if (exists $table{$family}{$taxon}) {
				print outB $table{$family}{$taxon} . " ";
				print outA $table{$family}{$taxon} . " ";}
			else {
				print outA "-1 ";
				print outB "-1 ";}
		}
		print outA "\n";
		print outB "\n";
	}
}



sub GetFullGis {
	my $ref = shift(@_);
	my %txgi = %{$ref};
	my @gis= values %txgi;
	return \@gis;
}


sub GetEmptyTaxa {
	my $ref = shift(@_);
	my %txgi = %{$ref};
	my $ref = shift(@_);
	my %txgi_all = %{$ref};
	my @taxa=();
	my @allTaxa = keys %txgi_all;
	
	foreach my $taxon (@allTaxa) {
		if (!(exists $txgi{$taxon})) {
			push @taxa, $taxon;
		}
	}
	
	return \@taxa;
}


sub GetConsensusBestHit {
	my $ref = shift(@_);
	my %txgi = %{$ref};
	my $taxon = shift(@_);
	my $ans=-1;
	
	#print "GetConsensusBestHit  taxon is $taxon    here is txgi";
	#print Dumper \%txgi;
	
	my $queryFastaPath = $sourceDirPathA . "/temp-query.fasta";
	my $blastResultsPath = $sourceDirPathA . "/temp-blast.results";
	
	MakeQueryFasta(\%txgi, $queryFastaPath);
	#print "calling Blast\n";
	Blast($queryFastaPath, $taxon, $blastResultsPath);
	
	$ans = ParseAndAssign($blastResultsPath);
	my @temp=@{$ans};
	#print "answer\= \(" . join(" ", @temp) . "\)\n";
		
	return \@temp;
}




sub MakeQueryFasta {
	
	my $hashref=shift(@_);
	my %txgi = %{$hashref};
	
	#print out "(begin MakeQueryFasta) here is txgi\:\n";
	#print out Dumper \%txgi;
	
	my $queryFastaPath = shift(@_);
	
	my $cmd = "touch $queryFastaPath";
	system($cmd);
	my $cmd = "rm $queryFastaPath";
	system($cmd);
	
	my %lookup=();
	
	foreach my $taxon (keys %txgi) {
		if ($txgi{$taxon} > -1) {
			$lookup{$txgi{$taxon}} = $taxon;
		}
	}
	
	open (fastaFile, ">>$queryFastaPath");
	
	foreach my $gi (keys %lookup) {
		my $taxon = $lookup{$gi};
		if ($gi > 0) {
			my $cmd="$blastdbcmdPath -entry $gi -db " . $blastdbDir . "/$taxon -outfmt \"\%f\"";
			#print out "cmd is: " . $cmd . "\n";
			my $seq = qx($cmd);
			#print out $seq;
			print fastaFile $seq;
		}
	}
	
	close fastaFile;
	#print out "(MakeQueryFasta)\n";
}

sub Blast {								#Hopefully the result is ortholog assignments that are not too stringent and not too permissive.
	my $queryFastaPath = shift(@_);
	my $taxon = shift(@_);
	my $blastResultsPath = shift(@_);
	my $cmd = $blastpPath . " -db " . $blastdbDir . "/$taxon -query $queryFastaPath -evalue 0.00000005 -outfmt \"6 qgi sgi bitscore\" -num_descriptions 1 -num_alignments 1 -parse_deflines";
	#print "cmd is: $cmd\n";
	my $output = qx($cmd);
	open (bres, ">$blastResultsPath");
	print bres $output;
	#print "the output is: " . $output . "\n";
	close bres;
}

sub ParseAndAssign {
	my $blastResultsPath=shift(@_);
	
	open (results, $blastResultsPath);
	
	my %blastHash=();
	my $totalSeqs=0;
	
	while (<results>) {
		my @tmp = split;
		my $query = $tmp[0];
		my $subject = $tmp[1];
		my $bitscore = $tmp[2];
		if (!(exists $blastHash{$subject})) {
			my @temp = (1, $bitscore);
			$blastHash{$subject} = [@temp];
		}
		else {
			my @temp = @{$blastHash{$subject}};
			$temp[0]++;
			$temp[1] += $bitscore;
			$blastHash{$subject} = [@temp];
		}
		$totalSeqs++;
	}
	close results;
	
	#print "\(ParseAndAssign\) here is blastHash \n";
	#print Dumper \%blastHash;
	my @blastHash_u=();
	
	foreach my $subj (keys %blastHash) {
		my @temp = ($subj, $blastHash{$subj}->[0], $blastHash{$subj}->[1]);
		push @blastHash_u, [@temp];
	}
	
	my @blastHash_s = sort { $b->[1] <=> $a->[1] } @blastHash_u;  #sort @blastHash_u descending by the freq column
	
	#print "\(ParseAndAssign\) here is blastHash_s \n";
	#print Dumper @blastHash_s;
	
	my $vote = 0;
	$vote = (($blastHash_s[0][1]) / $totalSeqs) if (($totalSeqs >= 1) && ($totalSeqs <= 3));
	$vote = (($blastHash_s[0][1]) / 3) if (($totalSeqs >= 4));
	my @ans = (-1, -1);
	
	#print "\(ParseAndAssign\) vote $vote \n";
	if ($vote > 0.5) {	
		$ans[0] = $blastHash_s[0][0];
		$ans[1] = ($blastHash_s[0][2] / $blastHash_s[0][1]);
	}
	#print "\(ParseAndAssign\) ans " . join(" ", @ans) . "\n";
	return \@ans;
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


sub FillPhylaHash {
	my %phylaHash=();
	my $taxaFile = shift(@_);
	my $phylaFilePath = "";
		$phylaFilePath = $STORIdir . "/" . $taxaFile . "\.txt";
		
	if (-e $phylaFilePath) {
		open phyla, $phylaFilePath;
		#print "opened $phylaFilePath\n";
	}
	else {
		print "could not find phyla file $phylaFilePath - program will not work.\n";
		return 0;
	}
	my @phylaArr=<phyla>;
	#print $#phylaArr;
	close phyla;
	my $count=0;
	foreach my $line (@phylaArr)
	{
		chomp($line);  #removes any newline characters at the end of the string
		$phylaArr[$count] = $line;
		#print $phylaArr[$count] . "\n";
		$count++;
	}
	foreach my $line (@phylaArr) {
		my @temp = split /\s/, $line;
		my $taxid = $temp[0];
		my $phylum = $temp[2];
		$phylaHash{$phylum}{$taxid} = 1;
		push @taxaArr, $taxid;
	}
	return \%phylaHash;
}