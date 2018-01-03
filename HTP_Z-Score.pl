#!usr/bin/perl -w

#Define filename
my $targetfile = $ARGV[0];
my $Randomizations = $ARGV[1];

open (INFILE, "$targetfile") || die "Can't open the infile!\n";
my @Fasta = <INFILE>;
close(INFILE) || die "can't close file";

my $N = @Fasta;


#Put FASTA seqs into a Hash

my %FastaHash = ();
my $CurrentSeq = "";

#Create array of names to keep track of order
my @SeqNames = ();

for ($i = 0; $i < $N; $i++) {

    $Line = $Fasta[$i];
    chomp $Line;

    #Place title and seq into their respective places in the FastaHash
	
	#Test if title
    if (substr($Line, 0, 1) eq ">") {

        my $SeqName = $Line;
        chomp $SeqName;
		$SeqName =~ s/\R//g;
        $SeqName =~ s/>//g;
        $SeqName =~ s/:/-/g;
        $CurrentSeq = $SeqName;
        push (@SeqNames, $CurrentSeq);
        }
    
	#Test if sequence data (must be IUPAC code)
    if ( $Line =~ m/(A|G|C|U|T|Y|R|K|M|B|D|H|V|N|S|W)/g && substr($Line, 0, 1) ne ">") {

    $FastaHash{$CurrentSeq} .= $Line;

    }
}

foreach my $SequenceName (@SeqNames) {

    chomp $SequenceName;   
    #print "$SequenceName\t";
	
    my $Sequence = $FastaHash{$SequenceName};
    chomp $Sequence;
    $Sequence =~ s/\s+//g;
    $Sequence =~ s/-//g;
    $Sequence = uc $Sequence;
	
	#Get nt frequencies
    my $NTFreqs = &NucFreqs($Sequence);	
	
	#Generate prescrambled control sequence
    my $ScrSequence = Scramble1X ($Sequence);

    #Get energy and ensemble diversity of the native sequence
	my $Command = "echo " . $Sequence . " | RNAfold -p --noPS";
	my @Out = `$Command`;
    my $STR_EN =  $Out[1];
    my @STR_EN = split(/\s+\(/, $STR_EN);
    my $EN = $STR_EN[1];
	my $STR = $STR_EN[0];
    $EN =~ s/\(//g;
    $EN =~ s/\)//g;	
    $EN =~ s/\R//g;
	chomp $EN;

	my $Centroid_Data = $Out[3];
	my @Centroid_Data = split(/\s+/, $Centroid_Data);
    my $Centroid = $Centroid_Data[0]; 	
	my $ED_Data = $Out[4];
	my @ED_Data = split(/\s+/, $ED_Data); 
	my $ED = $ED_Data[10];
	
	#print  "$EN\t$ED\t";
	
    #Get energy and ensemble diversity of the scrambled control sequence
	my $ScrCommand = "echo " . $ScrSequence . " | RNAfold -p --noPS";
	my @ScrOut = `$ScrCommand`;
    my $ScrSTR_EN =  $ScrOut[1];
    my @ScrSTR_EN = split(/\s+\(/, $ScrSTR_EN);
    my $ScrEN = $ScrSTR_EN[1];
    $ScrEN =~ s/\(//g;
    $ScrEN =~ s/\)//g;	
	
	my $ScrED_Data = $ScrOut[4];
	my @ScrED_Data = split(/\s+/, $ScrED_Data); 
	my $ScrED = $ScrED_Data[10];
	$ScrEN =~ s/\R//g;
	chomp $ScrEN;
	
	#Scamble the native sequence and get energies and ensemble diversities
    my @ScrambledSeqs = &Scramble($Sequence);	
	my @DGArray = &Energy(@ScrambledSeqs);
	my @EDArray = &EnsembleDiversity(@ScrambledSeqs);
	
	#Generate copies of these arrays to add the prescrambled values to
	my @ScrDGArray = @DGArray;
	my @ScrEDArray = @EDArray;
	
	#foreach my $Energy (@DGArray) {print "EN: $Energy\n";}	
	#foreach my $EnsD (@EDArray) {print "ED: $EnsD\n";}	
	
	
	#Add the native dG to the scrambled energy array and calculate z-score and p-value
	unshift (@DGArray, $EN);		
	my $DGZscore = &ZScore(@DGArray);
	my $DGPValue = &PValue(@DGArray);
	
	#Add the native ED to the scrambled ED array and calculate z-score and p-value
	unshift (@EDArray, $ED);
	my $EDZscore = &ZScore(@EDArray);
	my $EDPValue = &PValue(@EDArray);


	#foreach my $Energy (@DGArray) {print "EN: $Energy\n";}	
	#foreach my $EnsD (@EDArray) {print "ED: $EnsD\n";}		

	#Add the prescrambled dG to the scrambled energy array and calculate z-score and p-value
	unshift (@ScrDGArray,$ScrEN);
	my $ScrDGZscore = &ZScore(@ScrDGArray);
	my $ScrDGPValue = &PValue(@ScrDGArray);
	
	#Add the prescrambled ED to the scrambled ED array and calculate z-score and p-value
	unshift (@ScrEDArray,$ScrED);
	my $ScrEDZscore = &ZScore(@ScrEDArray);
	my $ScrEDPValue = &PValue(@ScrEDArray);	  

    #Print the output	
    print "$SequenceName\t$Sequence\t$STR\t$EN\t$DGZscore\t$DGPValue\t$Centroid\t$ED\t$EDZscore\t$EDPValue\t$ScrEN\t$ScrDGZscore\t$ScrDGPValue\t$ScrED\t$ScrEDZscore\t$ScrEDPValue\t$NTFreqs\n"
}

`rm  *.ps`;

######Sub-routine to scramble RNAs################################
sub Scramble {

        my $InSeq = $_[0];
        my @Out = ();

        for (my $i = 0; $i < $Randomizations; $i++) {
                my $OutSeq = "";
                my @InSeq = split ("", $InSeq);

                while (@InSeq > 0) {
                        my $Rand = rand(@InSeq);
                        my $RandBase = splice(@InSeq, $Rand, 1);
                        $OutSeq .= $RandBase;
                        }
                push(@Out, $OutSeq);
        }
        return @Out;
}

######Sub-routine to scramble RNAs################################
sub Scramble1X {

    my $InSeq = $_[0];
    my $OutSeq = "";
    
	my @InSeq = split ("", $InSeq);

    while (@InSeq > 0) {
        my $Rand = rand(@InSeq);
        my $RandBase = splice(@InSeq, $Rand, 1);
        $OutSeq .= $RandBase;
        }

    return $OutSeq;
}

######Sub-routine to calculate MFEs using RNAfold#################
sub Energy {

    my @engarr = @_;
    my $k = 0;
    my @returnarray = ();

    foreach my $Sequence (@engarr) {

	my $Command = "echo " . $Sequence . " | RNAfold";
	my @Out = `$Command`;
    #print @Out;
    #print "\n";
    my $STR_EN =  $Out[1];
    my @STR_EN = split(/\s+\(/, $STR_EN);
    my $EN = $STR_EN[1];
    #print "$EN\n";
    $EN =~ s/\(//g;
    $EN =~ s/\)//g;
    $EN =~ s/\R//g;
	chomp $EN;    
    push (@returnarray, $EN);

}
return @returnarray;
@returnarray = ();
}

######Sub-routine to calculate Nucleotide frequencies#######
sub NucFreqs {

        my $InSeq = $_[0];
        $InSeq =~ s/T/U/g;

        $Gs = 0;
        $Cs = 0;
        $As = 0;
        $Us = 0;

        while ($InSeq =~ /G/g) {$Gs += 1;}
        while ($InSeq =~ /C/g) {$Cs += 1;}
        while ($InSeq =~ /A/g) {$As += 1;}
        while ($InSeq =~ /U/g) {$Us += 1;}
        
	my $length = $As + $Cs + $Gs + $Us;
        my $OutFreqs = "$As\t$Gs\t$Cs\t$Us";
        return  $OutFreqs;
}

######Sub-routine to calculate Z-scores########
sub ZScore {

	my @Array = @_;

	#Remove the first value from the array
	my $Value = shift (@Array); 
	
	my $Sum = 0;
	my $Count = @Array;
	
	foreach my $l (@Array){ $Sum += $l;}
	
	my $AverageScr = $Sum/$Count;


	my $Sigma = 0;
    my $SumSigma = 0;
	
	foreach my $m (@Array){

		$Sigma = ($m - $AverageScr)**2;
		$SumSigma += $Sigma;
}
	my $SD = sqrt($SumSigma/$Count);

	my $ZScore = (($SD ne 0) ? (($Value - $AverageScr)/$SD):"Undefined");
        
	my $Output = substr($ZScore, 0, 5);
	return $Output;
}

######Sub-routine to calculate P-value (fraction of scrambled dG < native dG########
sub PValue {

	my @arr = @_;
	my $BelowNative = 0;
	my $TotalCount = @arr;
	
	my $Native = $arr[0]; 
	foreach my $l (@arr) {
	    if ($l < $Native) {$BelowNative += 1;}	    
	}
	
	my $Fraction = ($BelowNative / $TotalCount);

	
	return $Fraction;
}

######Sub-routine to calculate EDs using RNAfold#################
sub EnsembleDiversity {

    my @engarr = @_;
    my $k = 0;
    my @returnarray = ();

    foreach my $Sequence (@engarr) {

    my $Command = "echo " . $Sequence . " | RNAfold -p";
	my @Out = `$Command`;
	my $Data = $Out[4];
	my @Data = split(/\s+/, $Data); 
	my $EnsembleDiversity = $Data[10];
    
    push (@returnarray, $EnsembleDiversity);
}

return @returnarray;
@returnarray = ();
}
