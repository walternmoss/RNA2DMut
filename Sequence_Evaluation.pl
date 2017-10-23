#!usr/bin/perl -w

# Sequence_Evaluation is a program for evaluating the effects of mutations on RNA structure and ensemble structural properties.
# This program was developed by Walter N. Moss at Iowa State Universtiy. 
#
# Usage:
#
# $ perl Sequence_Evaluation.pl input_fasta_file temperature > output_file
#
# If used in research please Cite:
#
#  
# R. Lorenz, S.H. Bernhart, C. Hoener zu Siederdissen, H. Tafer, C. Flamm, P.F. Stadler and I.L. Hofacker (2011), "ViennaRNA Package 2.0", Algorithms for Molecular Biology: 6:26
# J.S. McCaskill (1990), "The equilibrium partition function and base pair binding probabilities for RNA secondary structures", Biopolymers: 29, pp 1105-1119

#Read FASTA file and fold sequences.

#Define filename
$targetfile = $ARGV[0];

open (INFILE, "$targetfile") || die "Can't open the infile!\n";
my @Fasta = <INFILE>;
close(INFILE) || die "can't close file";

my $N = @Fasta;

#Define temperature 
$Temperature = $ARGV[1];
if (! defined $Temperature) { $Temperature = 37;}

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
    my $Sequence = $FastaHash{$SequenceName};
    chomp $Sequence;
    $Sequence =~ s/\s+//g;
    $Sequence =~ s/-//g;
    $Sequence = uc $Sequence;
    
	#Run RNAfold partition function calculation
	my $Command = "echo " . $Sequence . " | RNAfold -p -T " . $Temperature;
	my @Out = `$Command`;
	
	#Parse the MFE data 
	my $MFEData = $Out[1];
	my @MFE = split(/\s+\(/, $MFEData);
	$MFEStructure = $MFE[0];
	my $MFE = $MFE[1];
	chomp $MFE;
	$MFE =~ s/\(//g;
	$MFE =~ s/\)//g;	
	chomp $MFEStructure;
	
	#Parse the centroid folding data
	my $CentroidData = $Out[3];
	my @Centroid = split(/\s+/, $CentroidData);
	my $CentroidStructure = $Centroid[0];
	
	#Parse the data for the Ensemble Diversity
	my $Data = $Out[4];
	my @Data = split(/\s+/, $Data); 
	my $EnsembleDiversity = $Data[10];

	# Print results for all mutant sequences and structures

	open (OUT, ">>", "RNA2DMut_Evaluate_output.txt");
	print OUT "$SequenceName\t$Sequence\t$MFEStructure\t$MFE\t$CentroidStructure\t$EnsembleDiversity\n";
    close OUT;

    }
    
`rm *.ps`;
