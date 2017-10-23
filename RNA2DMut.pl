#!usr/bin/perl -w

# RNA2DMut is a program for introducing all possible point mutations into an RNA
# and evaluating their effects on RNA structure and ensemble structural properties.
# This program was developed by Walter N. Moss at Iowa State Universtiy. 
#
# Usage:
#
# $ perl RNA2DMut.pl input_sequence constraints bracket_structure temperature > output_file
#
# If used in research please Cite:
#
#  
# R. Lorenz, S.H. Bernhart, C. Hoener zu Siederdissen, H. Tafer, C. Flamm, P.F. Stadler and I.L. Hofacker (2011), "ViennaRNA Package 2.0", Algorithms for Molecular Biology: 6:26
# J.S. McCaskill (1990), "The equilibrium partition function and base pair binding probabilities for RNA secondary structures", Biopolymers: 29, pp 1105-1119

my $InSequence = $ARGV[0];


my @Mutants = Mutagenize ($InSequence);
my @Results = ();

#Create outfiles and print header info
open (OUT1, ">>", "RNA2DMut_output_1.txt");
print OUT1 "\tSequence\tMFE_structure\tdelta(G)\tCentroid_structure\tED\n";
close OUT1;

open (OUT2, ">>", "RNA2DMut_output_2.txt");	
print OUT2 "nt\tWT\tWT_ED\tMut\tMut\tMut_ED\tMut\tMut\tMut_ED\tMut\tMut\tMut_ED\n";
close OUT2;

for (my $i=0; $i < @Mutants; $i += 1) {
    
	my $Sequence = $Mutants[$i];
    chomp $Sequence;
    $Temperature = $ARGV[3];	
    if (! defined $Temperature) { $Temperature = 37;}
	
	my $Command = "echo " . $Sequence . " | RNAfold -p -T " . $Temperature;
	my @Out = `$Command`;
	my $Data = $Out[4];
	my @Data = split(/\s+/, $Data); 
	my $EnsembleDiversity = $Data[10];
	#my $FracMFE = $Data[7];
	#$FracMFE =~ s/\;//;
    my $MFEData = $Out[1];
	my @MFE = split(/\s+\(/, $MFEData);
	$MFEStructure = $MFE[0];
	my $MFE = $MFE[1];
	chomp $MFE;
	$MFE =~ s/\(//g;
	$MFE =~ s/\)//g;	
	chomp $MFEStructure;
	
	my $CentroidData = $Out[3];
	my @Centroid = split(/\s+/, $CentroidData);
	my $CentroidStructure = $Centroid[0];

	# Print results for all mutant sequences and structures
	open (OUT1, ">>", "RNA2DMut_output_1.txt");
	print OUT1 "Mutant_$i\t$Sequence\t$MFEStructure\t$MFE\t$CentroidStructure\t$EnsembleDiversity\n";
    close OUT1;
	
   	my $Result = "$Sequence\t$MFE\t$MFEStructure\t$EnsembleDiversity\t$CentroidStructure";
    push (@Results, $Result);	

        }

`rm *.ps`;

#Extract WT Results to compare mutants to
my $WTResult = shift (@Results);
my @SplitWTResult = split (/\t/, $WTResult);
my $WTSequence = $SplitWTResult[0];
my $WTStructure = $SplitWTResult[2];
my $WTED = $SplitWTResult[3];	
my $WTCentroid = $SplitWTResult[4];	

#Go through results and find the most perturbing mutation at each residue

#Store diffED scores for VARNA colormap
my $ColorMap = "";

my $Counter = 1;
#Loop through every position in mutants 
for (my $x = 0; $x < (length $InSequence); $x += 1) {
    $Position = $x + 1;
    #Grab the WT base and initialize variables to store best result

	my $WTBase = substr ($WTSequence, $x, 1);	
    my $BestMutBase = $WTBase;
	my $MaxMutED = $WTED;
	my $MinMutED = $WTED;
	
	my $MutResults = "$Position\t$WTBase\t$WTED";
	
    #Loop through every mutant sequence to find the most perturbing mutant at every position: comparing to the WT Mutant_0 sequence.

	
	foreach my $MutResult (@Results) {
	
	    my @SplitMutResult = split (/\t/, $MutResult);
	    my $MutantSequence = $SplitMutResult[0];
        my $MutantED = $SplitMutResult[3];
		
		#If mutation present check ED and get max ED for that site
		my $MutantBase = substr ($MutantSequence, $x, 1);
							
		if ($MutantBase !~ m/$WTBase/) {

            $MutResults .= "\t$MutantBase\tMutant_$Counter\t$MutantED";
			$Counter += 1;
		
		    if ($MutantED > $MaxMutED) { 			
			    $MaxMutED = $MutantED;
            }					
		    if ($MutantED < $MinMutED) { 			
			    $MinMutED = $MutantED;
            }	
		}
    }
    
    my $diffMaxED = $MaxMutED - $WTED;
    my $diffMinED = $WTED - $MinMutED;	
	# print All mutations
	open (OUT2, ">>", "RNA2DMut_output_2.txt");	
    print OUT2 "$MutResults\n";
	close OUT2;	
	
	$ColorMap_MaxED .= "$diffMaxED,";
	$ColorMap_MinED .= "$diffMinED,";
}

#Print output for VARNA 
#print $ColorMap;

my $OutStructure = $ARGV[2];
if (! defined $OutStructure) { $OutStructure = $WTCentroid;} 

my $VARNA_Script_MaxED = "java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd  -algorithm radiate -baseOutline" . ' "#FFFFFF" -sequenceDBN "' . $WTSequence . '" -structureDBN "' . $OutStructure . '" -colorMapStyle red -colorMap "' . $ColorMap_MaxED . '" -o RNA2DMut_MaxED_Fig.eps';
my $VARNA_Out_MaxED = "java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd  -algorithm radiate -baseOutline" . ' "#FFFFFF" -sequenceDBN "' . $WTSequence . '" -structureDBN "' . $OutStructure . '" -colorMapStyle red -colorMap "' . $ColorMap_MaxED . '"';

my $VARNA_Script_MinED = "java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd  -algorithm radiate -baseOutline" . ' "#FFFFFF" -sequenceDBN "' . $WTSequence . '" -structureDBN "' . $OutStructure . '" -colorMapStyle blue -colorMap "' . $ColorMap_MinED . '" -o RNA2DMut_MinED_Fig.eps';
my $VARNA_Out_MinED = "java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd  -algorithm radiate -baseOutline" . ' "#FFFFFF" -sequenceDBN "' . $WTSequence . '" -structureDBN "' . $OutStructure . '" -colorMapStyle blue -colorMap "' . $ColorMap_MinED . '"';


open (OUT3, ">", "VARNA_output_maxED.txt");
print OUT3 $VARNA_Out_MaxED;
close OUT3;

open (OUT4, ">", "VARNA_output_minED.txt");
print OUT4 $VARNA_Out_MinED;
close OUT4;

# Run VARNA script to generate figure

system $VARNA_Script_MaxED;
system $VARNA_Script_MinED;
############################################################################################################################

sub Mutagenize {

    my $InSeq = $_[0];
	my @OutaRRAY = ();
	chomp $InSeq;
	$InSeq =~ s/\R//g;
	chomp $InSeq;
	my $Mask = $ARGV[1];
	
	# If user does not enter mask, mutate all bases
	if (! defined $Mask) { $Mask .= "." x (length $InSeq);} 	
	my @Mask = split ("", $Mask);
	
	push (@OutArray, $InSeq); 
	
	
	for ($i = 0; $i < (length $InSeq); $i += 1) {
	
	    # Test for mask constraint: "x" means to skip mutating (freezes WT base), "." means no constraint, 
		# "Y" means to only mutate to pYrimidines, and "R" means to only mutate to puRines.		
		
		my $Constraint = $Mask[$i];
		
		if ($Constraint =~ m/x/) {next;}
		
	    if ($Constraint =~ m/\./) {
		
			my $Seq = $InSeq;
	    	my $Base = substr ($Seq, $i, 1);

			if ($Base =~ m/A/) { 
		    	substr($Seq, $i, 1) = G; 
				push (@OutArray, ($Seq)); 
		    	substr($Seq, $i, 1) = C; 
				push (@OutArray, ($Seq));
		    	substr($Seq, $i, 1) = U; 
				push (@OutArray, ($Seq));
            	next;			
        	}			
			if ($Base =~ m/G/) { 
		    	substr($Seq, $i, 1) = A; 
				push (@OutArray, ($Seq)); 
		    	substr($Seq, $i, 1) = C; 
				push (@OutArray, ($Seq));
		    	substr($Seq, $i, 1) = U; 
				push (@OutArray, ($Seq));
            	next;			
        	}	
			if ($Base =~ m/C/) { 
		    	substr($Seq, $i, 1) = A; 
				push (@OutArray, ($Seq)); 
		    	substr($Seq, $i, 1) = G; 
				push (@OutArray, ($Seq));
		    	substr($Seq, $i, 1) = U; 
				push (@OutArray, ($Seq));
            	next;			
			}
			if ($Base =~ m/U/) { 
		    	substr($Seq, $i, 1) = A; 
				push (@OutArray, ($Seq)); 
		    	substr($Seq, $i, 1) = C; 
				push (@OutArray, ($Seq));
		  	    substr($Seq, $i, 1) = G; 
				push (@OutArray, ($Seq));
                next;			
            }		
		}
		
	    if ($Constraint =~ m/Y/) {
		
	        my $Seq = $InSeq;
	        my $Base = substr ($Seq, $i, 1);

		    if ($Base =~ m/A/) { 

		        substr($Seq, $i, 1) = C; 
			    push (@OutArray, ($Seq));
		        substr($Seq, $i, 1) = U; 
			    push (@OutArray, ($Seq));
                next;			
            }			
		    if ($Base =~ m/G/) { 

		        substr($Seq, $i, 1) = C; 
		    	push (@OutArray, ($Seq));
		        substr($Seq, $i, 1) = U; 
		    	push (@OutArray, ($Seq));
                next;			
            }	
		    if ($Base =~ m/C/) { 
		        substr($Seq, $i, 1) = U; 
		    	push (@OutArray, ($Seq));
                next;			
		    }
		    if ($Base =~ m/U/) { 
		        substr($Seq, $i, 1) = C; 
		    	push (@OutArray, ($Seq));
                next;			
            }				
		}
		
	    if ($Constraint =~ m/R/) {
		
	        my $Seq = $InSeq;
	        my $Base = substr ($Seq, $i, 1);

		    if ($Base =~ m/A/) { 
		        substr($Seq, $i, 1) = G; 
		    	push (@OutArray, ($Seq)); 
                next;			
            }			
		    if ($Base =~ m/G/) { 
		        substr($Seq, $i, 1) = A; 
			    push (@OutArray, ($Seq)); 
                next;			
            }	
		    if ($Base =~ m/C/) { 
		        substr($Seq, $i, 1) = A; 
			    push (@OutArray, ($Seq)); 
		        substr($Seq, $i, 1) = G; 
			    push (@OutArray, ($Seq)); 
                next;			
		    }
		    if ($Base =~ m/U/) { 
		        substr($Seq, $i, 1) = A; 
			    push (@OutArray, ($Seq)); 
		        substr($Seq, $i, 1) = G; 
			    push (@OutArray, ($Seq));
                next;			
            }				
		}
    }	
 
return @OutArray; 
		
}	
