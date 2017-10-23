#!usr/bin/perl -w

# Sequence_Deletion takes a string (e.g. biological sequence data like RNA/DNA/Protein) and makes deletion mutations. 
# Usage:
#
# $ perl Sequence_Deletion.pl input_sequence step_size indel_input > output.fasta
#
# Here, the input_sequence is whatever sequence you wish to manipulate. It must be a continuous string with no linebreaks. 
# The step_size defines the steps to take in making the serial "mutations" to the sequence: e.g. a value of  will make a change every 5 characters (nts)
# The indel_input defines what change to make. Numerical data (an integer value) will delete that number of residues from the sequence at every step.
# The output is in Fasta format where a title is given. The first sequence is the input_sequence, all mutants follow where the change is defined.    

my $InSequence = $ARGV[0];			# The input sequence
my $Step = $ARGV[1];				# The step size for making mutations 
my $Indel = $ARGV[2];				# Define the change size of the deletion to be made
chomp $InSequence;

my $Counter = 1;					# Set counter to keep track of mutants

#Print the input (e.g. native) sequence first
print ">Input sequence\n$InSequence\n";

#Scan through the input sequence using the defined step size
for ($i = 0; $i < ((length $InSequence) - $Indel); $i += $Step) {
    
	#Check if deletion (defined by a number). If true sequentially delete the fragment
    if ($Indel =~ m/[0-9]/) {
    my $OutSequence = $InSequence;
	substr ($OutSequence, $i, $Indel, "");
	print ">Mutant_$Counter deleted $Indel from position $i\n$OutSequence\n";
	
    $Counter += 1;
    }

}
