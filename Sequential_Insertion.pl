#!usr/bin/perl -w

# Sequential_Insertion takes a string (e.g. biological sequence data like RNA/DNA/Protein) and adds insertions. 
# Usage:
#
# $ perl Sequential_Insertion.pl input_sequence step_size indel_input > output.fasta
#
# Here, the input_sequence is whatever sequence you wish to manipulate. It must be a continuous string with no linebreaks. 
# The step_size defines the steps to take in making the serial "mutations" to the sequence: e.g. a value of  will make a change every 5 characters (nts)
# The indel_input defines what change to make. When indel_input is a string of characters (e.g. nts) that string is inserted at the step.
# The output is in Fasta format where a title is given. The first sequence is the input_sequence, all mutants follow where the change is defined.    

my $InSequence = $ARGV[0];			# The input sequence
my $Step = $ARGV[1];				# The step size for making mutations 
my $Indel = $ARGV[2];				# Define the change size of the deletion to be made
chomp $InSequence;

my $Counter = 1;					# Set counter to keep track of mutants

#Print the input (e.g. native) sequence first
print ">Input sequence\n$InSequence\n";

#Scan through the input sequence using the defined step size
for ($i = 0; $i < (length $InSequence); $i += $Step) {
    
	#Check if insertion (defined by a string). If true sequentially insert the fragment
    if ($Indel =~ m/[A-Za-z]/) {
    my $OutSequence = $InSequence;
	substr ($OutSequence, $i, 0, "$Indel");
	print ">Mutant_$Counter added $Indel to position $i\n$OutSequence\n";
	
    $Counter += 1;
    }	

}
