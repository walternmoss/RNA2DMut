#!usr/bin/perl -w

# Sequential_Substitution takes a string (e.g. biological sequence data like RNA/DNA/Protein) and sequentially substitutes a sequence fragment. 
# Usage:
#
# $ perl Sequential_Substitution.pl input_sequence step_size substituted_sequence > output.fasta
#
# Here, the input_sequence is whatever sequence you wish to manipulate. It must be a continuous string with no linebreaks. 
# The step_size defines the steps to take in making the serial "mutations" to the sequence: e.g. a value of 5 will make a change every 5 characters (nts)
# The substituted_sequence defines what sequence to substitute at each position. 
# The output is in Fasta format where a title is given. The first sequence is the input_sequence, all mutants follow where the change is defined.    

my $InSequence = $ARGV[0];			# The input sequence
my $Step = $ARGV[1];				# The step size for making mutations 
my $SubSeq = $ARGV[2];				# Define the sequence to substitute 
chomp $InSequence;
my $SubLen = length $SubSeq;		# Get subseq length to replace fragment in InSequence

my $Counter = 1;					# Set counter to keep track of mutants

#Print the input (e.g. native) sequence first
print ">Input sequence\n$InSequence\n";

#Scan through the input sequence using the defined step size
for ($i = 0; $i < (length $InSequence); $i += $Step) {
    
	my $OutSequence = $InSequence;
	substr ($OutSequence, $i, $SubLen, $SubSeq);
	print ">Mutant_$Counter substituted $SubSeq at position $i\n$OutSequence\n";
	
    $Counter += 1;
    }


