# RNA2DMut
RNA 2D structure mutation and analysis tools

This repository consists of 4 perl scripts used to analyze RNA 2D structure mutations:

1. RNA2DMut.pl - given an input sequence, this program will generate all possible point mutant at every position, then evaluate their predicted effect on folding. 
2. Sequence_Evaluation.pl - given an input fasta file with a sequence, or set of sequences, this program will evaluate their structural properties. 
4. Sequential_Substitution.pl - given an input sequence, a step size and replacement string, this program will substitute the replacement string at every step position on the sequence. 
5. Sequence_Deletion.pl - given an input sequence, a step size and deletion length, this program will delete sequences of the defined length at every step position on the sequence.
6. Sequence_Insertion.pl - given an input sequence, a step size and insertion string, this program will instert the insertion string at every step position on the sequence.

RNA2DMut.pl and Sequence_Evaluation.pl both require that the The ViennaRNA Package be installed. This set of programs is available here: https://www.tbi.univie.ac.at/RNA/

All scripts are written in Perl. Available here: https://www.perl.org/get.html

Instructions on running each script are annotated in the source files. 
