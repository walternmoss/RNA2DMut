use Cwd;
$dir = getcwd;

$infile = $ARGV[0];
$i = $ARGV[1];
$j = $ARGV[2];

open (INFILE, "$infile") || die "Can't open the infile!\n";
my @BigFasta = <INFILE>;
close(INFILE) || die "can't close file";

for (my $x = 0; $x < @BigFasta; $x += 2) {
    my $InHeader = $BigFasta[$x];
	chomp $InHeader;
    my $InSequence = $BigFasta[$x + 1]; 
	chomp $InSequence;
	my $Length = length $InSequence;
	if (($Length >= $i) and  ($Length <=  $j)) {print "$InHeader\n$InSequence\n"};
	}
