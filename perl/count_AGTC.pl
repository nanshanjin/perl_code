use strict;
use Bio::SeqIO;

my $infastafile = $ARGV[0];
my $is = Bio::SeqIO -> new(-format=>'fasta', -file => $infastafile);
open OUT,">result.txt";

my($count_A,$count_T,$count_G,$count_C);
while(my $seqobj = $is -> next_seq()){
	my $id = $seqobj -> display_id();
	my $seq =  $seqobj ->seq;
	$count_A=($seq=~tr/A//);  
	$count_T=($seq=~tr/T//);  
	$count_G=($seq=~tr/G//);  
	$count_C=($seq=~tr/C//);  
	print OUT"$id -> $count_A\t$count_T\t$count_G\t$count_C\n";
}

