use strict;
use Bio::SeqIO;

my $infastafile = $ARGV[0];
my $is = Bio::SeqIO -> new(-format=>'fasta', -file => $infastafile);

open OUT,">$ARGV[1]";

my %hash;
my %hash2;
while(my $seqobj = $is -> next_seq()){
	my $id=$seqobj -> display_id();
	my $seq=$seqobj -> seq;
	$hash2{$id}=$seq;
	my $lens=length($seq);
	$hash{$id}=$lens;
}
#my @keys=(sort { $hash{$b} <=> $hash{$a} } keys %hash)[0..5];#对valuse排序
#for(@keys){
#	print OUT">$_\n$hash2{$_}\n";
#}

foreach my $key ( (sort { $hash{$b} <=> $hash{$a} } keys %hash)[0..4]) {#对value排序，取最大的前六个值
			print OUT ">$key\n$hash2{$key}\n";
}








