#!usr/bin/perl -w
use strict;
use Bio::SeqIO;
open IN,"$ARGV[0]"or die "$!";
my @array;
while(<IN>){
	chomp;
	push @array,$_;
	print "$_\n";
}
close IN;
my $ins = Bio::SeqIO -> new(-file => $ARGV[1], -format => "fasta");
my $outs = Bio::SeqIO -> new(-file => ">".$ARGV[2], -format => "fasta");
while(my $seq = $ins -> next_seq()){
	#my $desc=$seq->desc;
	#my ($id)=$desc =~ /gene:(.*)/;
	my $id=$seq->id;
	print "$id\n";
	if(grep{$id =~ /$_/}@array){
	 	$outs -> write_seq($seq);
	}
}
