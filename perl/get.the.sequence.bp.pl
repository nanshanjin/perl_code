use strict;
use Bio::SeqIO;

my %hash;my %hash2;


my $infastafile = $ARGV[0];
my $is = Bio::SeqIO -> new(-format=>'fasta', -file => $infastafile);
open GFF,"<$ARGV[1]";
open OUT,">$ARGV[2]";

while(my $seqobj = $is -> next_seq())
{
   my $id = $seqobj -> display_id();
#   print $id,"\n";
   my $seq = $seqobj->seq;
   $seq=~s/\s+//g; $seq=~s/\r//g;$seq=~s/\n//g;
   $hash{$id}=$seq;
}

$/="\n";
	while (<GFF>) 
	{
		chomp; 
        next if /^Chrom/;
		next if ($_=~m/^$/);
		my ($chr,$pos,$ref,@others)=split/\t/,$_;
        my $pos1=$pos-300-1;
        my $pos2=$pos+300;
        my $tmp=substr($hash{$chr},$pos1,600);
        print OUT "$chr\t$pos\t$ref\t$tmp\n";
	}
close GFF;

