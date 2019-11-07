use strict;
use Bio::SeqIO;

my $cdsfile = $ARGV[0];

my $gcfile  = $ARGV[1];
open(OUT,">$gcfile") or die "can not open $gcfile due to $!\n ";

my $stream = Bio::SeqIO -> new(-file => $cdsfile, '-format' => 'fasta');

while(my $seqobj = $stream -> next_seq)
{
   my @seq =split(//, uc($seqobj->seq));#uc():change a to A;
   my $id = $seqobj -> display_id();

   my ($gc, $gc1, $gc2, $gc3) = (0, 0, 0, 0);
   for(my $i = 0; $i <=$#seq-3; $i+=3)
   {
     if($seq[$i] eq "G" || $seq[$i] eq "C")
     {$gc1 ++; $gc++;}
     if($seq[$i+1] eq "G" || $seq[$i+1] eq "C")
     {$gc2 ++; $gc++;}
     if($seq[$i+2] eq "G" || $seq[$i+2] eq "C")
     {$gc3 ++; $gc++;}
   }
   #print $id."\n";
if($#seq <1){next;}   
   $gc1 = 3*$gc1/($#seq + 1);
   $gc2 = 3*$gc2/($#seq + 1);
   $gc3 = 3*$gc3/($#seq + 1);
   $gc  = $gc/($#seq+1);

   my $tmpstr = $gc." ".$gc1." ".$gc2." ".$gc3;
   print OUT "$id\t$tmpstr\n";
}
