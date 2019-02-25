use strict;

use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Seq;

use Bio::Align::Utilities qw(aa_to_dna_aln);

my $idfile = $ARGV[0];
open(ID,$idfile) or die "can not open id file $idfile due to $!.\n";

my %seqid;

while(<ID>)
{
	$_ =~s/[\n\r]//g;
	my @array = split(/\t/,$_);
	$seqid{$array[1]} = $array[0];
}

     my $consensustree =$ARGV[1];

#### recover the original ids in the tree files
     open(TRI, $consensustree) or die"cannot open $consensustree due to $!\n";
     open(TRO, ">".$consensustree.".nwk") or die "cannot open new tree file\n";

     while(<TRI>)
     {
      
	     # my $_ =~s///
	     
	     my $line = $_;
      my $line2;

      my $line3 = $line;
      $line3 =~ s/\)\b\d+/\)/g;
      $line3 =~ s/:\b\d\.\d+//g;

      foreach my $key(sort(keys(%seqid)))
      {
      $line =~s/$key:/$seqid{$key}:/;
      $line2 = $line;
     $line2 =~ s/\)\b\d+/\)/g;
      $line2 =~ s/:\b\d\.\d+//g;
      }
      print TRO $line;
     }
