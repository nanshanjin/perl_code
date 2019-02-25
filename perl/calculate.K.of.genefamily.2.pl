use strict;
use Bio::SeqIO;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::Seq::EncodedSeq;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Clustalw;
use SeqStatistics;


my $input = $ARGV[0]; ### cds sequence file

my $is = Bio::SeqIO -> new(-file => $input, -format => 'fasta');

my $output = $input.".K.txt";
open(OUT, ">".$output) or die "cannot open $output due to $!.\n";

my %id2seqobj;


while(my $seqobj = $is -> next_seq())
{
   my @a = split(/\s/, $seqobj->display_id());
   my $id = $a[0];
   $id2seqobj{$id} = $seqobj;
}

my @id = sort(keys(%id2seqobj));

for(my $i=0; $i<$#id; $i++)
{
   for(my $j=$i+1; $j<=$#id; $j++)
   {
       my @kv = calculateK($id2seqobj{$id[$i]}, $id2seqobj{$id[$j]});
       print OUT "$id[$i] $id[$j] @kv[0..$#kv]\n";#$i $j 
   }
} 

close($output);

sub calculateK
{
   my $seqobj1 = $_[0];
   $seqobj1 -> seq(uc($seqobj1 -> seq));
   my $id1 = $seqobj1 -> display_id();

   my %dna_hash;
   $dna_hash{$seqobj1 -> display_id} = $seqobj1;

   my $seqobj2 = $_[1];
   $seqobj2 -> seq(uc($seqobj2 -> seq));
   my $id2 = $seqobj2 -> display_id();

   $dna_hash{$seqobj2 -> display_id} = $seqobj2;

   my @prot_arr;
   push @prot_arr, $seqobj1 ->translate();
   push @prot_arr, $seqobj2 ->translate();

   my $os_prot = Bio::SeqIO -> new(-file=> ">prot.fasta", -format=>"fasta");
   $os_prot -> write_seq($seqobj1 -> translate());
   $os_prot -> write_seq($seqobj2 -> translate());

   system("clustalw -infile=prot.fasta");

   my $is_prot_aln = Bio::AlignIO -> new(-file=>"prot.aln", -format=>"CLUSTALW");
   my $prot_aln = $is_prot_aln -> next_aln();
   system("rm prot.fasta prot.aln prot.dnd");

   my $dna_aln = &aa_to_dna_aln($prot_aln, \%dna_hash);
      
   my $stats = new SeqStatistics;
   my $result = $stats->calc_all_KaKs_pairs($dna_aln);

   my ($Da, $Ds, $Dn, $N, $S, $S_d, $N_d);
   for my $an (@$result)
   {
     #print "comparing ". $an->{'Seq1'}." and ". $an->{'Seq2'}. " \n";
     for (sort keys %$an )
     {
          next if /Seq/;
          #printf("%-9s %.4f \n",$_ , $an->{$_});
          if($_ eq "D_n"){$Dn = $an->{$_}};
          if($_ eq "D_s"){$Ds = $an->{$_}};
          if($_ eq "S_d"){$S_d = $an->{$_};}
          if($_ eq "N_d"){$N_d = $an->{$_};}
          if($_ eq "S"){$S = $an->{$_};}
          if($_ eq "N"){$N = $an->{$_};}

     }
     #print "\n\n";
   }
  
#   my $length1= length($seqobj1 -> seq);
#   my $length2= length($seqobj2 -> seq);

   if($Dn !~ /\d/){$Dn = -2;}
   if($Ds !~ /\d/){$Ds = -2;}

   return ($Dn, $Ds);
}
