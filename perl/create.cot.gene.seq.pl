use strict;
use Bio::SeqIO;
#### near to be finished
#    but it seems there is something of genes output by the sb-group
#
#
my $gffinfo = $ARGV[0];
open(GFF, $gffinfo) or die "cannot open the gff file $gffinfo due to $!.\n";        

my $outfile = $ARGV[2];
open(OUT,">",$outfile) or die "can not open out file $outfile due to $!.\n";

my ($thisgene, $lastgene) = ("", "");
my ($thischro, $lastchro) = ("", "");
my ($this_orient, $last_orient) = ("", "");
my ($thisid, $lastid)   = ("", "");
my ($ids, $function);
my ($this_gene_start, $this_gene_end);
my ($last_gene_start, $last_gene_end);

my ($last_exon_start, $last_exon_end) = (-1, -1);
my ($this_exon_start, $this_exon_end)  = (-1, -1);
my $ncds = "";
my $exon_num;
my $chroseq;
my $chroseqlen;
my $chro;
my $exons;
my $introns;
my $intergenic;
my ($utr3, $utr5);

while(<GFF>)
{
   $_ =~s/[\n\r]//g;
   my @array = split(/\s+/, $_);
#print "@array[0..$#array]\n";   

   $thischro = $array[0]; 
   if($thischro !~/^gr/){next;} ### only the genes assembled into pseudomolecules will be analyzed

   if($thischro ne $lastchro)
   {
print "reading chroseq .... $thischro.. \n";
      $chro = $thischro;
      $chroseq = get_chro($chro);
      $chroseqlen = length($chroseq);
print $chro."###############chro\n";
   }

   if($array[2] eq "mRNA")  ## start a new gene
   {
print "in mRNA \n";
       #### update the last gene
       $last_orient=$this_orient;
       $lastgene = $thisgene;
       $last_gene_start = $this_gene_start;
       $last_gene_end   = $this_gene_end;
       #### new gene
       $this_gene_start = $array[3];
       $this_gene_end   = $array[4];
       $this_orient= $array[6];
       if($this_gene_start > $last_gene_end)
       {
          $intergenic = $lastgene.",".substr($chroseq, $last_gene_end+1, $this_gene_start-$last_gene_end);
       }
       if($lastgene ne "")
       {
          my @arr = split(/,/, $intergenic);
          my $len = length($arr[1]);
          my $gc = calculateGCinNoncodingRegion(uc($arr[1]));
          my $intergenicgc = $gc;
          my @arr = split(/,/, $exons);
          my $str = join("", @arr[1..$#arr]);
          my $len = length($str);
          my $gc = calculateGCinCodingRegion(uc($str));
	  print OUT">$lastgene $len $thischro\n".uc($str)."\n";
       $thisgene = $array[8];
       $exons = $thisgene;
       $introns = $thisgene;
       $utr5 = $thisgene;
       $utr3 = $thisgene;
    
       $lastchro = $thischro;
       $exon_num = 0;
   }
   elsif($array[2] eq "CDS") ## a new exon
   {
print "in CDS\n";
       $this_exon_start = $array[3];
       $this_exon_end   = $array[4];

       $exon_num ++;
       if($exon_num > 1)
       {
       ### store introns
           if($this_orient eq "+")
           {
               $introns .= ",".substr($chroseq, $last_exon_end, $this_exon_start-$last_exon_end);
#               print $last_exon_end." ".$this_exon_start."\n";
           }
           else
           {
               $introns .= ",".negative(substr($chroseq, $this_exon_end, $last_exon_start-$this_exon_end));
#               print $this_exon_end." ".$last_exon_start."\n";
           }
       }

       ### store exons
       my $longer;
       if($this_orient eq "+")
       {
print $chroseqlen." ".$this_exon_start." ".$this_exon_end."\n";
          $exons .= ",".substr($chroseq, $this_exon_start - 1, $this_exon_end-$this_exon_start +1);
          #$longer = substr($chroseq[$chro] -> seq(), $this_exon_start-99, 102);
       }
       else
       {
print $chroseqlen." ".$this_exon_start." ".$this_exon_end."\n";
          $exons .= ",".negative(substr($chroseq, $this_exon_start -1 , $this_exon_end-$this_exon_start +1));
          #$longer = negative(substr($chroseq[$chro] -> seq(), $this_exon_end-2, 100));
       }
       ##########################################
       $last_exon_start = $this_exon_start;
       $last_exon_end   = $this_exon_end;
   }
}
close($gffinfo);
  
sub negative()
{
   my $seq = $_[0];
   $seq = reverse($seq);
   $seq = uc($seq); ### this is quite importand for some genome sequence may compose both uppercase and lowercase
   $seq =~ s/A/B/g;
   $seq =~ s/T/A/g;
   $seq =~ s/B/T/g;
   $seq =~ s/C/B/g;
   $seq =~ s/G/C/g;
   $seq =~ s/B/G/g;
   return $seq;
}

sub calculateGCinCodingRegion()
{
      my @seq =split(//, $_[0]);

      my ($gc1, $gc2, $gc3) = (0, 0, 0);
      if($#seq >= 0)
      {
      for(my $i = 0; $i <=$#seq-3; $i+=3)
      {
          if($seq[$i] eq "G" || $seq[$i] eq "C")
          {$gc1 ++;}
          if($seq[$i+1] eq "G" || $seq[$i+1] eq "C")
          {$gc2 ++;}
          if($seq[$i+2] eq "G" || $seq[$i+2] eq "C")
          {$gc3 ++;}
      }      
      $gc1 = 3*$gc1/($#seq + 1);
      $gc2 = 3*$gc2/($#seq + 1);
      $gc3 = 3*$gc3/($#seq + 1);
      }
      my $tmpstr = $gc1." ".$gc2." ".$gc3;

      return $tmpstr;
}

sub calculateGCinNoncodingRegion() ### here for intergenic region
{
     my @seq = split(//, $_[0]);
     my $GC = 0;
     for(my $i=0; $i<=$#seq; $i++)
     {
        if($seq[$i] eq "G" || $seq[$i] eq "C" || $seq[$i] eq "g" || $seq[$i] eq "c")
        {$GC++;}
     }
     if($#seq eq 0){return 0;}
     $GC /= $#seq;
     if($seq[0] =~ /^\w/)
     {return $GC;}
     else
     {return -1;}
}

sub get_chro()
{
   my $chro = $_[0];
   my $stream = Bio::SeqIO -> new(-file => $ARGV[1], -format => 'fasta');
   while(my $chroseq = $stream -> next_seq())
   {
      my $id = $chroseq -> display_id();
      if($id eq $chro){
      return $chroseq -> seq();
      }
   }
}

