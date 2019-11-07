use strict;
use Bio::SeqIO;

my %hash;my %hash2;


my $infastafile = $ARGV[0];
my $is = Bio::SeqIO -> new(-format=>'fasta', -file => $infastafile);
open GFF,"<$ARGV[1]";
open ID,"<$ARGV[2]";
open OUT,">$ARGV[3]";


while(my $seqobj = $is -> next_seq())
{
   my $id = $seqobj -> display_id();
   my $seq = $seqobj->seq;
   $seq=~s/\s+//g; $seq=~s/\r//g;$seq=~s/\n//g;
   $hash{$id}=$seq;
}

$/="\n";
	while (<GFF>) 
	{
		chomp; 
		next if ($_=~m/^$/);
		#my ($chr,$type,$start,$end,$strand,$gene)=(split(/\s+/,$_))[0,2,3,4,6,8];
		my @a=split;
		if ($a[1] eq "mRNA") 
		{
			#$gene=(split(/\;/,$gene))[0];
			#$gene=(split(/\=/,$gene))[1];
			$hash2{$a[5]}="$a[0]\t$a[2]\t$a[3]\t$a[4]\n"
		}
	}
close GFF;

$/="\n";
	while (<ID>) 
	{
		chomp; 
		next if ($_=~m/^$/);
		$_=~s/\s+//g; $_=~s/\r//g;
		my ($chr2,$start2,$end2,$strand2)=(split(/\s+/,$hash2{$_}))[0,1,2,3];
		my $gene_len=$end2-$start2+1;
		
		if ($strand2 eq "-") 
		{	
			my $tmp=substr($hash{$chr2},$end2,1000);
			$tmp=~tr/ATCGN/TAGCN/;
			$tmp=reverse $tmp;
			print OUT ">$_\n$tmp\n";
		}
		elsif ($strand2 eq "+") 
		{	
			my $tmp=substr($hash{$chr2},$start2-1001,1000);
			print OUT ">$_\n$tmp\n";
		}

	}
close ID;