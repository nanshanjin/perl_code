use strict;

my $blastfile = $ARGV[0];
open(IN,$blastfile) or die "can not open blast file $blastfile due to $!.\n";
my $outfile = $ARGV[1];
open(OUT,">",$outfile) or die "can not open out file $outfile due to $!.\n";

my $num = "";
my $num1 = "";
my $order=1;

while(<IN>)
{
	$_ =~s/[\n\r]//g;
	my @array = split(/\t/,$_);
	$num=$array[0];

   if ($num !~$num1)
   {
	$order=1;   
	print OUT $array[0]."\t".$array[1]."\t".$order."\n";
	
   }
   else
   {
	$order++;
	print OUT $array[0]."\t".$array[1]."\t".$order."\n";
	}
	
   	   
   	$num1 =$num;  
   	
   }
   	
   
   	


