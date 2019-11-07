use strict;
open IN,"<old.gff";
open OUT,">new.num.gff";
open OUT2,">lens";

my %hash;
my $chra="";
my $num=1;
while(<IN>){
	$_=~s/[\n\r]//g;
	my $chr=(split/\t/,$_)[0];
	$hash{$chr}+=1;
	if($chr eq $chra){		
		print OUT"$_\t$num\n";
		$num++;
	}else{
		$num=1;
		print OUT"$_\t$num\n";
		$chra=$chr;
		$num++;
	}
}
foreach my $key (sort  keys %hash){
	print OUT2"$key\t$hash{$key}\n";
}


