use strict;
open IN,"<gh_nn_gh_nn.block.rr.txt";
open OUT,">gh_nnD_gh_nnD.block.rr.txt";


my @line=();
while(<IN>){
	$_=~s/[\n\r]//g;
	@line[$#line+1]=$_;
}

for(my $i=0;$i<=$#line;$i++){
	if($line[$i] =~/MAX/)
	{
		print OUT$line[$i]."\n";
	}
	if($line[$i]=~/^the/ && $line[$i+1]=~/AD/){
		my @a=split/\s+/,$line[$i];
		print OUT"$line[$i]\n";
		for(my $j=1;$j<=$a[4];$j++){
			my @b=split/\s+/,$line[$i+$j];
			$b[0]=~/AD(\d+)g(\d+)/;
			#	print "$b[0]\n";
                        my $num1=$1;
			$b[2]=~/AD(\d+)g(\d+)/;
			my $num2=$1;
			if($num1 > 13 && $num2 > 13){
				print OUT"$line[$i+$j]\n";
			}
		}
		print OUT"$line[$i+1+$a[4]]\n";
	}
}


