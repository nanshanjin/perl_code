use strict;
my @arr=("D","P","C","O","F","G","S");
my %type=("D"=>"1","P"=>"1","C"=>"1","O"=>"1","F"=>"1","G"=>"1","S"=>"1");
#my %hash;
my @file=glob("*Bacteria.xls");
for my $file (@file){
    #print $file,"\n";
    my $filename=(split/\./,$file)[0];
    open IN,$file;
    open OUT,">$filename\_class.xls";
    my ($d,$p,$c,$o,$f,$g,$s);
    my %hash;
    my $start=""; ##用于记录上一次的级别
    while(<IN>){
        $_=~s/[\n\r]//g;
        chomp;
        my @array=split/\t/,$_;
        $array[5]=~s/\s+(.*)/$1/g;
        if(!exists $type{$array[3]}){
            next;
        } else {
            my $end=0; ##用于记录是否达到或越过上一次的级别
            $hash{$array[3]}="$array[3]__$array[5]";
            for(my $i=0;$i<=$#arr;$i++){
                ## 判断是否到达上一次的级别
                if($start eq $arr[$i]){
                    $end = $i;
                }
                ## 如果已经达到当前行所在的级别
                if($array[3] eq $arr[$i]){
                    print OUT $hash{$arr[$i]},"|\n";
                    $end = 0;
                    last;
                } else {
                    ## 如果已经到达或超过上一次的级别
                    if( $end ne 0 and ($i - $end) ge 1 ){
                        $hash{$arr[$i]} = "";
                    }
                    print OUT $hash{$arr[$i]},"|";
                }
            }
            $start = $array[3];
        }
    }
    close IN;
}

