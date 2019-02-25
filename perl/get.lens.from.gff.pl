#!/usr/bin/perl -w
##get lens file for dot_plot
use strict;
my $ingff=$ARGV[0];
open IN,$ingff or die "$!";
my $outlens=$ARGV[1];
open OUT,">",$outlens or die "$!";

my $i=1;
my %hash;

while(<IN>){
        $_=~s/[\n\r]//g;
	my @array=split;
        $array[0]=~s/Chr//g;
	
	if(defined $hash{$array[0]}){
		$hash{$array[0]}++;
	}else{
		$hash{$array[0]}=1;
	}
	$i++;
}
close(IN);
foreach my $element (sort {$a <=> $b} keys %hash){
	print OUT $element."\t".$hash{$element}."\n";
}
#哈希排序说明：
#按哈希键的数值大小排序：
#my @key =sort {$a <=> $b} keys %hash;
#@key里头存的是按哈希键的数值大小排序后的键。
#按哈希值的数值大小排序：
#my @key =sort {$hash{$a} <=> $hash{$b}} keys %hash;
#@key里头存的是按哈希值的数值大小排序后的键

      
