#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw(max);
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$fout,$windows);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fin,
	"o:s"=>\$fout,
	"windows:s"\$windows,
			) or &USAGE;
&USAGE unless ($fin,$fout,$windows);

open IN,$fin;
open OUT,">",$fout;

my %hash;
my $win=$windows;

while(<IN>){
    $_=~s/[\n\r]//g;
    my @array=split;
    next if /^#/;
    my $win_num=int($array[1]/$win)+1;
    #print "$win_num\n";
    $hash{$array[0]}{$win_num}++;
}
my $col=1;
foreach my $keys (sort keys %hash){
   # print "$keys\n";
    foreach my $num (sort keys %{$hash{$keys}}){
        
        my $start=1+$win*($num-1);
        my $end=$win*$num;
        my $snp=$hash{$keys}{$num}/$win;
       #print OUT "$keys\t$num\t$hash{$keys}{$num}\n";
        
        print OUT "$keys\t$start\t$end\t$snp\tfill_color=col$col\n";
        
        }
        $col++;
    }
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang;
Script:			$Script
Description:
	this program is to slide the windows of snp or indel
	eg:
	perl $Script -i -o -windows -c

Usage:
  Options:
  -i	<file>	input vcf  file
  -o	<file>	output windows to draw circos
  -windows windows numbers
  -h         Help

USAGE
        print $usage;
        exit;
}
