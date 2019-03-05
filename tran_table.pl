#!/home/SoftWare/perl-5.26/bin/ -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$fout);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fin,
	"o:s"=>\$fout,
			) or &USAGE;
&USAGE unless ($fin );
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang;
Script:			$Script
Description:
    script to tran table data	
	eg:
	perl $Script -i <infile.table> -o <outfile.table> -h
Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	
  -h         Help

USAGE
        print $usage;
        exit;
}
my $i=0;
my @data;
open IN,"$fin";
while(<IN>){
    $_=~s/[\n\r]//g;
    my @line=split/\t/,$_;
    #my @line=split/\t/ ? split/\t/ : split/\s/;
    for my $j(0..$#line){
        $data[$j][$i]=$line[$j];
        }
    $i++;
    }
close IN;
open OUT,">$fout";
for(@data){
    print OUT join("\t",@{$_}),"\n";
    }
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

