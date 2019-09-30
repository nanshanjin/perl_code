#!/home/SoftWare/perl-5.26/bin/ -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$fout,$ft);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fin,
    "t:s"=>\$ft,
	"o:s"=>\$fout,
			) or &USAGE;
&USAGE unless ($fin );
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@sinotechgenomics.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -i -o -h

Usage:
  Options:
  -i	<file>	input file name
  -t
  -o	<file>	
  -h         Help

USAGE
        print $usage;
        exit;
}
open IN,"$fin";
open OUT,">$fout";
my @counts;

while(<IN>){
    $_=~s/[\n\r]//g;
    if(/^circRNA_id/){#
        my @array=split/\t/,$_;
        print OUT join("\t",@array),"\n";
        for(my $i=0;$i<=$#array;$i++){
            if($array[$i]=~/$ft/){
                push(@counts,$i);
                }
            }
        }else{
            my @array=split/\t/,$_;
            my $num=0;
            for(my $i=0;$i<=$#counts;$i++){
                if($array[$counts[$i]]>10){#最大值大于10
                    $num++;
                    }
                }
            if($num >= 9){#至少在9组中大于10
                my $max=$array[$counts[0]];
                foreach(@array[$counts[0]..$#counts]){
                    $max =$_ if $_ > $max;
                }
                if($max > 30){
                    print OUT join("\t",@array),"\n";
                    }
                }
        }
    }
close IN;
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

