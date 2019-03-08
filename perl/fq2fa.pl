#!/home/SoftWare/perl-5.26/bin/ -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fq,$fa);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fq:s"=>\$fq,
	"fa:s"=>\$fa,
			) or &USAGE;
&USAGE unless ($fq and $fa );
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@gogo.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -fq -fa -h

Usage:
  Options:
  -fq	<file>	input fq file
  -fa	<file>	output fa file
  -h         Help

USAGE
        print $usage;
        exit;
}
if($fq=~/gz$/){
    open IN,"zcat $fq|";
    }else{
        open IN,"$fq";
        }
open OUT,"|bgzip  > $fa.gz";
my $n=1;
while(<IN>){
    $_=~s/[\n\r]//g;
    if ($n%4==1) {
        my $name=$_;
        $name=~s/^@//;
        print OUT ">$name\n";
        }
        elsif ($n%4==2) {
                my $seq=$_;
                print OUT "$seq\n";
       }
   $n++;
}
close IN;
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

