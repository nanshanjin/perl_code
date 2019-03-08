#!/home/SoftWare/perl-5.26/bin/ -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fin,
			) or &USAGE;
&USAGE unless ($fin );
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@gogo.com;
Script:			$Script
Description:
	this program is to check md5
	eg:
	perl $Script -i -h

Usage:
  Options:
  -i	<dir>	input dir
  -h         Help

USAGE
        print $usage;
        exit;
}
my @fq=glob("$fin/*.gz");
for(@fq){
    my $fq_file=basename($_);
    my $out=`md5sum $_`;
    my ($out_md5,undef)=split/\s+/,$out;
    my $md5_file=$_.".md5";
    my $md5s="";
    open IN,"$md5_file";
    while(<IN>){
        $_=~s/[\n\r]//g;
        my ($md5,$file)=split;
        $md5s=$md5;
        }
    close IN;
    if($out_md5 ne $md5s){
        print "there is some problem with $fq_file,you need check!!!\n";
        }else{
            print "the $fq_file is ok,you can analysis\n";
            }
    #print "$fq_file\t$md5s\n";
    }
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

