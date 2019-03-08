#!/home/SoftWare/perl-5.26/bin/ -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$fout,$g);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fin,
    "g:s"=>\$g,
	"o:s"=>\$fout,
			) or &USAGE;
&USAGE unless ($fin );
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@gogo.com;
Script:			$Script
Description:
    rename the header line
	eg:
	perl $Script -i -o -h

Usage:
  Options:
  -g	<file>
  -i
  -o	<file>	
  -h         Help

USAGE
        print $usage;
        exit;
}
open IN,"$g";
my %hash;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($sample,$group)=split/\t/,$_;
    $hash{$sample}=$group;
    }
close IN;
open IN,"$fin";
open OUT,">$fout";
my $line=0;
while(<IN>){
    $_=~s/[\n\r]//g;
    if($line==0){
        my @array=split/\t/,$_;
        my @new;
        for(my $i=0;$i<=$#array;$i++){
            if(exists $hash{$array[$i]}){
                push @new,$hash{$array[$i]};
                }else{
                    push @new,$array[$i];
                    }
            }
            print OUT join("\t",@new),"\n";
        }else{
            print OUT "$_\n";
            }
    $line++;
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

