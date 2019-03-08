#!/home/SoftWare/perl-5.26/bin/ -w
#use strict;
#use warnings;
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
Contact:        nanshan.yang\@gogo.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -i -o -h

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	
  -h         Help

USAGE
        print $usage;
        exit;
}
open IN,$fin;
open OUT,">$fout";
my %hash;
my @snps;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($sample,$snp,$type)=split/\t/,$_;
    push(@snps,$snp);
    $hash{$sample}{$snp}=$type;
    }
close IN;
my %sort;
my @snp_new = grep { ++$sort{$_} < 2 } @snps;
print OUT "sampleid";
for my $snp (sort @snp_new) {
    print OUT "\t$snp";
}
print OUT "\n";
for my $sample(sort keys %hash){
    print OUT "$sample";
    for my $snp (sort @snp_new){
        if(exists $hash{$sample}{$snp}){
            print OUT "\t$hash{$sample}{$snp}";
            }else{
                print OUT "\tNA";
                }
    }
    print OUT "\n";
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

