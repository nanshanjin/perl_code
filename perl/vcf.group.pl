#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw/sum/;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($length,$vcf,$out,$vpos,$windows,$group,$type,$eff,$region,$depth);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
    "group:s"=>\$group,
    "region:s"=>\$region,
    "eff:s"=>\$eff,
    "vpos:s"=>\$vpos,
    "depth:s"=>\$depth,
    "type:s"=>\$type,
    "length:s"=>\$length,
    "windows:s"=>\$windows,
			) or &USAGE;
&USAGE unless ($vcf and $out);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -vcf -out -h

Usage:
  Options:
  --vcf	<file>	input vcf file;must be gz
  --group <string>;must be given;option:pop1,sample1,sample2;pop2,sample3,sample4
  --region <string>the genome region option;chr,start,end
  --eff <string> varient eff;High,Moderate,Low,Modifier
  --vpos <string> varient pos;upstream_gene_variant,3_prime_UTR_variant,.....
  --depth <number> sample depth;start1,end1,start2,end2
  --type <string> variant type;option:SNP INDEL or all
  --length <number>    indel length;
  --windows <number>
  --out	<dir>	out dir;
  --h         Help

USAGE
        print $usage;
        exit;
}

mkdir $out if(!-d $out);
$region||="all";
$eff||="all";
$type||="0";
$vpos||="all";
$depth||="all";
$length||="0,10000000";
my %groups;
my %group_list;
my @Indi;
my @pop=split/\;/,$group;
my @gro;
for(@pop){
    my ($pop_group,@samples)=split/\,/,$_;
    push(@gro,$pop_group);
    for(@samples){
        $groups{$_}=1;
        $group_list{$_}=$pop_group;
        }
    }

if($vcf=~/gz$/){
    open(IN,"zcat $vcf |") or die "Cannot open $vcf!\n";
}else{
    open(IN,$vcf) or die "Cannot open $vcf!\n";
    }
open OUT,">$out/variant.txt";
my %hash;

while (<IN>) {
    $_=~s/[\n\r]//;
    next if ($_ eq "" || /^$/ ||/^##/);
    if (/^#/) {
        (undef,undef,undef,undef,undef,undef,undef,undef,undef,@Indi)=split(/\t/,$_);        
        print  join("\t","#Pos","Chr","Reference","Type",@gro,"Annotation"),"\n";
    }else{
        my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\t/,$_);
        next if ($type eq "INDEL" and ((length($ref) == "1" or ($alt !~/,/ and length($alt) == "1") or ($alt=~/,/ and length($alt) == "3"))));
        next if ($type eq "SNP" and (length($ref) > "1" or ($alt !~/,/ and length($alt) > "1") or ($alt=~/,/ and length($alt) > "3")));
        my $lens=length($ref)-length($alt);
        my ($len1,$len2)=split/\,/,$length;
        print join("\t",$ref,$alt),"\n";
        next if ($type eq "INDEL" and ($len1 > $lens or $lens > $len2));
        my @formats=split(/:/,$format);
        for (my $i=0;$i<@indi;$i++) {
            next if !exists $groups{$Indi[$i]};
            for(my $j=0;$j<@formats;$j++){{
                if($format[$j] eq "GT"){
                    my ($g1,$g2)=split(/\:/,$indi[$j]);

                    }
                if($format[$j] eq "AD"){
                    
                    }
                }
            my $geno=(split(/\:/,$indi[$i]))[0];
            }
    ###anno
    my @annotation;
    if($info=~/ANN=([^\;]*)/g){
        my @ann=split(/\,/,$1);
        for (my $i=0;$i<@ann;$i++) {
            my @str=split(/\|/,$ann[$i]);
            $str[0]||="--";
            $str[1]||="--";
            $str[2]||="--";
            $str[3]||="--";
            $str[4]||="--";
            my $ann=join("|",$str[0],$str[1],$str[2],$str[3],$str[4]);
            push @annotation,$ann;
            $hash{"$chr\t$pos\t$ref\t$alt"}{anno}=join(";",@annotation);
            }
        }

    }
}

open IN,"$out/variant.txt";
open OUT,">$out/depth.txt";
my %courage;

while(<IN>){
     $_=~s/[\n\r]//g;
     next if /^#/;
     my ($chr,$pos,undef,undef,undef,$dep1,undef,$dep2,@others)=split/\t/,$_;
     my $win_num=int($pos/$windows)+1;
     my $deps=($dep1+$dep2)/2;
     $courage{$chr}{$win_num}{dep}+=$deps;
     $courage{$chr}{$win_num}{num}++;
}
close IN;
foreach my $keys (sort keys %courage){
    foreach my $num (sort {$a<=>$b} keys %{$courage{$keys}}){
        my $start=1+$windows*($num-1);
         my $end=$windows*$num;
          my $deps=$courage{$keys}{$num}{dep}/$courage{$keys}{$num}{num};
          print OUT "$keys\t$num\t$deps\n";
    }
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub deps{
    my ($dep)=@_;
    my @return;
    foreach my $pop (sort keys %$dep){
        my @value=values %$dep;
        my $d=sum(\@value)/(scalar @value);
        push(@return,"$pop:$d");
    }
    return @return;
}
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

