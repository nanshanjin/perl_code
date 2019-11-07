#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw/sum/;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($length,$vcf,$out,$vpos,$windows,$sample,$type,$eff,$diff,$region,$depth);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
    "sample:s"=>\$sample,
    "diff:s"=>\$diff,
    "region:s"=>\$region,
    "eff:s"=>\$eff,
    "vpos:s"=>\$vpos,
    "depth:s"=>\$depth,
    "type:s"=>\$type,
    "length:s"=>\$length,
    "windows:s"=>\$windows,
			) or &USAGE;
&USAGE unless ($vcf and $sample and $out);
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
  --sample <string>;must be given;option:sampe1,sample2
  --diff <string>;option:diff or equal;
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
$diff||="all";
$region||="all";
$eff||="all";
#my %hash_eff;
#if($eff ne "all"){
#    my @varient_t=split/\,/,$eff;
#    for(@varient_t){
#        $_=uc($_);
#        $hash_eff{$_}=1;
#        }
#}
$type||="all";
$vpos||="all";
#my %hash_vpos;
#if($vpos ne "all"){
#    my @varient_p=split/\,/,$vpos;
#    for(@varient_p){
#        $hash_vpos{$_}=1;
#        }
#    }
$depth||="all";
my %groups;
my @Indi;
my @smaples=split/\,/,$sample;
for(@smaples){
    $groups{$_}=1;
    }
if($vcf=~/gz/){
    open(IN,"zcat $vcf |") or die "Cannot open $vcf!\n";
}else{
    open(IN,$vcf) or die "Cannot open $vcf!\n";
    }
open OUT,">$out/variant.txt";
#open OUT,">$out";
my %hash;
while (<IN>) {
    $_=~s/[\n\r]//;
    next if ($_ eq "" || /^$/ ||/^##/);
    if (/^#/) {
        (undef,undef,undef,undef,undef,undef,undef,undef,undef,@Indi)=split(/\t/,$_);
        my @sams;
        foreach my $indi (@Indi) {
            next if !exists $groups{$indi};
            push (@sams,$indi,"$indi\_depth");
            }
        print OUT join("\t","#Pos","Chr","Reference","Type",@sams,"Annotation"),"\n";
        #print OUT "#Pos\tChr\tReference\tType\tSample1 Genotype\tdepth\tSample2 Genotype\tSample2_depth\tAnnotation\n";
    }else{
        my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\t/,$_);
        my @formats=split/\:/,$format;
        my $dp_position;
        for(my $i=0;$i<=$#formats;$i++){
            if($formats[$i] eq "AD"){
                $dp_position=$i;
                }
            }
        #print $dp_position,"\n";
        #my (undef,$variant_pos,$variant_type,@others)=split/\|/,$info;
        my $check_type=0;
        if($eff ne "all"){
            my @variant_type=split/\,/,$eff;
                for(my $i=0;$i<=$#variant_type;$i++){
                    $variant_type[$i]=uc($variant_type[$i]);
                    #print $variant_type[$i],"\n";
                    if($info=~/$variant_type[$i]/){
                        $check_type=1;
                        next;
                        }#else{
                         #   $check_type=1;
                         #   }
                    #print $t,"\n";
                }
        }
        next if $check_type==0 and $eff ne "all";
        my $check_pos=0;
        if($vpos ne "all"){
            my @variant_pos=split/\,/,$vpos;
            for my $p (@variant_pos){
                if($info=~/$p/){
                    $check_pos=1;
                    next;
                    }#else{
                     #   $check_pos=1;
                     #   }
                }
            }
        next if $check_pos==0 and $vpos ne "all";
        ###variant chr start end
        if($region ne "all"){
            my ($region_chr,$region_start,$region_end)=split/\,/,$region;
            next if $chr ne $region_chr;
            next if $pos le $region_start or $pos ge $region_end;
            }
        ########################
        my @ale=split(/\,/,join(",",$ref,$alt));
        my @getype;
        for (my $i=0;$i<@indi;$i++) {
            if(exists $groups{$Indi[$i]}){
                my $geno=(split(/\:/,$indi[$i]))[0];
                my $dep=(split/\:/,$indi[$i])[$dp_position];
                my @dep_AD=split/\,/,$dep;
                my $depall=sum(@dep_AD);
                my ($g1,$g2)=split(/\//,$geno);
                if ($geno eq "./.") {
                    push(@getype,"N/N",$depall);
                }else{
                    #push(@getype,"$ale[$g1]/$ale[$g2]",$dep);
                    push(@getype,join("/",sort($ale[$g1],$ale[$g2])),$depall); 
                    }
        }
    }
    my ($sam1,$dep1,$sam2,$dep2)=(split/\t/,join("\t",@getype));
    next if $sam1 eq "N/N" or $sam2 eq "N/N";
    #print "$sam1\t$dep1\t$sam2\t$dep2\n";
    #####sample diff option
    if($diff ne "all"){
        if($diff eq "diff"){
            next if $sam1 eq $sam2;
            }else{
                next if $sam1 ne $sam2; 
                }
        }
    ##########
    ##########
    ##########
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
            }
    }
    ##########
    ###depth filter
    if($depth ne "all"){
        my ($s1,$e1,$s2,$e2)=split/\,/,$depth;
        if(($dep1 > $s1 and $dep1 < $e1) and ($dep2 > $s2 and $dep2 < $e2)){
            $hash{"$chr\t$pos\t$ref\t$alt"}=join("\t",$sam1,$dep1,$sam2,$dep2)."\t".join(";",@annotation);
            #print OUT "$chr\t$pos\t$ref\t$alt\t",join("\t",$sam1,$dep1,$sam2,$dep2),"\t","$variant_pos\t$variant_type","\n";
            }
    }else{
        $hash{"$chr\t$pos\t$ref\t$alt"}=join("\t",$sam1,$dep1,$sam2,$dep2)."\t".join(";",@annotation);
        #print OUT "$chr\t$pos\t$ref\t$alt\t",join("\t",$sam1,$dep1,$sam2,$dep2),"\t","$variant_pos\t$variant_type","\n";
        }
        #print OUT "$chr\t$pos\t$ref\t$alt\t",join("\t",$sam1,$dep1,$sam2,$dep2),"\t","$variant_pos\t$variant_type","\n";
    }
}
foreach my $keys (sort keys %hash){
    my ($chr,$pos,$ref,$alt)=split/\t/,$keys;
    if($type eq "SNP"){
        if((length($ref) eq "1" and length($alt) eq "1") or (length($ref) eq "1" and $alt=~/\,/ and length($alt) eq "3")){
            print OUT "$chr\t$pos\t$ref\t$alt\t$hash{$keys}\n";
            }
        }elsif($type eq "INDEL"){
            if(length($ref) > "1" or ($alt !~/,/ and length($alt) > "1") or ($alt=~/,/ and length($alt) > "3")){
                my $lens=length($ref)-length($alt);
                if($length){
                    my ($len1,$len2)=split/\,/,$length;
                    if($len1 < $lens and $lens < $len2){
                        print OUT "$chr\t$pos\t$ref\t$alt\t$hash{$keys}\n";
                        }
                    }else{
                        print OUT "$chr\t$pos\t$ref\t$alt\t$hash{$keys}\n";
                        }
                }
            }else{
                print OUT "$chr\t$pos\t$ref\t$alt\t$hash{$keys}\n";
               }
    #print OUT "$keys\t$hash{$keys}\n";
    }
close IN;
close OUT;

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

