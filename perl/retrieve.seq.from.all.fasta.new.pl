#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($id,$infa,$oufa);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"id:s"=>\$id,
    "infa:s"=>\$infa,
	"oufa:s"=>\$oufa,
			) or &USAGE;
&USAGE unless ($id and $infa and $oufa);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@sinotechgenomics.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -id -infa -oufa

Usage:
  Options:
  -id	<file>	input id
  -infa   <file>  input fa
  -oufa	<file>	out fa
  -h         Help

USAGE
        print $usage;
        exit;
}

my $idfile=$id;

open(ID, $idfile) or die "cannot open $idfile due to $!\n";

my %hash;
while(<ID>)
{
   $_ =~ s/[\n\r]//g;
   $hash{$_}=1;
}

my $infastafile = $infa;
my $is = Bio::SeqIO -> new(-format=>'fasta', -file => $infastafile);

my $outfastafile = $oufa;
my $os = Bio::SeqIO -> new(-format=>'fasta', -file => ">".$outfastafile);

while(my $seqobj = $is -> next_seq())
{
   my $id = $seqobj -> display_id();

   #print $id."\n";
   if(exists $hash{$id}){
       $seqobj -> display_id($id);
       $os -> write_seq($seqobj);
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

