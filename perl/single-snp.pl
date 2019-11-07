#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
use Data::Dumper;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$group);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"g:s"=>\$group,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $group);
my %group;
my %rgroup;
$rgroup{Total}=1;
if ($group){
	open In,$group;
	while (<In>){
		chomp;
		next if ($_ eq ""||/^$/ || /^#/);
		my ($id,$gr,undef)=split(/\s+/,$_);
		$group{$id}=$gr;
		$rgroup{$gr}=1;
	}
	close In;
}
open In,$fIn;
open MAF,">$fOut";
print MAF "#chr\tpos\tref\talt\t","\tformat\t",join("\t",sort keys %rgroup),"\n";
my @Indi;
while (<In>){
	chomp;
	next if ($_ eq ""||/^$/ ||/^##/);
	if (/^#/){
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\t/,$_);
		push @Indi,@indi;
	}else{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@Geno)=split(/\t/,$_);
		my %Geno;
		my %Alle;
		my @format=split(/\:/,$format);
		my @All=split(",",join(",",$ref,$alt));
		for (my $i=0;$i<@Geno;$i++){
			my @geno=split(/\:/,$Geno[$i]);
			my $gid="Total";
			if (exists $group{$Indi[$i]}){
				$gid=$group{$Indi[$i]};
			}
			for (my $j=0;$j<@geno;$j++){
				if ($format[$j] eq "GT"){
					next if ($geno[$j] eq "./.");
					my @alle=split(/\//,$geno[$j]);
					$Alle{Total}{$alle[0]}++;
					$Alle{Total}{$alle[1]}++;
					$Alle{$gid}{$alle[0]}++;
					$Alle{$gid}{$alle[1]}++;
					$Geno{Total}{$geno[$j]}++;
					$Geno{$gid}{$geno[$j]}++;
				}	
			}
		}
		my $maf=maf_calc(\%Alle);
		my $hw=hw_calc(\%Alle,\%Geno);
		my $shanon=shanon_calc(\%Alle);
		my $pic=pic_calc(\%Alle);
		my $home=home_calc(\%Alle);
		my @maf=split(/\s+/,$maf);
		my @hw=split(/\s+/,$hw);
		my @shanon=split(/\s+/,$shanon);
		my @pic=split(/\s+/,$pic);
		my @home=split(/\s+/,$home);
		my @out;
		my $n=0;
		foreach my $grid (sort keys %rgroup) {
			my @alle=sort{$Alle{$grid}{$b} <=> $Alle{$grid}{$a}} keys %{$Alle{$grid}};
			my @value=sort{$b<=>$a} values %{$Alle{$grid}};
			my @geno=sort{$Geno{$grid}{$b} <=> $Geno{$grid}{$a}} keys %{$Geno{$grid}};
			my @gvalue=sort{$b<=>$a} values %{$Geno{$grid}};
			push @{$out[$n]},join(":",join(",",@alle),join(",",@value),join(",",@geno),join(",",@gvalue));
			$n++;
		}
		for (my $i=0;$i<@maf;$i++) {
			push @{$out[$i]},$maf[$i];
			push @{$out[$i]},$hw[$i];
			push @{$out[$i]},$shanon[$i];
			push @{$out[$i]},$pic[$i];
			push @{$out[$i]},$home[$i];
		}
		my @Out;
		for (my $i=0;$i<@out;$i++) {
			push @Out,join("\:",@{$out[$i]});
		}

		print MAF join("\t",$chr,$pos,$ref,$alt,"alle:count:geno:count:maf:fixtion:hw:shanon:pic:hete:homo",@Out),"\n";
		
	}
}
close In;
close MAF;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub home_calc{
	my ($Geno)=@_;
	my @return;
	foreach my $pop(sort keys %$Geno){
		my @alle=sort{$$Geno{$pop}{$b} <=> $$Geno{$pop}{$a}} keys %{$$Geno{$pop}};
		my @value=sort{$b<=>$a} values %{$$Geno{$pop}};
		my $J=0;
		my $H=1;
		for (my $i=0;$i<@value;$i++){
			my $p=$value[$i]/sum(\@value);
			$J+=$p*$p;
		}
		$H=1-$J;
		push @return,"$H:$J";
	}
	return join("\t",@return);
}
sub pic_calc{
	my ($Geno)=@_;
	my @return;
	foreach my $pop(sort keys %$Geno){
		my @alle=sort{$$Geno{$pop}{$b} <=> $$Geno{$pop}{$a}} keys %{$$Geno{$pop}};
		my @value=sort{$b<=>$a} values %{$$Geno{$pop}};
		my $PIC=1;
		for(my $i=0;$i<@value;$i++){
			for(my $j=$i+1;$j<@value;$j++){
				my $p=$value[$i]/sum(\@value);
				my $q=$value[$j]/sum(\@value);
#				print $p*$p*$q*$q*2m
				$PIC-=$p*$p*$q*$q*2;
			}
			my $p=$value[$i]/sum(\@value);
			$PIC-=$p*$p;
		}
	#	print $PIC;die;
		push @return,$PIC;
	}
	return join("\t",@return);
}
sub shanon_calc{
	my ($Geno)=@_;
	my @return;
	foreach my $pop(sort keys %$Geno){
		 my @alle=sort{$$Geno{$pop}{$b} <=> $$Geno{$pop}{$a}} keys %{$$Geno{$pop}};
		 my @value=sort{$b<=>$a} values %{$$Geno{$pop}}; 
		 my $H=0;
		 for (my $i=0;$i<@value;$i++){
			$H+=-1* $value[$i]/sum(\@value) * log($value[$i]/sum(\@value));
		 }
		 push @return,$H;
	}
	return join("\t",@return);
}
sub maf_calc{
	my ($Geno)=@_;
	my @maf;
	foreach my $pop(sort keys %$Geno){
		my @alle=sort{$$Geno{$pop}{$b} <=> $$Geno{$pop}{$a}} keys %{$$Geno{$pop}};
		my @value=sort{$b<=>$a} values %{$$Geno{$pop}}; 
		my $info;
		if(scalar @alle == 1){
			push @maf,0
		}elsif(scalar @alle == 0){
			push @maf,"-";
		}else{
			push @maf,$alle[1]/sum(\@value);
		}
	}
	return join("\t",@maf);
}
sub hw_calc{
	my ($Alle,$Geno)=@_;
	my @HW;
	my @return;
	foreach my $pop(sort keys %$Geno){
		my @alle=sort{$$Alle{$pop}{$b} <=> $$Alle{$pop}{$a}} keys %{$$Alle{$pop}};
		my @value=sort{$b<=>$a} values %{$$Alle{$pop}};
		my @Gvalue=values %{$$Geno{$pop}};
		my %HWstat;
		my $sumE=0;
		my $sumO=0;
		for (my $i=0;$i<@alle;$i++){#now only for two
			for (my $j=$i+1;$j<@alle;$j++){
				my $p=$$Alle{$pop}{$alle[$i]}/sum(\@value);
				my $q=$$Alle{$pop}{$alle[$j]}/sum(\@value);
				my $geno=join("/",sort($alle[$i],$alle[$j]));
				if (!exists $$Geno{$pop}{$geno}) {
					$$Geno{$pop}{$geno}=0.0001;
				} 
				$HWstat{ob}{$geno}=$$Geno{$pop}{$geno}/sum(\@Gvalue);
				$HWstat{ex}{$geno}=$p*$q*2;
				$sumE+=$p*$q*2;
				$sumO+=$$Geno{$pop}{$geno}/sum(\@Gvalue);
			}
			my $geno=join("/",$alle[$i],$alle[$i]);
			if (!exists $$Geno{$pop}{$geno}) {
				$$Geno{$pop}{$geno}=0.0001;
			} 
			$HWstat{ob}{$geno}=$$Geno{$pop}{$geno}/sum(\@Gvalue);
			my $p=$$Alle{$pop}{$alle[$i]}/sum(\@value);
			$HWstat{ex}{$geno}=$p*$p;
		}
		my @order=sort {$HWstat{ob}{$a}<=>$HWstat{ob}{$b}}keys %{$HWstat{ob}};
		my @od;my @ex;
		for(my $i=0;$i<@order;$i++){
			push @od,$HWstat{ob}{$order[$i]};
			push @ex,$HWstat{ex}{$order[$i]};
		}#
		my $od=join(":",@od);
		my $ex=join(":",@ex);
		my $seg=Segregation($od,$ex,1);
		my $F;
		if ($sumE == 0) {
			$F=1;
		}else{
			$F=($sumE-$sumO)/$sumE;
		}
		push @return,join(":",$F,$seg);
	}
	return join("\t",@return);

}
sub sum{
	my ($a)=@_;
	my $sum=0;
	foreach my $n(@$a){
		$sum+=$n;
	}
	return $sum;
}
sub Segregation {#
		my ($theoretical_segregation,$segregation,$all)=@_;
		my @a=split ":",$theoretical_segregation;
		my @b=split ":",$segregation;
		if (scalar @a == 1 && scalar @b == 1) {
			return 0;
		}
		return "0.01" if (scalar @a != scalar @b || $all == 0) ;
		my @theoretical;
		my $a_sum=0;
		$a_sum+=$_ foreach (@a);
		push @theoretical,$_/$a_sum*$all foreach (@a);
		my $df=scalar @a -1;
		my $X2=0;
		if ($df == 1) {
			for (my $i=0;$i<@a ;$i++) {
				$X2+=X2df2($b[$i],$theoretical[$i]);
			}
		}else{
			for (my $i=0;$i<@a ;$i++) {
				$X2+=X2df1($b[$i],$theoretical[$i]);
			}
		}
		my $p=0;
		$p=Statistics::Distributions::chisqrprob($df,$X2);
		return int($p*10000)/10000;
	}

	sub X2df1 {#
		my ($A,$T)=@_;
		return ($A-$T)**2/$T;
	}

	sub X2df2 {#
		my ($A,$T)=@_;
		return (abs($A-$T)-0.5)**2/$T;
	}

sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
	-i	<file>	input vcf file 
	-g	<file>	input group file
	-o	<file>	output file

USAGE
	print $usage;
	exit;
}

