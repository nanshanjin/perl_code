#!/usr/bin/perl -w
use strict;
use SVG;

#explanation:this program is edited to 
#edit by hewm;   Wed Jan 20 11:04:17 CST 2016
#Version 1.0    hewm@genomics.org.cn 

die  "Version 1.0\t2016-01-20;\nUsage: $0 <InPut><Out>\n" unless (@ARGV ==2);

#############Befor  Start  , open the files ####################

open (IA,"$ARGV[0]") || die "input file can't open $!";
open (IB,"$ARGV[1]") || die "input file can't open $!";

#open (OA,">$ARGV[1]") || die "output file can't open $!" ;

################ Do what you want to do #######################
my %scaf_length=();
my %scaf_Col=();
my %chr_max_CM=();
my %chr_max_PY=();

my %chr_scaffold=();

	while(<IA>) 
	{ 
		chomp ; 
		next if  ($_=~s/#/#/g);
		my @inf=split ;
		$scaf_length{$inf[0]}=$inf[2];
		$scaf_Col{$inf[0]}=$inf[3];
#		$chr_max_CM{$inf[4]}=$inf[-1];
		if (!exists $chr_max_PY{$inf[4]})
		{
			 $chr_max_PY{$inf[4]}=$inf[-1];
		}
		elsif ( $chr_max_PY{$inf[4]}  < $inf[-1] )
		{
			$chr_max_PY{$inf[4]}=$inf[-1];
		}
		
		if  (!exists  $chr_scaffold{$inf[4]})
		{
	    	$chr_scaffold{$inf[4]}=$inf[0];
		}
		else
		{
		  $chr_scaffold{$inf[4]}=$chr_scaffold{$inf[4]}."\t".$inf[0];
		}
	}
close IA;

#MC07    scaffold140_3171214     259.593 scaffold140     3171214 +       MC07    30671539
$chr_max_CM{"MC01"}=213.393;	
$chr_max_CM{"MC02"}=162.413;	
$chr_max_CM{"MC03"}=184.852;	
$chr_max_CM{"MC04"}=186.191;	
$chr_max_CM{"MC05"}=175.982;	
$chr_max_CM{"MC06"}=257.19;	
$chr_max_CM{"MC07"}=135.263;	
$chr_max_CM{"MC08"}=275.343;	
$chr_max_CM{"MC09"}=243.139;	
$chr_max_CM{"MC10"}=199.521;	
$chr_max_CM{"MC11"}=170.662;	

my $max_CM=0;
my $max_PY=0 ;
#34593500;
foreach my $k  (keys  %chr_max_CM )
{
	if  ($chr_max_CM{$k} > $max_CM )
	{
		$max_CM=$chr_max_CM{$k};
	}
	if  ($chr_max_PY{$k} > $max_PY )
	{
		$max_PY=$chr_max_PY{$k};
	}
}

my ($up, $down, $left, $right) = (120, 100, 150, 100);

my $X_bin=120;

my $Xlength=1320 ;
my $YLength=1000;

my $Y1_Bin=$YLength/$max_CM;
my $Y2_Bin=$YLength*1.0/$max_PY;
   $Y2_Bin= $Y2_Bin*0.7;	

my $height = $up + $down +  $YLength;
my $width = $left + $right + $Xlength ;
my $fontsize = 25;
my $cirsize=2;
my $spacelc = 'grey';# space line color


#close OA ;

my %Chr_X=();


my $svg = SVG->new('width',$width,'height',$height);

my $count=1;
foreach my $chr  (sort keys  %chr_max_CM )
{
	my $X_This=$right+$X_bin*$count;
	my $HH=$Y1_Bin*$chr_max_CM{$chr};
	$svg->text('text-anchor','middle','x',$X_This,'y',$up-20,'-cdata',"$chr",'stroke','black','font-family','Arial','font-size',$fontsize/2);
    $svg->rect('x',$X_This,'y',$up,'width',5,'height',$HH,'fill',$spacelc,,'stroke',$spacelc);
	$Chr_X{$chr}=$X_This;
	$count++;
}

my %Sca_X=();
my %Sca_Y=();
my $sca_bin=10;

$count=1;

foreach my $chr  (sort keys  %chr_max_PY )
{
	my $X_This=$right+$X_bin*$count+$X_bin/3 ;
	my $YY_pp=$up;
	#print "$chr\n"; next;
	$svg->text('text-anchor','middle','x',$X_This,'y',$up-20,'-cdata',"scaffold",'stroke','black','font-family','Arial','font-size',$fontsize/2);
#	print "$chr\t",$chr_scaffold{$chr},"\n";
	my @scaName=split /\s+/,$chr_scaffold{$chr};
	foreach my $kk (0..$#scaName)
	{
		my $ID=$scaName[$kk];
		my $Y_this=$YY_pp;
		if  ($kk!=0)
		{
			$Y_this+=$sca_bin;
		}
    	
		my $HH=$Y2_Bin*$scaf_length{$ID};
		my $col="blue";
		if ($scaf_Col{$ID} eq "+"  ||  $scaf_Col{$ID} eq "-")
		{
			$col="green";
		}
		#$svg->text('text-anchor','middle','x',$X_This,'y',$up-20,'-cdata',"$chr",'stroke','black','font-family','Arial','font-size',$fontsize);
	    $svg->rect('x',$X_This,'y',$Y_this,'width',5,'height',$HH,'fill',$col,'stroke',$spacelc);
		$Sca_X{$ID}=$X_This;
		$Sca_Y{$ID}=$Y_this;

		$YY_pp=$Y_this+$HH;
	}
	$count++;
}

#exit(1);
#<IB>;
	while(<IB>)
	{ 
		chomp ; 
		my @inf=split ;
		my $sca_site=$inf[4];
		my $sca_namne=$inf[3];
		my $chr=$inf[6];
		if ($scaf_Col{$sca_namne} eq "-")
		{			
			$sca_site=$scaf_length{$sca_namne}-$sca_site;
		}
		
		my $Thix_Chr_X=$Chr_X{$chr}+5;
		my $Thix_Sca_X=$Sca_X{$sca_namne};

		my $Thix_Chr_Y=$up+$Y1_Bin*$inf[2];

		my $Thix_Sca_Y=$Sca_Y{$sca_namne}+$Y2_Bin*$sca_site;
		$svg->line('x1',$Thix_Chr_X,'y1',$Thix_Chr_Y,'x2',$Thix_Sca_X,'y2',$Thix_Sca_Y,'stroke','black','stroke-width',0.5);
	
	}
close IB ;

print $svg->xmlify();

######################swimming in the sky and flying in the sea ###########################
