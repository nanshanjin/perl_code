use strict;
open IN,"at.vs.at.blastp.1e-3.m8";
my %hash;
open OUT,">out.best.txt";
while(<IN>){
    $_=~s/[\n\r]//g;
    my @array=split/\t/,$_;
    push(@{$hash{$array[0]}},$array[-1]);
    }
close IN;
foreach  my $keys (keys %hash){
    my $max=0;
    for(@{$hash{$keys}}){
        #print $_,"\n";
        if($_ > $max){
            $max = $_;
            }
        #print "$keys\t$max\n";
        }
    print OUT "$keys\t$max\n";
    }
