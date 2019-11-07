#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

open (OUT, ">snp.file")|| die "can not open the text.file\n";

my %data;
my @files = glob "*.vcf";
for my $file (@files) {
        my $file_name = basename ($file);
		print "$file_name\n";
		open my $fh, '<', $file or die $!;
		while (<$fh>) {
                next if /^\s*$/;
                next if /^gene_number/;
				next if /^#/;
				
			    $_ =~s/[\n\r]//g;
	            my @array = split(/\t/,$_);
				my ($num, $pos, $snp) =($array[0],$array[1],$array[3]."|".$array[4]);
				#print "$num\t$pos\t$snp\n";
                chomp;
                #my ($num, $pos, $snp) = split;
                $data{$num}{$pos}{$file_name} = $snp;
        }
        close $fh;
}
  print OUT"reference\tposition";
for my $file (sort @files) {
        my $file_name = basename ($file);
		my @new_file_name=split/_/,$file_name;
        print OUT"\t$new_file_name[0]";
}
print OUT"\n";

   for my $num (sort keys %data) {
        for my $pos (sort keys %{$data{$num}}) {
                print OUT"$num\t$pos";
                for my $file (sort @files) {
                        my $file_name = basename ($file);
                        if (exists $data{$num}{$pos}{$file_name}) 
						{
                                print OUT"\t$data{$num}{$pos}{$file_name}";
                        } else 
						{
                                print OUT"\t.";
                        }
                }
                print OUT"\n";
        }
}      