#!/usr/bin/perl
#useage:perl get_the.....pl blast.out > out.file
my %DATA;

while (<>) {
	$_=~s/[\n\r]//g;
	my @arr=split/\t/,$_;
	my $C3=$arr[0];
	my $C4=$arr[11];
	#my ( $C3, $C4 ) = /^\S+\s+\S+\s+(\S+)\s+\(\s*(\S+)/;
    if ( exists $DATA{$C3} ) {
        $DATA{$C3} = [ $C4, $., $_] if $DATA{$C3}[0] < $C4;
    }
    else {
        $DATA{$C3} = [ $C4, $., $_];#[ $C4, $., $_ ]
    }
}

my @LINES = map { $_->[2] } sort { $a->[1] <=> $b->[1] } values %DATA;

for my $ele (@LINES){
  print "$ele\n";
}
#print @LINES;

