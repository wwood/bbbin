#!/usr/bin/perl -w


# Take in a parsed PWM output, and return the hits that are in the last part of it
# STDOUT - number and ratio of hits
# STDERR - the hits themselves

$length = $ARGV[0];
$top = $ARGV[1];


if ($top > $length){
#    $top = $length;
}


$cutoff = $length-$top;

foreach $line (<STDIN>){
    @splits = split ' ', $line;

    if ($splits[0] > $cutoff && $splits[0]<$length){
	push @hits, $line;
    }
}


foreach (@hits){
    print STDERR "$_";
}


$num = $#hits+1;
$ratio = $num/$top;
print "$num $ratio\n";
