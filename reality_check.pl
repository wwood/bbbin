#!/usr/bin/perl -w

# This script takes a number of pwm hit records, in the form:

# <start> <finish>
# And then it converts the list of pwm hits to a bigger list of
# numbers to be plotted as a frequency plot

$FIRST = 0;
$LAST = 2500;

foreach ($FIRST..$LAST){
  $points[$_] = 0;
}

foreach $line (<STDIN>){
  @splits = split "\t",$line;
  
  $start = $splits[0];
  $finish = $splits[1];

  foreach $i ($start..$finish){
    $points[$i]++;
  }
}


foreach $i ($FIRST..$LAST){
  print "$points[$i]\n";
}
