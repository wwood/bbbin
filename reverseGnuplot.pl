#!/usr/bin/perl -w

# Given a sequence of newline separated records of the 
# form "<index> <value>", reverse the order of them, so
# that the largest index becomes the first index, and the
# first index becomes the largest index. It is assumed that
# the first index is the smallest.

foreach (<STDIN>){
  @splits = split ' ',$_;
  push @indices, $splits[0];
  push @values, $splits[1];
}

$smallest = $indices[0];
$largest = $indices[$#indices];

#print "$largest . $smallest\n";
for ($i=$#indices; $i>=0; $i--){
  print $indices[$i]-$largest+$smallest+2*($#indices-$i)." $values[$i]\n";
}
