#!/usr/bin/perl -w


# Given an SNP csv file and a start and end, count the number of SNPs as a percentage

$usage = "$0 csv start stop";

if ($#ARGV != 2){
  die $usage;
}


$csv = $ARGV[0];
$start = int $ARGV[1];
$stop = int $ARGV[2];



open CSV, "$csv";

$count = 0;

foreach $line (<CSV>){
  @splits = split ',',$line;
  $loc = int $splits[0];
  #print ":$loc.$start.$stop:";

  if ($loc > $start && $loc < $stop){
    $count++;
  } else {
    print STDERR "outside: $loc\n";
  }
}


$percent = $count / ($stop-$start+1);
print "$percent $count $start $stop $csv\n";



