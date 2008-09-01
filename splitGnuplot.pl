#!/usr/bin/perl -w

# Given a gnuplot input and a cutoff, split it into 2 parts:
# those before the cutoff and those on or after the cutoff

$cutoff = $ARGV[0];
$out1 = $ARGV[1];
$out2 = $ARGV[2];


open OUTA, ">$out1";
open OUTB, ">$out2";


foreach $line (<STDIN>){
  @splits = split ' ',$line;

  if ($splits[0]<$cutoff){
    print OUTA $line;
  } else {
    print OUTB $line;
  }
}
