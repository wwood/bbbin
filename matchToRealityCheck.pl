#!/usr/bin/perl -w

# Convert the match parsed output into something that can be used by realityCheck.pl


#Example:
#3	12	V$CREL_01	ggagaATTCC	c	(+)	0.862	Match	A,8;T,9;T,10;C,11;C,12
#3	17	F$HSF_03	ggagaattccCTTCC	H	(-)	0.876	Match	C,13;T,14;T,15;C,16;C,17
#8	23	V$ELK1_01	attccCTTCCtgttca	E	(-)	0.876	Match	C,13;T,14;T,15;C,16;C,17


# For the sake of the common denominator with MatInspector, JASPAR, I will not use the end point directly
foreach $line (<STDIN>){
  if ($line =~ m/^\s*$/){next;}

  #print $line;

  @splits = split "\t",$line;
  
  $start = $splits[0];
  $consensus = $splits[3];

  #print "\n>>>>>$start.$consensus\n";
  $finish = $start + length($consensus);

  print "$start\t$finish\n";
}
