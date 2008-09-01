#!/usr/bin/perl -w

# Parse input from Match

#Example:
#V$OCT1_Q6                    1 (+)  0.797  0.674  ttaaagGTGAAtatc                    Oct-1
#V$AREB6_02                   1 (-)  1.000  0.987  ttaaAGGTGaat                       AREB6


$usage = "usage: $0 <csvFile>\n".
  "output is piped to STDOUT\n\n";

use promoterAnalysis;

open INPUT, "<$ARGV[0]" or die "could not open input file\n\n$usage";


foreach $line (<INPUT>) {
  chomp $line;
  if (!($line =~ m/^(\S+)\s+(\d+)\s+(\([+-]\))\s+(\d.\d+)\s+(\d.\d+)\s+(\w+)\s+(.*)/)){
    print STDERR "Bad line $line\n";
    #exit;
  }

  #print $line;
  $AccId = $1;
  $StartPos = $2;
  $Direction = $3;
  $Confidence = $5;
  $Consensus = $6;
  $Name = $7;
  $Program = "Match";

  $FinishPos = $StartPos + length($Consensus) -1;

  $sepChar = "\t";
  print $StartPos
    .$sepChar. $FinishPos
      .$sepChar. $AccId
	.$sepChar. $Consensus
	  .$sepChar. $Name
	    .$sepChar. $Direction
	      .$sepChar. $Confidence
		.$sepChar. $Program
		  .$sepChar;

  promoterAnalysis::printCore($Consensus, $StartPos, ',');
  print "\n";
  #exit;

}
