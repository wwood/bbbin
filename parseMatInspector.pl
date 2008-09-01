#!/usr/bin/perl -w

# This script parses in a MatInspector results table from csv created by an openoffice spreadsheet middleman.

$usage = "usage: $0 <csvFile> <SequenceName>\n".
  "output is piped to STDOUT\n\n";

use promoterAnalysis;

open INPUT, "<$ARGV[0]" or die "could not open input file\n\n$usage";
$seqName = $ARGV[1];


$started = 0;
foreach $line (<INPUT>) {

  # If not up to the correct table yet, look for the start
  if ($line =~ m/Inspecting sequence .* \[(.*)\]/) {

    if ($1 eq $seqName) {
      #print "started on \"$line\n\n";
      $started = 1;
      next;
    } else {
      #print "finished on \"$line\"\n";#stop reading out entries
      $started = 0;
    }
  }

  #print "$#_\n";

  if ($started==1) {

    $line =~ s/\n//;

    #parse values properly with comma separation
    @values = split ',',$line;

    if (defined($values[1]) && $values[1] =~ m/\$/) {
      #print "-----------------$values[1]---\n";

      foreach (@values) {
	s/\"//g;
      }


      $offset = $#values-16;	 # to correct for extra commas in the name

      $AccId = $values[1];
      $Name = $values[4+$offset];
      $Confidence = $values[6+$offset];

      $pos = $values[8+$offset];
      #print "pos: .$pos.\n";
      $pos =~ m/(\d+).*?(\d+)/;
      $StartPos = $1;
      $FinishPos = $2;

      $Direction = $values[11+$offset];

      $Consensus = $values[16+$offset];
      if ($Direction eq '(-)') {
	$Consensus = promoterAnalysis::revcom_keepcase($Consensus);
      }

      $Program = "MatInspector";

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

    }
  }
}
