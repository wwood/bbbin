#!/usr/bin/perl -w

# Given the output from pwmConsensus.pl, which combines the output from Match and MatInspector,
# add in the ones where JASPAR agrees

#MA0071 ,RORA ,8.125 ,0.827722151132543 ,1 ,10 ,1 ,TTAAAGGTGA 
#MA0020 ,Dof2 ,6.086 ,0.8684839414205 ,3 ,8 ,1 ,AAAGGT 


use promoterAnalysis;

$usage = "$0 <JASPARcsv> <Consensus>\n\n";
if ($#ARGV != 1){die $usage;}

$ConsensusFile = pop @ARGV;
$JASPARFile = pop @ARGV;

open JASPAR_IN, "$JASPARFile" or die "could not open JASPAR csv file: $JASPARFile";
$i=0;
foreach $line (<JASPAR_IN>){
  @splits = split ",",$line;

  foreach (@splits){
    s/\s*//g;
    s/\W*//g;
  }

  $JasparAccIds[$i] = $splits[0];
  $JasparNames[$i] = $splits[1];
  $JasparConfidences[$i] = $splits[3];
  $JasparStartPoses[$i] = $splits[4];

  $JasparCores[$i] = promoterAnalysis::getCore($splits[7],$JasparStartPoses[$i]);

  $i++;
}

print "$#JasparCores $JasparCores[1]\n";


#Read in the consensus stuff, and try to align them
open CONSENSUS_IN, "$ConsensusFile" or die "could not open Consensus file: $ConsensusFile";

#For each line of the consensus, try to match each Jaspar hit to it.
foreach $line (<CONSENSUS_IN>){
  @splits = split "\t",$line;
  

  $i=0;
  foreach $jline (@JasparCores){

    if (promoterAnalysis::coreMatches($jline, $line)){
      $sepChar = "\t";
      print $JasparNames
	.$sepChar .$line
	.$sepChar .$JasparConfidences[$i]
	.$sepChar .$JasparAccIds[$i]
	."\n";
    }
    $i++;
  }
}
