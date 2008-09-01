#!/usr/bin/perl -w


#use strict;

$usage = "$0 <MatchParsed> <MatInspectorParsed>\n".
  "Output is piped to STDOUT\n\n";

if ($#ARGV != 1){print $usage;exit;}
$MatInspectorFile = pop @ARGV;
$MatchFile = pop @ARGV;





#Read in the Match file into an array

#1	12	AREB6_02	ttaaAGGTGaat	A	(-)	0.987	Match	A,5,G,6,G,7,T,8,G,9
#9	20	N$SKN1_02	gaatATCATgcc	S	(+)	0.986	Match	A,13,T,14,C,15,A,16,T,17

open MATCH_IN, "$MatchFile" or die "could no open Match file: $MatchFile";
$i=0;


		 
foreach $line (<MATCH_IN>){

  if ($line =~ m/^\s*$/){next;}	# ignore blank lines

  #else read in a useful line
  $line =~ s/\n//;
  @splits = split "\t",$line;

  if ($#splits != 8){die "unrecognised line in Match file: $line";}
  
  #else {print $splits[8]."\n";}

  push @matchMatchesStartPoses, $splits[0];
  push @matchMatchesAccIds, $splits[2];
  push @matchMatchesConfidences, $splits[6];
  push @matchMatchesCores, $splits[8];
}
close MATCH_IN;






#Read in MatInspector and match them to the Match ones

#1	17	P$GTBX/GT1.01	ttaaagGTGAatatcat	GT1-Box...	(+)	0.85	MatInspector	G,7,T,8,G,9,A,10
#1	13	V$ZFHX/AREB6.02	ttaaaGGTGaata	AREB6....ctor 6)	(-)	0.97	MatInspector	G,6,G,7,T,8,G,9

open MATINSPECTOR_IN, "$MatInspectorFile" or die "could not open MatInspector file: $MatInspectorFile";
foreach $line (<MATINSPECTOR_IN>){
  if ($line =~ m/^\s*$/){next;}	# ignore blank lines

  #else read in a useful line
  $line =~ s/\n//;
  @splits = split "\t",$line;

  if ($#splits != 8){die "unrecognised line in MatInspector file: $line";}

  #search through the Match hits to see if any are the same
  $mati_core = $splits[8];
  @coreSplits = split ';',$mati_core;

  #if each element of the current line's cores is a single Match core
  foreach $i (0..$#matchMatchesCores){
    $errors = 0;

    $mcore = $matchMatchesCores[$i];
    foreach $elem (@coreSplits){
      if (!($mcore =~ m/$elem/)){
	$errors++;
      }
    }


    # Make sure the first letter of the Id's is the same
    $MatchOffset = 2;
    $MatInspectorOffset = 7;
    $MatchFirst = substr $matchMatchesAccIds[$i], $MatchOffset, 1;
    $MatInspectorFirst = substr $splits[2], $MatInspectorOffset, 1;
    if (!($MatchFirst eq $MatInspectorFirst)){
      $errors++;
    }


    #print $errors.".";
    if ($errors == 0){
      $sepChar = "\t";
      print $splits[2]
	."\n     ".$matchMatchesAccIds[$i]
	.$sepChar.$matchMatchesStartPoses[$i]
	.$sepChar.$splits[0]
#	.$sepChar.$matchMatchesNames[$i]
	.$sepChar.$matchMatchesConfidences[$i]
	.$sepChar.$splits[6]
	.$sepChar.$splits[4]
	  ."\n\n";
      #exit;
    }


  }
}
