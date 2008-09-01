#!/usr/bin/perl -w

open NAMES, "$ENV{HOME}/thesis/pwm/Match/parsedMatchReducedSet.txt";
@nam = <NAMES>;
chomp @nam;


#$a = 'FOXJ_';
#@b = ('FOXJ_','fdsa');
#print grep (/$a/, @b);
#print "\n..\n";
open APP, ">$ENV{HOME}/thesis/doc/thesis/reducedSetAppendix.tex";

print APP '\chapter{Reduced Set}'."\n"
  .'\label{sec:reducedSetAppendix}'."\n"
  .'\begin{longtable}{|c|c|}'."\n"
  #.'\caption{Summary of reduced TRANSFAC set}'."\n"
  #.'\centering'."\n"
  .'\hline'."\n"
  .'\bfseries Family & \bfseries Name\\\\\\hline'."\n\n";






foreach $line (@nam) {
  chomp $line;
  $line =~ s/\&/\\&/g;


  @splits = split "\t",$line;

  if ($#splits > 1 && !($line=~m/uperclass/)) {
    
    $name = $splits[1];
    $name =~ s/;/\\\\&/g;
    print APP "\n".$splits[0].'&'.$name.'\\\\\\hline';
  }


}


print APP '\hline'."\n\n"
  .'\end{longtable}'."\n";
