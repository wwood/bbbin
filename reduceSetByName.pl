#!/usr/bin/perl -w

open NAMES, "$ENV{HOME}/thesis/pwm/Match/parsedMatchReducedSet.txt";
@nam = <NAMES>;
chomp @nam;


#$a = 'FOXJ_';
#@b = ('FOXJ_','fdsa');
#print grep (/$a/, @b);
#print "\n..\n";




foreach $line (<STDIN>){
  chomp $line;


  @splits = split "\t",$line;
  #print $splits[]."mm\n";
  #print $mats[0]."nn\n";


  @matches = grep /$splits[4]/, @nam;

  if ($#matches != -1){
    print $line."\n";
    foreach (@matches){
      print "\t".'n.'.$_.".\n";
    }
  }
}
