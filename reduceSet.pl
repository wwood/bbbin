#!/usr/bin/perl -w

open MATRIX, "$ARGV[0]";
@mats = <MATRIX>;
chomp @mats;


#$a = 'FOXJ_';
#@b = ('FOXJ_','fdsa');
#print grep (/$a/, @b);
#print "\n..\n";

foreach $line (<STDIN>){
  chomp $line;

  #dollar signs mess things up with the matching and
  #I don't know how to make them be treated as the
  #characters that they are
  $line =~ s/\$/\\\$/g;

  @splits = split "\t",$line;
  #print $splits[2]."mm\n";
  #print $mats[0]."nn\n";

  grep s/\$/\\\$/g,@mats;

  @matches = grep /$splits[2]/, @mats;

  if ($#matches != -1){
    print $line."\n";
    foreach (@matches){
      print "\t".'n.'.$_.".\n";
    }
  }
}
