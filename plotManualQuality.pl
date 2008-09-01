#!/usr/bin/perl -w

# Given a quality file and a csv file telling which positions
# are spaces, output the quality in a form readable by gnuplot 
# showing the quality, making sure it is compatible with the
# new sequence positions


# Define that the coverage of a space is the same as the
# coverage of the previous element


$usage = "Usage: $0 csv qualityFile start stop";

if ($#ARGV != 3){ die "$usage";}

$csv = $ARGV[0];
$qual= $ARGV[1];
$start = $ARGV[2];
$stop = $ARGV[3];


open CSV, "$csv" or die "couldn't open $csv";
open QUAL, "$qual" or die "couldn't open $qual";


# Define the sequence position hash
$index = 0;
$hash = 0;
foreach $line (<CSV>){
  @splits = split ',',$line;

  while ($index+$start<=$splits[0]){
    $hasher{$start+$index+$hash} = $start+$index;
    #print "a: $start+$index\n";
    $index++;
  }


  if ($splits[1] eq 'between'){
    #print "index: $index\n";
    $hash++;
    #print "hash: $hash\n";
  }

  else {
    $hasher{$start+$index+$hash} = $start+$index;
    
    #print "b: $start+$index\n";
    $index++;
  }
}

# Finish off until we get to the end of the part required.
# Assume that there are no csv lines specifying between's outside the region
while ($index+$start<=$stop){

  #print "c: $start+$index\n";
  $hasher{$start+$index+$hash} = $index+$start;
  $index++;
}

#print ".$hash..\n";
#foreach $i ($start..$stop){
#  print "$i"."::";
#  print "$hasher{$i}\n";
#}
#print "\n";


# Write out the quality info as -(97-qual), so that higher
# values mean lesser quality
$first = <QUAL>;
$first = $first;
@all = <QUAL>;
chomp @all;
$qualslist = join ' ',@all;
#print $qualslist."\n";
@quals = split ' ',$qualslist;

foreach ($start..$stop+$hash){

  $i = $_;

  if (exists $hasher{$i}){

    $local = $hasher{$i};
    $pos = $i-$start;
    $q = -(97-$quals[$local-1]);

    # So the 0 values are visible
    if ($q == 0){
      $q = -0.1;
    }

  } else {
    $pos = $i-$start;
    $q = -(98-0)
  }
  print "$pos $q\n";


}


