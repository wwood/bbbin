#!/usr/bin/perl -w

# takes in a parsed pwm file, and returns an array, so that a graph of
# confidences and number of predictions can be drawn.

$CONFIDENCE_INDEX = 6;
$FIRST_CONF = 0.8;
$CONF_BINS = 0.0025;
$LAST_CONF = 1;

# Read the confidences into memory
foreach $line (<STDIN>){

  @splits = split "\t",$line;
  push @confidences, $splits[$CONFIDENCE_INDEX];
}



for ($conf=$FIRST_CONF; $conf<$LAST_CONF; $conf = $conf+$CONF_BINS){
  $count = 0;
  foreach $c (@confidences){
    if ($c > $conf){
      $count++;
    }
  }
  print "$conf $count\n";
}

