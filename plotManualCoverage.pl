#!/usr/bin/perl -w


# Given a coverage info file from consed and a csv telling which positions
# are spaces, output the coverage of each of those
# points


$usage = "Usage: $0 csv coverageFile start stop";

if ($#ARGV != 3){ die "$usage";}

$csv = $ARGV[0];
$cov= $ARGV[1];
$start = $ARGV[2];
$stop = $ARGV[3];


open CSV, "$csv" or die "couldn't open $csv";
open COV, "$cov" or die "couldn't open $cov";


# Define the sequence position hash
$index = 0;
$hash = 0;
@csvs = <CSV>;
$total_snps = $#csvs +1;
foreach $line (@csvs){
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
#foreach $i ($start..$stop+$hash){
#  if (exists $hasher{$i}){
#    print "$i"."::";
#    print "$hasher{$i}\n";
#  }
#}
#print "\n";


# Create reverse hash
while ((@h = each %hasher)){
  $reverse{$h[1]} = $h[0];
}

#print ".$hash..\n";
#foreach $i ($start..$stop){
#  if (exists $reverse{$i}){
#    print "$i"."::";
#    print "$reverse{$i}\n";
#  }
#}
#print "\n";



# Create an array of start and stop positions for the info
# files.

# Example:
# high quality segment from unpadded pos 1 to 5253
#number of reads:  35
#ti|958653299    length 912  102  1013 
#ti|859924275    length 1013  393  1405 
#ti|930653855 C  length 1059  1  1059 
#ti|866279140 C  length 917  188  11
$dummy = <COV>;
$dummy = <COV>;
$i = 0;
foreach $line (<COV>){
  #print $line;
  chomp $line;
  if (!($line =~ m/length\s+(\d+)\s+(-*\d+)\s+(-*\d+)\s+$/)){
    die "Bad line: $line";
  }


  $begin = $2;
  $end = $3;

  if ($begin>$stop || $end<$start){
    next;			# ignore out of bounds ones
  }


#  if (exists $reverse{$2}){
  $trace_starts[$i] = exists $reverse{$2}? $reverse{$2}: $reverse{$start};
#    else
  $trace_stops[$i] = exists $reverse{$3}? $reverse{$3}: $reverse{$stop};


  $i++;



#  if ($#splits == 5){
#    $start_ind = 2;
#    $stop_ind = 3;
#  }
#  elsif ($#splits == 6){
#    $start_ind = 3;
#    $stop_ind = 4;
#  }
#  else {
#    print STDERR "Bad line: $line\n";
#    next;
#  }

#  $trace_starts[$i] = $splits[$start_ind];
#  $trace_stops[$i] = $splits[$stop_ind];
}


#foreach (0..$#trace_starts){
#  print "t: $trace_starts[$_] $trace_stops[$_]\n";
#}
#print "\n";



# Foreach position, go through and work out the number of 
# traces that cover that position.
#
# If it was a space, then just ignore it
#$prev_cov = 0;
#print $start.":::::";
$num = 0;			# how many?
$total = 0;			# how much coverage total?
foreach $pos ($start..$stop+$hash) {

  #print "$pos: ";

  #count how many for each position
  $coverage = 0;
  foreach $i (0..$#trace_starts) {
    #print "s> $i <$trace_starts[$i] $trace_stops[$i]\n";
    #die;
    if ($pos >= $trace_starts[$i] &&
	$pos <= $trace_stops[$i]) {
      $coverage++;
    }
  }

  if (exists $hasher{$pos}){
    print $pos-$start." $coverage\n";

    $num++;
    $total += $coverage;

  }

}


print STDERR "coverage $num $total ".$total/$num
  ." $total_snps ".$total_snps/($stop-$start+1)."\n";
