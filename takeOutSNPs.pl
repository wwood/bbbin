#!/usr/bin/perl -w

# Given a gnuplot file corresponding to quality values,
# take out those qualities that are consistent with
# the assembler predicting a SNP and output where those
# positions so they can be graphed

$count = 0;
foreach (<STDIN>){
  @splits = split ' ';
  
  push @indices, $splits[0];
  
  $qual = $splits[1];
  if ($qual == -97){
    push @quals, -10;
  }
  else {
    push @quals, -1;
    $count++;
  }
}


foreach ($i=0; $i<=$#indices; $i++){
  print $indices[$i].' '.$quals[$i]."\n";
}

print STDERR "$count ".$#indices+1.' '.$count/($#indices+1)."\n";
