#!/usr/bin/perl -w


# Given an set of reduced matrices, and a file giving the relationships between
# the matrices and the names, output a list of matrices and names and groups and classifications
# that are ok to be used.


open REDUCED_SET, "$ARGV[0]";
@reduced_set = <REDUCED_SET>;
grep s/\$/\\\$/g,@reduced_set;

foreach $line (<STDIN>) {

  #Given a line from the predictions, it
  # is only good if the matrix greps to the one in the reduced set
  #print $line;

  #ignore certain lines
  if ($line =~ m/ignored:/){
    next;
  }

  chomp $line;
  $line =~ s/\$/\\\$/g;

  @splits = split "\t",$line;

  if ($#splits != 2){
    print STDERR "strange line: .$line.\n";
    exit;
  }


  @matches = grep /$splits[1]/, @reduced_set;




  print $line;
  if ($#matches != -1){
    print "\n";
    foreach (@matches) {
      print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\t$_";
    }
  }
  print "\n";
  #exit;

}
