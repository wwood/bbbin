#!/usr/bin/perl -w

# Removes the orientation from a gff file, primarily useful
# so that gff2ps plots everything together.

use Getopt::Std;
use Bio::FeatureIO;

our ($opt_m, $opt_r);
getopts('mr');




# Modify the script so that orientation info is encoded into the 
# feature part, as opposed to the orientation part, so it fits better
# with gff2ps
foreach $line (<>){

  if (&is_comment_line($line)){print $line; next;}


  @bits = split(/\s+/, $line);

  # if negative, make it positive and change the type
  if (!$opt_r and $bits[6] eq '-') {
    $bits[2] = $bits[2].'Reversed';
  }

  #opt_r makes everything positive but puts Reversed
  elsif ($opt_r and $bits[6] eq '+') {
    $bits[2] = $bits[2].'Reversed';
  }

  # Change all orientations to forward
  $bits[6] = '+';

  $main = join "\t", @bits[0..7];
  if ($#bits > 7) {
    $end = join ' ', @bits[8..$#bits];
    print $main."\t".$end."\n";
  } else {
    print $main."\n";
  }
}



# Return 1 if the line is a comment, else 0
sub is_comment_line {
  my ($line) = @_;

  if ($line =~ m/^\s*$/ or $line =~ m/^\s*\#/){
    return 1;
  } else {
    return 0;
  }
}
