#!/usr/bin/perl -w


# Generate a sequence of random fasta files with arbitrary length

if ($#ARGV != 0){
  print "usage: $0 length\n".
    "length is the number of base pairs you want to be in the sequence. Output is piped to STDOUT\n\n";
}

$num = $ARGV[0];
@bases = ('A','G','T','C');

print ">random_sequence_length_$num\n";

foreach (1..$num){
  $i = int rand 4;
  print $bases[$i];
}

print "\n";
