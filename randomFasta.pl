#!/usr/bin/perl -w


# Generate a sequence of random fasta files with arbitrary length

if ($#ARGV != 0 and $#ARGV != 1){
  print "usage: $0 length [num]\n".
    "length is the number of base pairs you want to be in the sequence, num is the number of them (default 1). Output is piped to STDOUT\n\n";
}

$num = $ARGV[0];
@bases = ('A','G','T','C');

$num_total = 1;
if ($#ARGV == 1){
    $num_total = $ARGV[1]+0;
}

foreach(1..$num_total){
    print ">random_sequence_length_$num"."_$_\n";
    foreach (1..$num){
        $i = int rand 4;
        print $bases[$i];
    }
    print "\n";
}
