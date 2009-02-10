#!/usr/bin/perl -w


# This file reads in a fasta sequence, and then pipes it out, with each of the sequences in it reverse complemented.

use Bio::SeqIO;


#read it in from STDIN
my $original;
if ($#ARGV == 0){
  $original = Bio::SeqIO->new(-file => $ARGV[0],
			      '-format' => 'fasta');
} else {
  $original = Bio::SeqIO->new(-fh => \*STDIN,
			      '-format' => 'fasta');
}

while(defined($cur = $original->next_seq()))
{
  #print out name, then reversed and complemented sequence
  print ">";
  print $cur->id();
  print "_revcom\n";
  print $cur->revcom()->seq();
  print "\n";
}
