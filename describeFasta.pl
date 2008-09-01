#!/usr/bin/perl -w

# This script takes in an input fasta file, and then describes it

use Bio::SeqIO;
#use Getopt::Std;

#read it in from STDIN
if ($#ARGV==-1){
  $original = Bio::SeqIO->new(-fh => \*STDIN,
			      '-format' => 'fasta');
} elsif ($#ARGV==0) {
  $original = Bio::SeqIO->new(-file => $ARGV[0],
			      '-format' => 'fasta');
}

$i=0;
while(defined($cur = $original->next_seq()))
{
  #print out name, then reversed sequence
  $seqs[$i] = $cur;
  $i++;
}

# Print the number of sequences in the input file
# Print in media-wiki table format
print '{| border="1" cellspacing="0" cellpadding="5" align="center"';
print "\n";
print "!Number\n";
print "!Name\n";
print "!Length\n";

foreach $i (0..$#seqs){
  print "|-\n";

  print "|";
  $j = $i+1;
  print "$j";
  print "\n";

  print "|";
  print $seqs[$i]->display_id();
  print "\n";

  print "|";
  print $seqs[$i]->length();
  print "\n";

  #print $seqs[$i]->subseq(1,10);
  #print "\n";
}

print "|}\n";

#{| border="1" cellspacing="0" cellpadding="5" align="center"
#!Name
#!Length
#|-
#|Contig0.1      
#|2693
#|-


