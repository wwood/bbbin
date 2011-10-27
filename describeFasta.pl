#!/usr/bin/perl -w

# This script takes in an input fasta file, and then describes it

use Bio::SeqIO;
#use Getopt::Std;

#read it in from STDIN
my $original = undef;
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
print "Number\tName\tLength\n";

my $separater_string = "\t";

foreach $i (0..$#seqs){
  $j = $i+1;

  print $j.$separater_string.$seqs[$i]->display_id().$separater_string.$seqs[$i]->length()."\n";
}

