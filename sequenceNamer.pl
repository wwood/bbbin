#!/usr/bin/perl -w

# Each line coming is assumed to be a nucleotide or protein sequence.
# The output is a fasta file with sequences named consecutively seq1, seq2,
# etc.

my $count = 1;
foreach (<>){
  print '>seq'.$count."\n".$_;
  $count = $count+1;
}
