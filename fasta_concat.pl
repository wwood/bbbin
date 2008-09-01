#!/usr/bin/perl -w

#Given a fasta file with multiple sequences, concatenate it to stdout
# and print where all the original sequences start in the new sequence

use strict;
use Bio::SeqIO;

my $start = 1;

my $seqio;
my $outfilename;
if ($#ARGV==0){
  $seqio = Bio::SeqIO->new(-fh => \*STDIN,
			      '-format' => 'fasta');
  $outfilename = $ARGV[0];
} elsif ($#ARGV==1) {
  $seqio = Bio::SeqIO->new(-file => $ARGV[0],
			      '-format' => 'fasta');
  $outfilename = $ARGV[1];
} else {
  print STDERR "usage: $0 <input_seq> <output_seq>\n or the input_seq can be STDIN as well";
  exit;
}


my $cumulative = '';
while(defined(my $cur = $seqio->next_seq()))
{
  print $cur->display_name.' '.$cur->description." $start\n";
  $start += $cur->length;
  $cumulative = $cumulative.$cur->seq;
}


#Write the sequence
open OUT, ">$outfilename" or die "failed to open '$outfilename' for writing";
print OUT '>'.$outfilename."_concatenated\n";
print OUT $cumulative."\n";
close OUT;

