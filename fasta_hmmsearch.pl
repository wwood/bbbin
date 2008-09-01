#!/usr/bin/perl -w

use Getopt::Std;

# A script to take a fasta file of protein sequences as a fasta file,
# create the HMM from that, and then apply that HMM to a seqfale,
# all with HMMER.

getopt('f');

my $fasta_file;
if ($opt_f) {
  $fasta_file = $opt_f;
  if ($#ARGV == 0){
    $database = $ARGV[0];
  } else {
    &usage();
  }
} elsif ($#ARGV == 1) {
  $fasta_file = $ARGV[0];
  $database = $ARGV[1];
} else {
  &usage();
}


#create the temp files
# I tried to use the perl tempfile library, but I can't seem to
# get just the name of it, I always have to actually create the
# file. And that isn't what I want.
my $tmp_filename = 'clustalw.temp.'.rand;
my $hmm_filename = 'hmm.temp.'.rand;
print `clustalw -infile=$fasta_file -outfile=$tmp_filename`;


#create the HMM
print `hmmbuild $hmm_filename $tmp_filename`;
print `hmmsearch $hmm_filename $database`;



sub usage {
  die "Usage: $0 <fasta file> <database>\n\n";
}
