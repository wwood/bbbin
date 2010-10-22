#!/usr/bin/perl -w

# This script takes in an input fasta sequence, and where to chop it.
# then it outputs the substring of this. I'm sure I've written this code before. eh.


use Bio::SeqIO;
use Getopt::Std;

our($opt_f, $opt_n);
getopts('f:n');

if ($#ARGV != 2) {
  print "usage: sequenceChop.pl <fasta> <start> <stop>\n";
  print "output is piped to command line. 1-2 are the first two
           bases of the sequence.\n";
  print "if <fasta> is '-', then input is from STDIN, not a file.\n";
  print "if <start> is '-', then 1 is used.\n";
  print "if <stop> is '-', then the max is used.\n\n";
  print "Negative values are also allowable, -1 indicates max (last position)\n\n";
  exit(0);
}

$in_fasta = $ARGV[0];
$start = $ARGV[1];
$stop = $ARGV[2];



if ($in_fasta eq '-') {
  $original = Bio::SeqIO->new(-fh => \*STDIN,
			      '-format' => 'fasta');

} else {
  $original = Bio::SeqIO->new(-file => $in_fasta,
			      '-format' => 'fasta');
}


# If opt_f, read in the file specifying where each sequence should be chopped.
my %seq_to_start_hash = ();
my %seq_to_stop_hash = ();
if ($opt_f){
  open IN, $opt_f or die "Could not open file specifying where each sequence should be chopped";
  while (<IN>){
    chomp;
    my @splits = split /\t/, $_;
    unless ($#splits == 2){
      die "file specifying where to chop is not defined correctly - should be name, start, stop split by whitespace. It is `$_'\n";
    }
    my $name = $splits[0];
    $seq_to_start_hash{$name} = int $splits[1];
    $seq_to_stop_hash{$name} = int $splits[2];
  }
}


$original_start = $start;
$original_stop = $stop;

# Print the chopped file
while ($cur_seq = $original->next_seq()) {

  #Reset the things
  $start = $original_start;
  $stop = $original_stop;

  if ($opt_f){
    my $name = $cur_seq->display_id;
    unless ($seq_to_start_hash{$name}){
      die "Unable to find display id `".$cur_seq->display_id."'!\n";
    }
    $start_f = $seq_to_start_hash{$name};
    $stop_f = $seq_to_stop_hash{$name};
    $start = $start+$start_f;
    $stop = $stop+$stop_f;
  }


  if ($start eq '-') {
    $start = 1;
  }
  if ($stop eq '-') {
    $stop = $cur_seq->length;
  }
  # Convert negative values into positive ones since bioperl can't handle it
  if ($start < 0) {
    if ($opt_n){
      $start = 1;
    } else {
      $start = $cur_seq->length+$start+1;
    }
  }
  if ($stop < 0) {
    if ($opt_n){
      $stop = $cur_seq->length;
    } else {
      $stop = $cur_seq->length+$stop+1;
    }
  }

  if ($cur_seq->description) {
    print ">".$cur_seq->display_id().' '.$cur_seq->description()."_chopped_$start-$stop\n";
  } else {
    print ">".$cur_seq->display_id()."_chopped_$start-$stop\n"
  }

  # If the sequence is less than the stop, then just use the whole sequence
  if ($stop > $cur_seq->length) {
    $stop = $cur_seq->length;
  }


  
  print $cur_seq->subseq($start,$stop);
  print "\n";
}
