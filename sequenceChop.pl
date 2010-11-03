#!/usr/bin/perl -w

# This script takes in an input fasta sequence, and where to chop it.
# then it outputs the substring of this. I'm sure I've written this code before. eh.


use Bio::SeqIO;
use Getopt::Std;

our($opt_f, $opt_n);
getopts('f:n');

if ($#ARGV != 2) {
  print "usage: sequenceChop.pl [-f <chopping_definitions>] <fasta> <start> <stop>\n";
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
      die "file specifying where to chop is not defined correctly - should be name, start, stop split by tabs. It is `$_'\n";
    }
    my $name = $splits[0];
    $seq_to_start_hash{$name} = int $splits[1];
    $seq_to_stop_hash{$name} = int $splits[2];
  }
}


$original_start = $start;
$original_stop = $stop;
$total_start_bumpers = 0;
$total_end_bumpers = 0;

# Print the chopped file
while ($cur_seq = $original->next_seq()) {

  #Reset the things
  $start = $original_start;
  $stop = $original_stop;
  #Record if the start and end have been hit
  $hit_start = 0;
  $hit_end = 0;

  if ($opt_f){
    my $name = $cur_seq->display_id;
    unless ($seq_to_start_hash{$name}){
      die "Unable to find display id `".$cur_seq->display_id."' in the chopping definition file!\n";
    }
    $start_f = $seq_to_start_hash{$name};
    $stop_f = $seq_to_stop_hash{$name};
    $start = $start+$start_f;
    $stop = $stop+$stop_f;
  }

  if ($start eq '-') {
    $start = 1;
    $hit_start = 1;
  }
  if ($stop eq '-') {
    $stop = $cur_seq->length;
    $hit_end = 1;
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

  # Round off the sequences so they begin and end within the sequence
  if ($start < 1){
    $start = 1;
    $hit_start = 1;
  }
  if ($stop >= $cur_seq->length){
    $stop = $cur_seq->length;
    $hit_end = 1;
  }
  $total_start_bumpers += 1 if $hit_start;
  $total_end_bumpers += 1 if $hit_end;

  print ">".$cur_seq->display_id();
  if ($cur_seq->description) {
    print ' '.$cur_seq->description();
  }
  print "_chopped_";
  print '^' if $hit_start;
  print "$start-$stop";
  print '$' if $hit_end;
  print "\n";


  # If the sequence is less than the stop, then just use the whole sequence
  if ($stop > $cur_seq->length) {
    $stop = $cur_seq->length;
  }


  
  print $cur_seq->subseq($start,$stop);
  print "\n";
}

if ($total_start_bumpers > 0){
  print STDERR "Hit the start on $total_start_bumpers sequence(s)\n";
}
if ($total_end_bumpers > 0){
  print STDERR "Hit the end on $total_end_bumpers sequence(s)\n";
}

