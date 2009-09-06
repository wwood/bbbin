#!/usr/bin/perl -w

use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::AlignIO;
use Getopt::Std;

getopt('fpie');

# Takes in a single fasta file with >1 fasta sequence in it and blasts each pair together
$bl2seq_program = 'blastn';
$MIN_EXPECTATION = 10;

$filtering = 0;			# no filtering low complexity regions by default
if ($opt_f){
  $filtering = $opt_f;
  print "filtering: $filtering\n";
}
if ($opt_p){
  $bl2seq_program = $opt_p;
}
if ($opt_e){
  $MIN_EXPECTATION = $opt_e;
}

my $original;
if ($opt_i){
  $original = Bio::SeqIO->new(-file => $opt_i,
			      '-format' => 'fasta');
} elsif ($#ARGV == 0){
  $original = Bio::SeqIO->new(-file => $ARGV[0],
			      '-format' => 'fasta');
} else {
  $original = Bio::SeqIO->new(-fh => \*STDIN,
			      '-format' => 'fasta');
}



while (defined($orig = $original->next_seq())) {
  push @seqs, $orig;
}


foreach $i (0..$#seqs){
  foreach $j ($i+1..$#seqs){

    $first = $seqs[$i];
    $second = $seqs[$j];


    $cur_file = 'blastYourself'.'.'.$i.'.'.$j.'-'.$first->display_id().'.'.
      $second->display_id();

    # Make sure there is no naughty characters in the filename
    $cur_file =~ s/[\|\*\[\]\{\}]//g;

    print "writing bl2seq output file: $cur_file...";
    print "filtering: $filtering\n";
    $factory = Bio::Tools::Run::StandAloneBlast->new('program' => 
						     $bl2seq_program,
						     -outfile => $cur_file,
						     #'F' => $filtering,
						     'e' => $MIN_EXPECTATION
						    );


    $o = length($first->seq());
    $c = length($second->seq());
    if ($o==0 || $c==0) {
      print STDERR "zero length sequence somewhere below:\n";
      print STDERR "orig: $o\n";
      print STDERR "curContig: $c\n\n";
      exit(1);
    }


    #$bl2seq_report = 
    $factory->bl2seq($first,$second);

    print "\n";
  }
}
