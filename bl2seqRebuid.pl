#!/usr/bin/perl -w

# THis compares 2 sequences and then produces a report of them
# input are the original fasta file, which is piped in,
# and the contigs file from the rebuilding stage

use Getopt::Std;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::AlignIO;
#use Bio::Tools::Run::RemoteBlast;
#use strict;




# PARSE ARGUMENTS
getopt('fspoeij');
#print "f: $opt_f\n";
#print "s: $opt_s\n";
#print "p: $opt_p\n";
#print "o: $opt_o\n";
#print "i: $opt_i\n";
#print "j: $opt_j\n";

if (!$opt_s && !opt_j) {
  &usage();
}

$firstfile = "../original_proteins.fasta";
if ($opt_f){
  $firstfile = $opt_f;
} elsif ($opt_i) {
  $firstfile = $opt_i;
}

if ($opt_j){
  $secondfile = $opt_j;
} else {
  $secondfile = $opt_s;
}

$bl2seq_output = "bl2seq.out";
if ($opt_o) {
  $bl2seq_output = $opt_o;
}

$bl2seq_program = "tblastn";
if ($opt_p) {
  $bl2seq_program = $opt_p;
}

$MIN_EXPECTATION = "1e-30";
if ($opt_e){
  $MIN_EXPECTATION = $opt_e;
}



print "1st fasta: $firstfile\n";
print "2nd fasta: $secondfile\n";
print "program: $bl2seq_program\n";
print "expectation: $MIN_EXPECTATION\n";



$original = Bio::SeqIO->new(-file => $firstfile,
			    '-format' => 'fasta');


# bl2seq each of the contigs against the original

$i = 1;				#original protein iteration number;
while (defined($orig = $original->next_seq())) {

  $j = 1;

  $contigs = Bio::SeqIO->new(-file=> $secondfile,
			     '-format' => 'fasta');

  #rebuilt sequence iteration number;
  while (defined($curContig = $contigs->next_seq())) {
    $cur_file = $bl2seq_output.'.'.$i.'.'.$j.'-'.$orig->display_id().'.'.
      $curContig->display_id();

    # Make sure there is no naughty characters in the filename
    $cur_file =~ s/[\|\*\[\]\{\}]//g;

    print "writing bl2seq output file: $cur_file...";
    $factory = Bio::Tools::Run::StandAloneBlast->new('program' => 
						     $bl2seq_program,
						     -outfile => $cur_file,
						     'F' => 0, # no filtering low complexity regions
						     'e' => $MIN_EXPECTATION
						    );


    $o = length($orig->seq());
    $c = length($curContig->seq());
    if ($o==0 || $c==0) {
      print STDERR "zero length sequence somewhere below:\n";
      print STDERR "orig: $o\n";
      print STDERR "curContig: $c\n\n";
      exit(1);
    }


    #$bl2seq_report = 
    $factory->bl2seq($orig,$curContig);
    print "done";


    # Use AlignIO.pm to create a SimpleAlign object from the bl2seq report
    $str = Bio::AlignIO->new(-file=> $cur_file,'-format' => 'bl2seq');

    # delete the output, whatever happens
    unlink $cur_file;

    # Get out the first align if possible
    $first = $str->next_aln();
    if (defined($first)){
      print " hit found";


      # Re-run the bl2seq, this time including all hits, not
      # just the ones above MIN_EXPECTATION
      $factory2 = Bio::Tools::Run::StandAloneBlast->new('program' => 
							$bl2seq_program,
							-outfile => $cur_file,
							'F' => 0, # no filtering low complexity regions
						       );

      $factory2->bl2seq($orig,$curContig);


    } else {
      # delete the output
      print " no hit";
    }


    $j++;
    print "\n";
  }

  $i++;
}

exit();







sub usage()
  {
    print "Usage parameters:\n";
    print " -f or -i (required) first input file\n";
    print " -s or -j (required) second input file\n";
    print " -p bl2seq blast program (eg. tblastn). Defaults to 'tblastn'.\n";
    print " -o bl2seq direct output file. Defaults to 'bl2seq.out'.\n";
    print " -e min expectation value for a the best hit in a file.\n";
    exit(0);
  }
