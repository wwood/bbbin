#!/usr/bin/perl -w

#Author: Benjamin Woodcroft: donttrustben at gmail DO.T com
#Last modified 26 Sept 2007

# This script takes in a XML output from blast (-m7 using blastcl3)
# and converts the first non-hypothetical hit's hsps
# into a exonpart hint that can be used by exonerate



use Bio::SearchIO;
use Getopt::Std;


# use getopts not getopt because you want -g to be boolean
our ($opt_g,$opt_r);
getopts('gr:');




if ($#ARGV != 0) {
  print STDERR "usage: $0 [-g] <blastxml_file>\n";
  print STDERR "\tg - make output more useful for gff2ps\n";
  exit;
}

my $xmlfile = $ARGV[0];


my $searchio = Bio::SearchIO->new(-file   => $xmlfile,
				  -format => 'blastxml') or die "parse failed\n";





if (my $result = $searchio->next_result()) {




  # get the name of the sequence and the biggest hit score

  my $hit = $result->next_hit();
  my $query_name = $result->query_name() or die "unknown hit sequence name\n";

  # remove junk from the filename if wanted. Bit of extraneous feature though.
  if ($opt_r){
      $query_name =~ s/$opt_r//;
  }



  if (!defined($hit)){
    print STDERR "no hits found";
    exit;
  }


  #keep searching until a decent hit is found
  my $name;
  do {
    $name = $hit->name() or die "could not parse out hit name\n";
    $desc = $hit->description or die "could not parse hit description\n";

    # Ignore genes that are hypothetical
    foreach $word ('hypothetical','unknown','predicted','like'){
	if ($desc =~ m/$word/i){
	    print STDERR "Found description with $word, ignoring: $desc ($name)\n";
	    $hit = undef;
	}
    }
  } while (!defined($hit) and $hit = $result->next_hit());


  if (!defined($hit)){
    print STDERR "some hits found, but all were unsuitable, so ignoring.\n";
    exit;
  }



  while (my $hsp = $hit->next_hsp()){
    my $lseq = $hsp->get_aln->get_seq_by_pos(1);
    my $strand = '.';
    if ($hsp->strand eq '1'){$strand = '+';}
    elsif ($hsp->strand eq '-1'){$strand = '-';}
    else {print STDERR "Strand Unknown in hit $name - ".$hsp->strand."\n";}

    if ($opt_g){
	print 'blast'."\t$query_name\texonpart\t".$lseq->start."\t".$lseq->end."\t10\t".$strand."\t.\tsource=M ; description: ".$hit->description."\n";
    } else {
	print $query_name."\tanchor\texonpart\t".$lseq->start."\t".$lseq->end."\t0\t".$strand."\t.\tsource=M ; description: ".$hit->description."\n";
    }
  }

}

