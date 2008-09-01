#!/usr/bin/perl -w

#Author: Benjamin Woodcroft: donttrustben at gmail DO.T com
#Last modified 16 Mar 2008

# This script takes in a XML output from blast (-m7 using blastcl3)
# The output a list of best hits, that can be opened nominally
# as a CSV file in excel, and then they can be added to the
# microarray excel document easily.




use Bio::SearchIO;
use Bio::SearchIO::Writer::HTMLResultWriter;
use Bio::Assembly::IO;

use Getopt::Std;


our ($opt_s);
getopts('s');




if ($#ARGV != 1) {
  print STDERR "usage: $0 [<options>] <blastxml_file> <ACE_file> \n".
      "Options:\n".
      "\t-s: contigs number found from name 'ContigX', otherwise singlet is\n".
      "\t    assumed. Default is ContigX and SingletX, otherwise fail.\n\n";
  exit;
}

my $xmlfile = $ARGV[0];


my $searchio = Bio::SearchIO->new(-file   => $xmlfile,
				  -format => 'blastxml') or die "parse failed";



my $SEP = "\t";




# Print the Column headings
my @headings = ('Name', 'Cluster', 'Autoblast Webpage', 'Autoblast E-value', 'Autoblast Accession', 'Autoblast Name', 'Autoblast Description'); #names of the output data columns
my $line = join $SEP, @headings;
print $line."\n";


#my %seq_to_contig = &create_seq_to_contig_hash($ARGV[1]);
#testing
#print 'final - '.$seq_to_contig{'CHR4'}."\n";


my $io = new Bio::Assembly::IO(-file=>$ARGV[1],
			       -format=>"ace");

my $assembly = $io->next_assembly;

#    my @seq_ids = $assembly->get_contig_by_id($contig)->get_seq_ids;

# Loop over the blast hits for each contig.
# Record information about the top hit;
# name
# accession
# e-value
while (my $result = $searchio->next_result()) {

  my $contig = $result->query_name;

  my @seqs;


  # Classify the input sequence as 
  if (!($contig =~ m/^Contig(\d+)$/)) {

      # opt_s means whole thing is singlet if not a contig
      if ($opt_s){
	  my $singlet_id = $contig;
	  @seqs = ($singlet_id);
      } else {

	  # Improper input
	  if (!($contig =~ m/^Singlet(\S+)$/)) {
	      die "Query name not 'ContigX' or SingletY like. This script can't handle anything else. It would need to be rewritten.";
	  }

	  # Singlet input
	  else {
	      my $singlet_id = $1;
	      @seqs = ($singlet_id);
	      #print STDERR "processing singlet $singlet_id: ".($#seqs+1)."\n";
	  }
      }
  }

  # Contig Input
  else {
    my $contig_id = $1;
    @seqs = $assembly->get_contig_by_id($contig_id)->get_seq_ids;
    #print STDERR "processing contig $contig: ".($#seqs+1)."\n";
  }

  my $autoblast_link = 'http://reefedge.sols.uq.edu.au/autoblast/AbaloneMicroarray1/blast_xml_hits/'.$contig.'.html';

  my $hit = $result->next_hit;
  if (defined($hit)) {
    my $hsp = $hit->next_hsp;

    #print sequence name, contig id, blast stuff
    # blast stuff is e-value accession name
    my $score = $hsp->expect;
    my $acc = $hit->accession;
    my $name = $hit->name;
    my $description = $hit->description;
    

    my $rest = $contig."$SEP".
      "$autoblast_link$SEP$score$SEP$acc$SEP$name$SEP$description";

    foreach my $s (@seqs) {
      print $s."$SEP$rest\n";
    }
  }

  # if there was no top hit at all, just print the contig name
  else {
    #print "not using seqs\n";
    foreach my $s (@seqs) {
      print $s."$SEP".
	$contig."$SEP".
	  $autoblast_link."$SEP\n";
    }
  }

  #die "reached the loop end";
}








# Create a hash that maps sequence names to contig ids
# in an assembly
# This is DEAD CODE I think.
sub create_seq_to_contig_hash {
  my $assembly_filename = $_[0];

  # Assembly loading methods
  my $io = new Bio::Assembly::IO(-file=>$assembly_filename,
				 -format=>"ace");

  my $assembly = $io->next_assembly;
  my @contig_ids = $assembly->get_contig_ids;


  foreach my $contig (@contig_ids) {
    my @seq_ids = $assembly->get_contig_by_id($contig)->get_seq_ids;
    foreach my $seq_id (@seq_ids) {
      if ($seq_to_contig{$seq_id}) {
	print STDERR "Warning: Sequence $seq_id occured multiple times in the assembly - this script will probably fail";
      } else {
	#print $seq_id.' - '.$contig."\n";
	$seq_to_contig{$seq_id} = $contig;
      }
    }
  }

  return %seq_to_contig;
}
