#!/usr/bin/perl -w

# JGI releases data in GFF format for the genes. The idea of this script
# is to take that GFF data and insert it into a MySQL Enxembl database.
# Originally, this was because the pseudogene.org pipeline
# required this.

#use strict;
use Bio::Tools::GFF;
use Bio::SeqIO;

use Getopt::Std;


# Check the inputs are correct
#my $opt_g;
getopt('gf');
my $gff;
my $seqio;
if (defined($opt_g)) {
  $gff = Bio::Tools::GFF->new(-file => $opt_g);
} else {
  &usage();
  exit;
}
if (defined($opt_f)) {
  $seqio = Bio::SeqIO->new(-file => $opt_f, -format=>'Fasta');
} else {
  &usage();
  exit;
}







# I can print the exons out in any order i want - just need to get the
# ranks correct in the exon_transcript table.
# BUT I do need to know the first and last exon, because these are in
# the Translation table.

# Exons do not seem to have a phase in the JGI script. What does this mean?

my $last_gene_name = undef;


while (my $f = $gff->next_feature) {

  my $cur_gene_name = ($f->annotation()->get_Annotations('name'))[0] or die "gene name not found";

  # first line of first gene only.
  if (!defined($last_gene_name)) {
    $last_gene_name = $cur_gene_name;
  }

  # if it is a change in genes
  if (!($cur_gene_name eq $last_gene_name)) {
#    print "found gene '$last_gene_name' with ".
#      ($#exons+1).' exons, '.
#	($#cds+1).' CDS,'.
#	  ($#start_codons+1).' start codons, and '.
#	    ($#stop_codons+1)." stop codons.\n";


    # Print out each of the exon hits, because they are the simplest.


    # Print out the experiment. That is, add
    # Another experiement. If you take the sum total of all the bases in a protein, then is this number divisible by 3? Or is the end like the start, having random bases at the end? I can't see why not, but doing the experiment will probably help me get the script together anyway, so I'll do that.
    my $total = 0;
    $total += $cds[0]->end - $cds[0]->start - $cds[0]->frame+1;
    foreach my $i (1..$#cds) {
      print $cds[$i]->frame."\t".($total%3)."\n";
      $total += $cds[$i]->end - $cds[$i]->start+1;
#      if ($cds[$i]->frame == (3-($total%3))%3){
#	print "yes\n";
#      } else {
#	print "no\n";
#      }
#      $total += $cds[$i]->end - $cds[$i]->start+1;
    }
    #print "frame: ".($total%3)." $cur_gene_name\n";



    #offset should be 3-(total%3)


    #exit;

    # reset the variables
    $last_gene_name = $cur_gene_name;
    @exons = ();
    @cds = ();
    @start_codons = ();
    @stop_codons = ();
  }

  # Else we are just using the current gene.
  elsif ($f->primary_tag eq 'exon') {
    push @exons, $f;
  } elsif ($f->primary_tag eq 'CDS') {
    push @cds, $f;
  } elsif ($f->primary_tag eq 'start_codon') {
    push @start_codons, $f;
  } elsif ($f->primary_tag eq 'stop_codon') {
    push @stop_codons, $f;
  } else {
    print STDERR "Unknown primary tag ".$f->primary_tag." found.\n";
  }
}

print "HERE";
exit;
{

  #print $cur_gene_name."\n";
  exit;


  my @anns = $f->annotation->get_all_annotation_keys;
  print $#anns."\n";
  #print $anns."\n";
  exit;

  if ($f->primary_tag eq 'CDS') {
    #print $f->frame."\n";
    #print $f->start."\t".(($f->end - $f->start) % 3)."\n";
  }
  #print $f->start()."\n";
  #print $f->end()."\n";
}







exit;


# Connect to the database
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
					      -user   => 'guest',
					      -dbname => 'ensembl_nematostella',
					      -host   => 'localhost',
					      -driver => 'mysql'
					     );



#my $f = $gff->next_feature;
#print $f->seq_id()."\n";
#print $f->start()."\n";
#print $f->end()."\n";
#print $f->primary_tag()."\n";
#print $f->strand()."\n";
#my $gs = ($f->annotation()->get_Annotations('name'))[0];
##my $gs = join '--',@fs;
#print $gs."\n";
#exit;


#&fasta_to_seq_region($seqio, $dba);
&gff_to_tables($gff, $dba);




# Convert the JGI GFF into a Ensembl Tables
# Assumes the seq_regions are already defined
# Usage gff_to_tables(gff, dba)
sub gff_to_tables {

  my $gff = shift or die "couldn't read gff from arguments";
  my $dba = shift or die "couldn't read dba from arguments";

  my $exon_adaptor = $dba->get_ExonAdaptor();
  my $slice_adaptor = $dba->get_SliceAdaptor();
  my $translation_adaptor = $dba->get_TranslationAdaptor();


  #my $f = $gff->next_feature;
  #print $f->()."\n";
  #print $f->start()."\n";
  #print $f->end()."\n";
  #print $f->primary_tag()."\n";
  #print $f->strand()."\n";
  #my $gs = ($f->annotation()->get_Annotations('name'))[0];
  ##my $gs = join '--',@fs;
  #print $gs."\n";
  #exit;

  # Create an analysis object (Needed to create a Gene)
  my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => 'JGI');
  my $analysis_adaptor = $dba->get_AnalysisAdaptor();
  $analysis_adaptor->store($analysis);


  # Keep track of what the current gene is
  my @cur_exons;
  my @cur_cds;
  # Keep track of waht the last start thing was
  my $cur_gene_name;
  my $cur_transcript_name;

  my $cur_strand;
  my $cur_gene_slice;


  # for each GFF feature
  while ($f = $gff->next_feature()) {

    # if it is an exon a new gene, the last transcript is
    #ready to be processed.
    my $gene_name = ($f->annotation()->get_Annotations('name'))[0] or die "couldn't find name";

    #if it is the first go, do nothing
    if (!defined($cur_gene_name)) {
      $cur_gene_name = $gene_name;


      # Else if it is the beginning of a new gene, the old gene needs to
      # be uploaded.
    } elsif (!($gene_name eq $cur_gene_name)) {
      $cur_gene_name = $gene_name;

      # process a whole gene
      &create_new_gene($dba, $cur_strand, $cur_gene_slice, \@cur_exons, $first_cds_exon, $last_known_cds_exon, $analysis);	

      # Reset all the variables
      $first_cds_exon = undef;
      $last_known_cds_exon = undef;
      @cur_exons = undef;
    }


    #what type of thing are we parsing
    if ($f->primary_tag eq 'exon') {


      # create exon
      my $slice = $slice_adaptor->fetch_by_region(undef, $f->seq_id()) or die "could not find the slice for ".$f->seq_id.".";
      $cur_gene_slice = $slice;
      $cur_strand = $f->strand;
      my $ex = new Bio::EnsEMBL::Exon(-START     => $f->start(),
				      -END       => $f->end(),
				      -STRAND    => 1, #$cur_strand,
				      -SLICE     => $slice,
				      -VERSION   => 1,
				      -PHASE     => 0,
				      -END_PHASE => 1 #TODO: make these right.
				     );

      # store the exon in the list of current gene exons
      $exon_adaptor->store($ex);
      $last_exon = $ex;
      push @cur_exons, $ex;
    } elsif ($f->primary_tag eq 'CDS') {




      # if it is a CDS, create a new transcript and add
      # it to the list.
      if (defined($first_cds_exon)) {
	$first_cds_exon = $last_exon;
      }
      $last_known_cds_exon = $last_exon;
    }
  }



  # process the last gene on the list.

  print "Exons length: $#cur_exons\n";
  print "Slice: $cur_gene_slice\n";
  &create_new_gene($dba, $cur_strand, $cur_gene_slice, \@cur_exons, $first_cds_exon, $last_known_cds_exon, $analysis);
}







# A subroutine to put a whole gene into the database
sub create_new_gene {
  my ($dba, $cur_strand, $cur_gene_slice, $cur_exons, $first_cds_exon, $last_known_cds_exon, $analysis) = @_ or die "bad usage of create_new_gene";

  # Get adaptors
  my $gene_adaptor = $dba->get_GeneAdaptor();
  my $transcript_adaptor = $dba->get_TranscriptAdaptor();

  # create a gene and a transcript (taken to be the same thing)
  # All exons are included here
  my $gene = Bio::EnsEMBL::Gene->new(-START => $$cur_exons[0]->start,
				     -END   => $$cur_exons[$#$cur_exons]->end,
				     -STRAND=> $cur_strand,
				     -SLICE => $cur_gene_slice,
				     -ANALYSIS => $analysis);
  $gene_adaptor->store($gene) or die "couldn't store gene..";
  my $transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $cur_exons);
  $transcript_adaptor->store($transcript, 1, 1); #don't care that the gene table isn't complete.

  # create a translation
  # Only exons that are part of CDS's are included here.
  my $translation = Bio::EnsEMBL::Translation(-START_EXON => $first_cds_exon,
					      -END_EXON => $last_known_cds_exon);
}



sub gff_strand_to_ensembl {
  my $gff_strand = shift or die "need a gff strand to start";
  print $gff_strand."\n";
  exit;
}




# Go through the fasta sequences, creating new seq_region entry for each one.
# Usage fasta_to_seq_region(seqio, dba)
sub fasta_to_seq_region {
  my $seqio = shift or die "need a valid seq io file";
  my $dba = shift or die "need a valid dba second";

  my $cs_adaptor = $dba->get_CoordSystemAdaptor() or die "could not create adaptor";
  my $slice_adaptor = $dba->get_SliceAdaptor() or die "could not create adaptor";

  # I'm debugging so this works, but when you create one for the first time
  # you'll need to create your own scaffold system.
  my $cs = $cs_adaptor->fetch_by_dbID(1) or die "no coord system defined.";
  #Bio::EnsEMBL::CoordSystem->new(-NAME    => 'scaffolds',
  #             -VERSION => 'JGI-ben2',
  #             -RANK    => 1);
  #$cs_adaptor->store($cs);

  while ($s = $seqio->next_seq) {

    #print $s->display_id()."\n";
    #print $s->length()."\n";
    my $len = $s->length();
    my $region = new Bio::EnsEMBL::Slice(-coord_system => $cs,
					 -start => 1,
					 -end => $len,
					 -strand => 1,
					 -seq_region_name => $s->display_id,
					 -seq_region_length => $len);
    $slice_adaptor->store($region) or die "failed to store slice.";
  }
}





sub usage {
  print "Usage: $0 -f <JGI_Fasta> -g <JGI_GFF>\n";
}




