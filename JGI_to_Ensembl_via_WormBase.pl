#!/usr/bin/perl -w

# An attempt at using the wormbase scripts from Ensembl to
# convert them into mysql Ensembl tables from the JGI Genes GFF
# file.


use lib '/home/uyen/bioinfo/ensembl/ensembl-pipeline/ensembl-pipeline/scripts/DataConversion/wormbase';
use WormBase;

use strict;
use Getopt::Std;

use Bio::EnsEMBL::Analysis;

our $opt_g;

getopt('g');
my $gff;
if (defined($opt_g)) {
  $gff = $opt_g;
} else {
  &usage();
  exit;
}

# Connect to the database
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
					      -user   => 'guest',
					      -dbname => 'ensembl_nematostella',
					      -host   => 'localhost',
					      -driver => 'mysql'
					     );

my $slice_adaptor = $dba->get_SliceAdaptor();



# Create an analysis object
my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => 'JGI');
my $analysis_adaptor = $dba->get_AnalysisAdaptor();
$analysis_adaptor->store($analysis);

# Get the sequence
my $seq = $slice_adaptor->fetch_by_region('scaffolds','scaffold_258') or die "couldn't get sequence";
&parse_gff($gff, $seq, $analysis);





sub usage {
  print "Usage: $0 -g <JGI_GFF>\n";
}
