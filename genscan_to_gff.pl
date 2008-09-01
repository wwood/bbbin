#!/usr/bin/perl -w


#Convert a genscan file into a GFF file, usiong bioperl libraries.
#does not actually run the genscan, just uses the results.


use Bio::Tools::Genscan;
use Bio::FeatureIO::gff;
use Bio::SeqFeature::Annotated;
use Getopt::Std;


getopt('n');
my $underlying_name = 'SEQ';
if ($opt_n){
  $underlying_name = $opt_n;
}


#load the genscan
$genscan = '';
if ($ARGV[0]){
  $genscan = Bio::Tools::Genscan->new(-file => $ARGV[0]);
} else {
  $genscan = Bio::Tools::Genscan->new(-fh => \*STDIN);
}
my $featureOut = Bio::FeatureIO->new(
				     -format => 'gff',
				     -version => 3,
				     -fh => \*STDOUT,
				     -validate_terms => 0
				    );


#@names = ('fgf',
#'hex',
#'hypothetical1',
#'hypothetical2',
#'hypothetical3',
#'hypothetical4',
#'inositol',
#'kenesin1',
#'kenesin2',
#'major_facilitator',
#'msx',
#'nk2',
#'nk567a',
#'nk567b',
#'nk567c',
#'rolling_stone',
#'tetracycline',
#'tlx',
#'vacuolar');
#$name_index = 0;




#write the gff
#Can't seem to get it to convert usefully. Annoying. Have to do it manually.
#<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
$i=1;
while($gene = $genscan->next_prediction()) {
  if ($opt_n){
    print "###################### $underlying_name #################\n";
  }

  #$f = Bio::SeqFeature::Annotated->new();
  @feats = $gene->features;
  foreach $f (@feats){
    $gff = $f->gff_string;

    $annotation = $underlying_name.'Gene'.$i;

    $gff =~ s/SEQ/$underlying_name/;
    $gff =~ s/(.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t).*/$1$annotation/;
    print $gff."\n";
  }
  $i++;
}
