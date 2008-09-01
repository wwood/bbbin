#!/usr/bin/perl -w
 
# This is similar to code example 4 in the Graphics-HOWTO
use strict;
#use lib "$ENV{HOME}/projects/bioperl-live";
use Bio::Graphics;
use Bio::SearchIO;
use Bio::Graphics::Glyph;
 
my $file = shift or die "Usage: $0 <blast file>\n";
 
my $searchio = Bio::SearchIO->new(-file   => $file,
                                  -format => 'blastxml') or die "parse failed";
 
my $result_index = 1;

while (my $result = $searchio->next_result()){
 
my $panel = Bio::Graphics::Panel->new(
                                      -length    => $result->query_length,
                                      -width     => 800,
                                      -pad_left  => 10, #sometimes the name goes for longer
                                      -pad_right => 100,
                                     );
 
my $full_length = Bio::SeqFeature::Generic->new(
                                                -start        => 1,
                                                -end          => $result->query_length,
                                                -display_name => $result->query_name,
                                               );
$panel->add_track($full_length,
                  -glyph   => 'arrow',
                  -tick    => 2,
                  -fgcolor => 'black',
                  -double  => 1,
                  -label   => 1,
                 );
 
my $track = $panel->add_track(
                              -glyph       => 'graded_segments',
			      -height      => 5,
                              -label       => 1,
			      #-linewidth   => 0,
				-min_score => 0,
				-max_score => 250,
                              -connector   => 'dashed',
                              -bgcolor     => 'red',
                              -font2color  => 'blue',
                              -sort_order  => 'high_score',
                              -description => sub {
                                my $feature = shift;
                                return unless $feature->has_tag('description');
                                my ($description) = $feature->each_tag_value('description');
                                my $score = $feature->score;
                                "$description, score=$score";
                               },
                             );

while( my $hit = $result->next_hit ) {
  #next unless $hit->significance < 1E-20;

  # The score here is not a simple thing. There does not appear to be any 'overall' score, but the individual parts of the hit
  # are scored differently. I'm going to use the hit score from the first hit, if this is not available.

  #print $hit->num_hsps."\n";

  my @hsps = $hit->hsps;
  my $first_hsp = shift @hsps;
  my $feature = Bio::SeqFeature::Generic->new(
                                              -score        => $first_hsp->score,
                                              #-display_name => $hit->description,
                                              -tag          => {
								#significance =>$first_hsp->evalue,
								description =>$hit->description
                                                               }
                                             );
#  my $feature = Bio::Graphics::Glyph::graded_segments->new(
#                                              -score        => $first_hsp->evalue,
#                                              -display_name => $hit->name,
#                                              -tag          => {
#                                                                description => $hit->description
#                                                               },
#						-max_score  => 40,
#						-min_score  => 0
#                                             );
#  while( my $hsp = $hit->next_hsp ) {
##    $feature->add_sub_SeqFeature($hsp,'EXPAND');
#    $feature->add_SeqFeature($hsp,'EXPAND');
#  }
  $feature->add_SeqFeature($first_hsp,'EXPAND');
 
  $track->add_feature($feature);
}
 


#Open a file with the graphics panel in it
my $png_name = "result$result_index.png";
open OUT, ">$png_name" or die "Couldn't open file for writing png: $png_name";
print OUT $panel->png;
close OUT;

# Reset for next time
$result_index++;
}
