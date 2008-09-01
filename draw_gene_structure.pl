#!/usr/bin/perl -w

# Create the panel and track corresponding to various gene structure elements.

use strict;


use Bio::Graphics;
use Bio::SeqFeature::Generic;
my $bsg = 'Bio::SeqFeature::Generic';

my $span         = $bsg->new(-start=>1,-end=>1000);
my $span2         = $bsg->new(-start=>1000,-end=>2000);
my $test_feature = $bsg->new(-start=>300,-end=>700,
                             -display_name=>'Test Feature',
                             -source_tag=>'This is only a test');
my $test_feature2 = $bsg->new(-start=>700,-end=>1000,
                             -display_name=>'Test Feature',
                             -source_tag=>'This is only a test2');
my $test_feature3 = $bsg->new(-start=>1000,-end=>1300,
                             -display_name=>'Test Feature',
                             -source_tag=>'This is only a test2');
 
my $panel        = Bio::Graphics::Panel->new(-width=>600,-length=>$span->length,
                                             -pad_left=>12,-pad_right=>12);
$panel->add_track($span,-glyph=>'arrow',-double=>1,-tick=>2);
 
$panel->add_track($test_feature,
                  -glyph   => 'generic',
                  #-bgcolor => 'orange',
                  -font2color => 'red',
                  -height  => 20,
                  -label   => 1,
                  -description => 1,
   );
$panel->add_track($test_feature2,
                  -glyph   => 'generic',
                  -bgcolor => 'black',
                  -font2color => 'red',
                  -height  => 2,
                  -label   => 1,
                  -description => 1,
   );
$panel->add_track($span2,-glyph=>'arrow',-double=>1,-tick=>2);
$panel->add_track($test_feature3,
                  -glyph   => 'generic',
                  -bgcolor => 'orange',
                  -font2color => 'red',
                  -height  => 20,
                  -label   => 1,
                  -description => 1,
   );
 
print $panel->png;
