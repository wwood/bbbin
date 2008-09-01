#!/usr/bin/perl

use lib "/home/leishcyc/lib/perl5/site_perl/5.8.5/";
use Bio::Perl;

use strict;
use Bio::SearchIO;
use Getopt::Std;

our ($opt_s, $opt_p);
getopt('sp');

my $usage = "usage: blast_parse_aw.pl [-s output_separator] [-p max_number_of_hsps_per_hit] <blast report> <max # of hits>\n";
if (@ARGV != 2) { die $usage; } 

my $infile  = $ARGV[0];
my $maxhits  = $ARGV[1]; # to specify in the command line how many hits
                        # to report for each query

# Set the default output separator
my $separator = "\t";
if ($opt_s){
 $separator = $opt_s;
}

# Set the number of HSPs found
my $maxhsps = 1000;
if ($opt_p){
 $maxhsps = int $opt_p;
}

open(  IN, $infile )     || die "Can't open inputfile '$infile'! $!\n";


my $report = new Bio::SearchIO(-file => $infile);

my $titles = join($separator,
  ('QueryName',
'QueryDesc',
'Len',
'HitName',
'Desc',
'Len',
'Num_Hsps',
'Start_Query',
'Start_Hit',
'PerIden',
'HspLen',
"Evalue\n"));
print $titles;

while( my $result = $report->next_result ) {


  
my $i = 0;

while( ($i < $maxhits ) && (my $hit = $result->next_hit ) ) {

$i++;

my $num_hsps = 0;

while ( ($num_hsps < $maxhsps) && ( my $hsp = $hit->next_hsp ) ) {
$num_hsps++;
           print $result->query_name,        $separator;
print "\"";
           print $result->query_description, "\"", $separator;
           print $result->query_length,      $separator;
           print $hit->name,                 $separator;
print "\"";
           print $hit->description, "\"",    $separator;
           print $hit->length,               $separator;
           print $hit->num_hsps,             $separator;
           print $hsp->start('query'),       $separator;
           print $hsp->start('hit'),         $separator;
           print $hsp->percent_identity,     $separator;
           print $hsp->length('total'),      $separator;
           print $hsp->evalue,               "\n"; 
 }
}

if ($i == 0) { 
print $result->query_name,        $separator;
print "no hits\n"; 
}

}
