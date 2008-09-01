#!/usr/bin/perl -w

#Author: Benjamin Woodcroft: donttrustben at gmail DO.T com
#Last modified 23 June 2008

# This script takes in a XML output from blast (-m7 using blastcl3)
# The output is a series of HTML files, 1 for each query sequence.
# It also outputs a tableofcontents.html file that lists all the sequences and
# hyperlinks to them.




use Bio::SearchIO;
use Bio::SearchIO::Writer::HSPTableWriter;





if ($#ARGV != 0) {
  print STDERR "usage: $0 <blastxml_file>\n";
  exit;
}

my $xmlfile = $ARGV[0];


my $searchio = Bio::SearchIO->new(-file   => $xmlfile,
				  -format => 'blastxml') or die "parse failed";



# blast -m8 Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
my @array = qw(
	       query_name
	       hit_name
	       frac_identical_query
	       length_aln_query
	       gaps_total
	       start_query
	       end_query
	       start_hit
	       end_hit
	       bits
	       hit_description
	      );

# print headings
my $str = join "\t",@array;
print $str."\t";

# print results
my $writer = Bio::SearchIO::Writer::HSPTableWriter->new(-columns => [@array]);

my $out = Bio::SearchIO->new( -writer => $writer,
                                  -fh   => \*STDOUT );

while ( my $result = $searchio->next_result() ) {
  $out->write_result($result);
}
