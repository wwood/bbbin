#!/usr/bin/perl -w

#Author: Benjamin Woodcroft: donttrustben at gmail DO.T com
#Last modified 26 Sept 2007

# This script takes in a XML output from blast (-m7 using blastcl3)




use Bio::SearchIO;
use Bio::SearchIO::Writer::HTMLResultWriter;





if ($#ARGV != 0) {
  print STDERR "usage: $0 <blastxml_file>\n";
  exit;
}

my $xmlfile = $ARGV[0];


my $searchio = Bio::SearchIO->new(-file   => $xmlfile,
				  -format => 'blastxml') or die "parse failed";

my $first = 1;

my $DIR = 'blast_xml_hits';	#name of the directory where the individual hit files go
my $i = 1;


my @headings = ('Name', 'Cluster', 'Top Autoblast Hit'); #names of the output data columns
#print join ',', @headings;
#exit;



while (my $result = $searchio->next_result()) {

  #if this is the first time, create the index and the new directory
  #If this was done outside the loop, then if no results were found then
  #The directories would still be created - so this is not something I
  #want to do.
  if ($first) {
    #create new directory where the individual hits are shown
    mkdir "$DIR" or die "failed to create hits directory: $!";
    open INDEX, '>index.html' or die 'failed to open contents file';
    print INDEX "<html><body><table>\n";
  }
  $first = 0;
  $result->rewind;


  # write to the index
  # get the name of the sequence and the biggest hit score
  my $name = $result->query_name() or die 'unknown name';
  my $hitfile = "$DIR/$name.html";
  my $hit = $result->next_hit();
  my $score = '-';
  my $hsp;
  if (defined($hit)){$score = $hit->next_hsp()->expect(); $hit->rewind;}
  print INDEX "<tr><td><a href='$hitfile'>$name</a><td>$score</tr>\n";



  # next lines stolen from http://www.bioperl.org/wiki/HOWTO:SearchIO
  my $writerhtml = new Bio::SearchIO::Writer::HTMLResultWriter();
  my $outhtml = new Bio::SearchIO(-writer => $writerhtml,
				  -file   => ">$hitfile");
  # get a result from Bio::SearchIO parsing or build it up in memory
  $outhtml->write_result($result);



  $i += 1;
}

#close the index file if one was created
if ($first==0) {
  print INDEX "</table></body></html>";
  close INDEX;
}

exit;
#!/usr/bin/perl -w

use Bio::SearchIO;
print 1;flush STDOUT;
my $in = new Bio::SearchIO(-file   => $ARGV[0],
			   -format => 'blastxml');
print 2;flush STDOUT;
my $out = new Bio::SearchIO(-output_format  => 'GbrowseGFF',
			    -output_cigar   => 1,
			    -output_signif  => 1,
			    -file           => ">result.gff");
print 3;flush STDOUT;
while ( my $r = $in->next_result ) {
  #print 'yo';
  $out->write_result($r);
}
