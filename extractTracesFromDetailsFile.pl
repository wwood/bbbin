#!/usr/bin/perl -w

# The purpose is extract the trace ID's from the details file that comes with the rebuilding step

use Getopt::Std;


$usage = "usage: extractTracesFromDetails.pl <detailsFile> <options>\n".
  "OPTIONS:\n".
  " -e EXCEPT - get all contigs except the ones given in the comma separated list argument following this option\n".
  " -s SPECIC - only get the traces from the contigs given in the comma separated list following this argument\n".
  "only traces are printed to STDOUT, metadata goes to STDERR\n\n";
$details = shift or die "$usage";

#What are we doing here?
getopt('es');

$excepting = 0;
$specificing = 0;

if ($opt_e){
  $excepting = 1;
  @contigList = split ',', $opt_e;
}
elsif ($opt_s){
  $specificing = 1;
  @contigList = split ',', $opt_s;
}
else {die "$usage"};


open DETAIL, "$details" or die "failed to open \"$details\"";
@dets = <DETAIL>;

# for each contig, extract the traces (Just from the top section)

foreach $i (0..$#dets) {
  $line = $dets[$i];

  # if end of top section
  if ($line =~ m/DETAILED DISPLAY OF CONTIGS/) {
    last;
  }

  # if start of new contig
  #print "lineRaw: $line";
  if ($line =~ m/\*\*\*\*\*\*\*\*\*\*\*\*\* (.*) \*\*\*\*\*\*\*/) {
    #print "gotcha\n";
    $contigName = $1;

    $i++;
    @traces = ();
    while (!(($line = $dets[$i]) =~ m/\*\*\*\*\*\*\*\*/)) {

      #print "newline: $line";
      if ($line =~ m/^(\S+)\s+\d+/) {
	$trace = $1;
	#print "trace1: $trace\n";

      } else {
	
	if ($line =~ m/^\s+(\S+) is in (\S+)/) {
	  $trace = $1;
	  #print "trace2: $trace\n";
	} elsif ($line =~ m/^\s*\n$/ ||
		 $line =~ m/^\s*The contigs above and below are linked by 2 constraints.\s*\n$/ ||
		 $line =~ m/^\s*DETAILED DISPLAY OF CONTIGS\s*\n$/ ||
		 $line =~ m/^\s*The contigs above and below are linked by 6 constraints.\s*$/
		) {
	  ;			# ignore this line
	} else {
	  print STDERR "unrecognized line: $line";
	  $i++;
	  next;
	}
      }

      push @traces, $trace;


      $i++;
    }
    $i--;


    #print the new contig if:
    # -we are specificing and the contig is in the contig list -or-
    # -we are excepting and the contig is not in the contig list
    @matches = grep m/$contigName/, @contigList;
    if ($#matches != -1){$matched = 1}
    else {$matched = 0;}
    #print ">>>>Contig: $contigName ^^$matches[0]^^ $#matches $matched\n";
    if (($specificing == 1 && $matched==1) ||
	($excepting == 1 && $matched!=1)) {

      print STDERR ">>>>>>>>>$contigName\n";
      foreach (@traces) {
	s/\+//g;
	s/-//g;
	print "$_\n";
      }
    }



  }
}

