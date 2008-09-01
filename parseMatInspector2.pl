#!/usr/bin/perl -w

# THis is an attempt to not have to copy out from the HTML, because that takes ages.

use  promoterAnalysis;


$given_seq_name = $ARGV[1];


open IN, "<$ARGV[0]" or die "HTML file opening failed: $ARGV[0]";

@everything = <IN>;
$everythingString = join "", @everything;

@bigParts = split "<h2>Inspecting sequence ", $everythingString;
#print $#bigParts.">>>";

foreach $part (@bigParts){

  # Find the sequence name
  $part =~ m/\[(.*)\]/;
  $seq_name = $1;

  #print "Test one: $seq_name>$given_seq_name\n";
  if (!$seq_name || !($seq_name eq $given_seq_name)){next;}
  #print "Found one: $seq_name\n";


  @pieces = split "<tr align=\"center\">", $part;
  #print $#pieces.">>>";

  #skip the first one because it is not a real value
  foreach $i (1..$#pieces){
    # Read in a single putative site
    $pieces[$i] =~ s/\n\n/\n/;
    #print "\n\n$pieces[$i]\n";
    @site_lines = split "\n", $pieces[$i];


    #EXAMPLE TABLE ENTRY
    #<tr align="center">
    #<td bgcolor="#7F5BE9">&nbsp;</td>
    #<td align="left"><a href="/cgi-bin/matinspector_prof/matrix_help.pl?s=c4c5ecb8a2d5558bda6f8e8b0e555605;NAME=V$ATATA.01" onclick="lnk('/cgi-bin/matinspector_prof/matrix_help.pl?s=c4c5ecb8a2d5558bda6f8e8b0e555605;NAME=V$ATATA.01','MatrixList');return false">V$TBPF/ATATA.01</a></td>
    #<td align="left">Avian C-type LTR TATA box               </td><td>0.78</td>
    #<td>6&nbsp;-&nbsp;22</td><td>(-)</td>
    #<td>1.000</td><td bgcolor="lightgreen">0.906</td><td align="left">tag<font color="red">t</font><font color="red">a</font><font color="red">t</font><font color="red">c</font><font color="red">T</font><font color="red">A</font><font color="red">A</font><font color="red">G</font><font color="red">c</font>tgccg</td>
    #</tr>


    #TABLE HEADINGS
    #Family/matrix 	Further Information 	Opt.	Position	Str.	Core sim. 	Matrix sim.	Sequence

    $j = 2;

    if (!($site_lines[$j++] =~ m/return false">(.*?)<\/a>/))   {$j--; die "bad line: $j: $site_lines[$j]";}
    $AccId = $1;

    if (!($site_lines[$j++] =~ m/^<td align="left">(.*?)<\/td><td>(.*?)<\/td>/))
      {
	$j--;
	$line = join " ",($site_lines[$j++],$site_lines[$j++],$site_lines[$j++]);
	#print STDERR $line."\n";
	if (!($line =~ m/^<td align="left">(.*?)<\/td><td>(.*?)<\/td>/)){
	  $j = $j-3;
	  die "bad line: $j: $line";
	}
	else {
	  $Name = $1;
	  chomp $Name;
	  #$optimal = $2;
	}
      }
    else {
      $Name = $1;
      chomp $Name;
      #$optimal = $2;
    }

    if (!($site_lines[$j++] =~ m/^<td>(.*?)\&nbsp;-&nbsp;(.*?)<\/td><td>(...)<\/td>/))   {$j--; die "bad line: $j: $site_lines[$j]";}
    $StartPos = $1;
    $FinishPos = $2;
    $Direction = $3;

    if (!($site_lines[$j++] =~ m/^<td>(.*?)<\/td><td bgcolor="lightgreen">(.*?)<\/td><td align="left">(.*)$/))   {$j--; die "bad line: $j: $site_lines[$j]";}
    #$core_sim = $1;
    $Confidence = $2;
    $Consensus = $3;
    $Consensus =~ s/<.*?>//g;
    if ($Direction eq '(-)') {
	$Consensus = promoterAnalysis::revcom_keepcase($Consensus);
      }


    #EXAMPLE OUTPUT LINE
    #9	21	V$CHRF/CHR.01	gtatTTGAatctg	Cell cycle gene homology region (CDE/CHR tandem elements regulate cell cycle dependent repression)	(+)	0.92	MatInspector	T,13;T,14;G,15;A,16
    $Program = "MatInspector";
    $sepChar = "\t";
    print $StartPos
      .$sepChar. $FinishPos
	.$sepChar. $AccId
	  .$sepChar. $Consensus
	    .$sepChar. $Name
	      .$sepChar. $Direction
		.$sepChar. $Confidence
		  .$sepChar. $Program
		    .$sepChar;

    promoterAnalysis::printCore($Consensus, $StartPos, ',');
    print "\n";
  }




}

