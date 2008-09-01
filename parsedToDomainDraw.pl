#!/usr/bin/perl -w


# Convert a parsed file into a format recognisable by Nick's domain draw format

$sepChar = '&';


if ($#ARGV != 0){
  print "usage: $0 <seqName>\n\n";
  exit;
}
$seqName = $ARGV[0];


#if ($format eq "Match"){
  $startIndex = 0;
  $nameIndex = 2;
  $consensusIndex = 3;
#} elsif ($format eq "MatInspector"){
#  $startIndex = 0;
#  $nameIndex = 2;
#  $consensusIndex = 3;
#} else {
#  print "ERROR: Unrecognised inFormat: $format\n";
#}


$circleIndex = 0;
$borderIndex = 0;

@colourCircle = ('#228b22','#002266','#CD3333','#FF7F00');
@borderCircle = ('Y','R','G');

# protein&domain&range&class&color&Desc.&Method&Veracity&length&url
$i = 0;
foreach $line (<STDIN>){
  @splits = split "\t", $line;

  $finish = $splits[$startIndex]+length($splits[$consensusIndex]);
  print $seqName
    .$sepChar.$splits[$nameIndex].'_'.$i
      .$sepChar.$splits[$startIndex]
	.':'.$finish
	  .'&S&'.
	    $colourCircle[$circleIndex]
	      .'&NULL&MOP&'
		.$borderCircle[$borderIndex]
		  .'&2500&'
	    ."\n";

  #print STDERR $circleIndex;

  #increment, then wrap colour if necessary
  if (++$circleIndex > $#colourCircle){
    $circleIndex = 0;
  }
  if (++$borderIndex > $#borderCircle){
    $borderIndex = 0;
  }

  $i++;
}

