#!/usr/bin/perl -w

# Takes in an HTML conreal output, and parses it into a recognisable output format



#</pre>
#  <h3>CONREAL ALIGNED HITS</h3><pre><table border="1">
#<tbody><tr>
#  <td>Matrix</td>
#  <td>Identity</td>
#  <td>Seq1 from</td>
#  <td>Seq1 to</td>
#  <td>Seq2 from</td>
#  <td>Seq2 to</td>
#  <td>Strand</td>
#  <td>Score1</td>
#  <td>Score2</td>
#  <td>Rel.Score1</td>
#  <td>Rel.Score2</td>
#  <td>Factor</td>
#  <td>Found by:<br>C = CONREAL<br>L = LAGAN<br>M = MAVID<br>B = BLASTZ<br></td>
#</tr>
#
#<tr>
#  <td><a href="http://jaspar.cgb.ki.se/cgi-bin/jaspar_db.pl?rm=present&amp;ID=MA0033">MA0033</a></td>
#  <td>52.63</td>
#  <td>74</td>
#  <td>81</td>
#  <td>361</td>
#  <td>368</td>
#  <td>1</td>
#  <td>5.104</td>
#  <td>4.654</td>
#  <td>0.84</td>
#  <td>0.82</td>
#  <td>FREAC-7</td>
#  <td>C only</td>
#</tr>



$mode = 'BEGIN';
@fields = ('Matrix','Identity',
	   'Seq1 from','Seq1 to'
	   ,'Seq2 from'
	   ,'Seq2 to'
	   ,'Strand'
	   ,'Score1'
	   ,'Score2'
	   ,'Rel.Score1'
	   ,'Rel.Score2'
	   ,'Factor'
	   ,'Found by');
$fieldCount = 0;
$recordCount = 0;

foreach $line (<STDIN>) {


  if ($mode eq 'BEGIN') {
    if ($line =~ m/<h3>CONREAL ALIGNED HITS</) {
      $mode = 'CONREAL_ENTRY';
    }
  }

  elsif ($mode eq 'CONREAL_ENTRY') {
    if ($line =~ m/<\/tr>/) {
      $mode = 'CONREAL';
    }
  }

  elsif ($mode eq 'CONREAL') {
    if ($line =~ m/<td>(.*)</) {
      $records[$recordCount][$fieldCount++] = $1;
      if ($fieldCount > $#fields) {
	$fieldCount = 0;
	$recordCount++;
      }

      elsif ($line =~ m/ALIGNED HITS</) {
	$mode = 'LAGAN';
      }
    }
  }
}
