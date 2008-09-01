#/usr/bin/perl -w
#/opt/local/bin/perl -w


# given a csv file documenting SNPs, check for typographical or other such errors



#Example
#4021,snp,t,*,1,1
#4831,snp,g,*,2,1
#4846,snp,a,c,0,2
#4848,snp,a,a,0,2
#4872,snp,t,a,0,3
#4877,snp,a,*,1,1
#4885,snp,t,g,0,2
if ($#ARGV != 2){
  print "Usage: $0 snpCsvFile beginNum endNum\n\n";
  exit;
}

$csvFile = $ARGV[0];
$last_num = $ARGV[1];
$end_num = $ARGV[2];

open IN, "$csvFile" or die "couldn't open file: $csvFile";
print "ANALYSING: $csvFile\n";

foreach $line (<IN>){
  chomp $line;

  #skip blank lines
  if ($line eq ""){
    next;
  }

  $line =~ s/[^\w,\*]//g;
  #print "fds: $line\n";
  @splits = split ',', $line;
  

  if ($#splits != 5 &&
      $#splits != 6){
    print "Not enough commas: $line\n";
    next;
  }

  $num = $splits[0];
  $type = $splits[1];
  $out_base = $splits[2];
  $in_base = $splits[3];
  $top = $splits[4];
  $bottom = $splits[5];


  # Debugging
  #foreach (0..5){
  #  print "split[$_]: $splits[$_]\n";
  #}


  # Check that the current number is possible
  if ($num < $last_num || $num > $end_num){
    print "Bad number detected: $line\n";
  }

  # check that the 2nd field is either 'between' or 'snp'
  if (!($type =~ m/^(snp|between)$/)){
    print "bad 2nd col: $line\n";
  }

  # Check the 3rd col and the 4th col are a,t,g,c, or *, 
  # and that they are not the same
  if (!($in_base =~ m/^[atgc\*]$/) ||
      !($out_base =~ m/^[atgc\*]$/) ||
      $in_base eq $out_base){
    print "Bad line in/out: $line\n";
  }

  # Check there were at least 2 hits
  if ($top + $bottom <2){
    print "Bad line how manies: $line\n";
  }


  #exit;
}


print "\n";
