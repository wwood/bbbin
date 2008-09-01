#/usr/bin/perl -w


# For plotting SNPS directly from the csv file


@csvs = <STDIN>;

$first = $ARGV[0];
$second = $ARGV[1];

# Go through each of the lines, looking for between ones.
# If a between one is found, push it onto the stack. Keep track
# of the number of betweens and add accordingly to the ones after that
$num_betweens = 0;
foreach $line (@csvs){
  $line =~ m/^(\d*)/;
  $num = $1;

  if ($line =~ m/between/){
    $num_betweens++;
  }

  push @snps, $num+$num_betweens;
}


print "0 20\n";


print STDERR ".$first.$second.\n";

if ($first < $second) {
  print STDERR "normal mode\n";
  foreach $i ($first..($second+$num_betweens)) {
    $num = $i-$first;
    print "$num ";


    if (grep m/^$i\s*$/, @snps) {
      print 10;
    } else {
      print 0.01;
    }

    print "\n";
  }
}

else {
#  for ($i=$first;$i>=$second;$i--)
  print STDERR "reverse mode\n";
  foreach $i ($second..($first+$num_betweens)){
    $num = $i-$second;
    print "$num ";

    if (grep m/^$i\s*$/, @snps) {
      print 1;
    } else {
      print 0.5;
    }

    print "\n";
  }
}
