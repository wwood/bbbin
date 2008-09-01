#/usr/bin/perl -w


@snps = <STDIN>;

print "0 1.5\n";

foreach $i ($ARGV[0]..$ARGV[1]){
  print "$i ";

  if (grep m/^$i\s*$/, @snps){
    print 1;
  }else {
    print 0.5;
  }

  print "\n";
}
