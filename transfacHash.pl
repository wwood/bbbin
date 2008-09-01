#!/usr/bin/perl -w


#<A HREF="/pub/databases/transfac/doc/factor2.html#MX" target="_new">MX</A>   <A HREF="/cgi-bin/pub/databases/transfac/getTF.cgi?AC=m00251">M00251</A> V$XBP1_01.


## THIS BIT WAS MISGUIDED, BUT WORKS, if I want the matrix names from the html files
#opendir DIRECTORY, ".";

#while (defined($file = readdir DIRECTORY)) {
#  open HTML, $file;

#  if (!$file =~ m/.*\.html$/) {
#    next;
#  }

#  foreach $line (<HTML>) {
#    chomp $line;
#    # If it is a matrix line, extract the ID
#    if ($line =~ m/MX</) {
#      #print $line."\n\n";
#      if ($line =~ m/getTF.cgi\?AC=(.*)\">(.*)<\/A> (.*)\./) {
#      #if ($line =~ m/getTF.cgi?AC=(.*)\">(.*)</) {
#	print "$file $3 $2 $1\n";
#      } else {
#	print STDERR "PROBLEM $file: $line\n";
#      }
#      #exit;
#    }
#  }
#}



opendir DIRECTORY, ".";

while (defined($file = readdir DIRECTORY)) {
  open HTML, $file;

  if (!$file =~ m/.*\.html$/) {
    next;
  }

  foreach $line (<HTML>) {
    chomp $line;
    # If it is a matrix line, extract the ID
    if ($line =~ m/>FA<\/A>(.*)/) {

      $name = $1;
      $name =~ s/\s*//g;

      $file =~ s/\.html//;

      print "$name $file";
    }

    # If it is a matrix line, extract the ID
    if ($line =~ m/MX</) {
      if ($line =~ m/getTF.cgi\?AC=(.*)\">(.*)<\/A> (.*)\./) {
	print " $3";
      } else {
	print STDERR "PROBLEM $file: $line\n";
      }
    }
  }

  print "\n";
}
