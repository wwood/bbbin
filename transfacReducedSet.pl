#!/usr/bin/perl -w

# This script is for parsing in from pwmMatching.doc the TRANSFAC entries that 
# we want to keep in the reduced set.
# Output is a list of group descriptions, classification numbers and descriptions of entries that should stay in the set.


#INPUT FILE STRUCTURE
#
# Blank lines OK
# Comment lines start with '#'
# Groupings start with '*'
# All other lines must start with a number dot number dot etc and then description.
# Lines with Class: or Family: or Subfamily: are ignored.


$group = "Unassigned";

@classificationLineParts = ('Class:','Family:','Subfamily:');



$line_number = 1;
foreach $line (<STDIN>) {
  chomp $line;
	
  $line =~ s/^\s*//;
  $line =~ s/\s*$//;

  # Blank lines and comment lines
  if ($line eq "" or $line =~ m/^\#/) {
    # Do nothing
  }

  # Group lines
  elsif ($line =~ m/^\*(.*)/) {
    $group = $1;
  }

  #Classification lines
  elsif ($line =~ m/$classificationLineParts[0]/ ||
	 $line =~ m/$classificationLineParts[1]/ ||
	 $line =~ m/$classificationLineParts[2]/) {
    print "ignored: $line\n";
  }

  #TF match line
  elsif ($line =~ m/(\d\.)*\d\s+(.*)$/) {
    $desc = $+;
    $line =~ m/^(\S*)/;
    $id = $1;

    print "$group\t$desc\t$id\n";
  } else {
    print STDERR "Bad line:$line_number: \"$line\"\n";
  }

  $line_number++;
}
