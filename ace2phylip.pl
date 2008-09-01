#!/usr/bin/perl -w

use Bio::Assembly::IO;

if ($#ARGV != 0)
  {
    print "usage: ace2phylip <aceFile>\n";
    print "output file is piped out\n\n";
    exit(0);
  }
$ace_file = $ARGV[0];

$io = new Bio::Assembly::IO(-file=>$ace_file,
			    -format=>"ace");

print $io;

$assembly = $io->next_assembly;

#@singlets = $assembly->all_singlets();

foreach $i (0..2)
  {
#    print $singlets[$i];
  }
