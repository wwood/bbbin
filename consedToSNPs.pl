#!/usr/bin/perl -w

# This takes the sequence exported from consed (export consensus (with options)) and to convert it into a format readable by gnuplot so I can investigate the SNP density.

$usage = "$0 <qualFile> <start> <stop>\n\n";
if ($#ARGV != 2){die "$usage";}

$qual = $ARGV[0];
$start = $ARGV[1];
$stop = $ARGV[2];


open QUAL, "$qual" or die "Failed to open quality file: $qual";

#$fastaName = <QUAL>;
<QUAL>;

@qualities = <QUAL>;
$all = join ' ', @qualities;
@quals = split ' ',$all;

if ($start<$stop){
  foreach $i ($start-1..$stop-1){
    print "$quals[$i]\n";
  }
} else {
  for ($i=$start; $i>=$stop; $i--){
        print "$quals[$i]\n";
  }
}

