#!/usr/bin/perl -w

# This is a script to take the unparsed Match data and convert it into a
# histogram showing the numbers in each confidence bin


@promoters = ('SoxB2','SoxB','SoxC','SoxF');


# parse into raw histogram data
foreach $p (@promoters){
    $cmd = 'parseMatch.pl 3_Match_'.$p.'_minSum.txt |simpleCutoffHistogram.pl >3_'.$p.'_simpleCutoffHistogram.txt';
    #print $cmd;
    #exit;

    print `$cmd`;
}



# Plot the histograms
$gnuplotScript = "/tmp/histograms.gnuplot"; 
open GNUPLOT, ">$gnuplotScript";
print GNUPLOT 'set terminal png'."\n";
foreach $p (@promoters){
    print GNUPLOT 'set output "3_'.$p.'_simpleCutoffHistogram.png"'."\n";
    print GNUPLOT 'plot "3_'.$p.'_simpleCutoffHistogram.txt" with histeps'."\n";
}


`gnuplot $gnuplotScript`;