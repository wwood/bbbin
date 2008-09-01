#!/usr/bin/perl -w

# A script to take a HMM library, and then make multiple single files
# out of it.


#EXAMPLE:
# ...
#     -      *      *      *      *      *      *      *      *      *      *      *      *   
#   *      *      *      *      *      *      *      * 
#     -      *      *      *      *      *      *      *      *      0 
#//
#HMMER2.0  [2.3.1]
#NAME  2-Hacid_dh
#ACC   PF00389.21
# ...



open LIB, "<$ARGV[0]" or die "could not open HMM library";
my $count = 0;

open OUTPUT, '>current.pfam' or die "could not open dummy output file";
while ($line = <LIB>){
    print OUTPUT $line;
    #if the line tells you the accession number of the current HMM
    if ($line =~ m/^ACC\s+(.+)/){
	$line =~ m/ACC\s+(PF.+)\.(\d+)/ or die "line didn't match: '$line'";
	$cur_acc = $1;
    } elsif ($line =~ m/^\/\/$/){
	#if a separator between HMMs
	close OUTPUT;
	if (!defined($cur_acc)){
	    die "undefined accession in the last HMM, or this script is too stupid";
	}
	print `mv current.pfam $cur_acc`;
	open OUTPUT, '>current.pfam' or die "couldn't open dummy file after the first time";
	$cur_acc = undef;
	$count++;
    }
}

#don't have finish the last one off, because it ends with a '//' already


#checking mechanisms
`rm current.pfam`;
print  "NEXT ONES SHOULD BE EQUAL IF THE CWD WAS EMPTY\n";
print "total processed: $count\n";
print "num files in directory: ".`ls |wc -l`."\n";



close LIB;
