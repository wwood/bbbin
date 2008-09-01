#!/usr/bin/perl -w
# A script for downloading the relevant information from tracer


$logfile = "wget_log.txt";
open LOG, ">>$logfile" or die "could not append to log file\n";

print LOG "\n\n\n*******************************************************\n";
print LOG "*******************************************************\n";
print LOG "*******************************************************\n";
print LOG `date`;

if ($#ARGV != 0)
{
    print "usage: rget URL\n\n";
    exit();
}


$cmd = 'wget -r --level=3 --http-user=b.woodcroft --http-passwd=ben129 '.$ARGV[0].' &>'.$logfile;


print LOG "$cmd\n";
print LOG "*******************************************************\n";
close LOG;

print `$cmd`;

