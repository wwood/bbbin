# Joins lines in two files based on a shared column value. Will also return
# the lines not joined from the first file. Program takes 4 arguments: file1, file2, column to match in file 1, column to match in file2

$f1 = $ARGV[0];
$f2 = $ARGV[1];
$col1 = $ARGV[2];
$col2 = $ARGV[3];

$usage = "to run type: merge.pl
<file1> <file2> <file1 column no.> <file2 column no.>\n";
if (@ARGV != 4) { die $usage; }

# enter lines from 2nd file into hash where keys are the values in the
# column to be matched
open(F2,$f2) || die "Can't open inputfile $f2! $!\n";
while (<F2>) {
    s/\r?\n//;
    @F=split /\t/, $_;
    $seen{$F[$col2]} .= "$_\n"
};

# match values from column specified in 1st file to values in 2nd file
open(F1,$f1) || die "Can't open inputfile $f1! $!\n";
while (<F1>) {
    s/\r?\n//;
    @F=split /\t/, $_;
    $x = $seen{$F[$col1]};
    if ($x) {
        $x =~ s/^/$_\t/gm;
        print "$x\n"; #print matched lines from file 1 and file2
    }
    else {print "$_\n";}  # if no match to file 1 value print unmatched value  
}

