#!/usr/bin/perl -w
#this is a program that counts the no. of basic and acidic aas in a file of seqs
#seq that has the highest no. of bases and lowest basic/acidic aa ratio
#INPUT:       perl window.pl --infile aaatest --winsize 20
#OUTPUT:      SeqID Desc WindowSeq PosOfWindowInSeq(start) No.BasicAAs
#             No.AcidicAAs Basic/AcidicRatio

#use lib "/home/leishcyc/lib/perl5/site_perl/5.8.5/";
use Bio::Perl;

use Bio::SeqIO;
use Getopt::Long;

my($infile, $winsize, $seqio, $seqobj, $id, $len, $i, $window, $bcount);
my($acount, $aa, $ratio, $result, @prot, @resultarray, @sortedresult);

my $infile = "";
my $winsize = 10;  #default
GetOptions('infile=s' => \$infile, 'winsize=i' => \$winsize);
die("Usage = window.pl --infile  <fasta file> [--winsize <int>]\n") if(!$infile);

my $seqio = Bio::SeqIO->new('-file' => $infile, '-format' => fasta);

print "ID\tDesc\tWindow\tStart\tBasic\tAcidic\tB/A Ratio\n";

while(my $seqobj = $seqio->next_seq) {
    my $id  = $seqobj->display_id;
    my $desc  = $seqobj->desc;
    my $len = $seqobj->length;
    my @resultarray = ();

# for the length of the seq make windows of size: winsize 
for(my $i = 1; $i <= ($len-($winsize-1)); $i++) {
    my $window = $seqobj->subseq($i,$i+($winsize-1));
    if($window =~/[KRH]/ig) {  # select only windows that have a basic aa 
        foreach ($window) {
        my $bcount = 0;  # sets count of basic aas to zero for every line
        my $acount = 0;  # sets count of acidic aas to zero for every line
        my @prot=split("", $window);
        while(@prot) {
           my $aa=shift(@prot); 
            if($aa=~/[KRH]/ig)
            {$bcount++;} # counts basic aas in window
            if($aa=~/[DE]/ig) 
            {$acount++;} # counts acidic aas in window

        }
# if window has acidic aas calculate the basic/acidic ratio 
       if ($acount != 0){
       $ratio=$bcount/$acount;
       } else {
       $ratio = 0;
       }
       $result = "$id\t$desc\t$window\t$i\t$bcount\t$acount\t$ratio\n"; 
 
# put the result for every window into an array 
      push @resultarray, $result; 
    }
}
}
# extract the top-scoring window for the seq = highest no. of bases, lowest
# basic/acidic ratio and closest to start of seq
@sortedresult = map { $_->[0] } # whole line
        sort {
              $b->[5] <=> $a->[5]  # no. bases
                      ||
              $a->[7] <=> $b->[7]  # b/a ratio
                      ||
              $a->[4] <=> $b->[4]  # seq pos for start of window
           }
      map  { [ $_, (split /\t/) ] }
      @resultarray;
      print $sortedresult[0]; # print top-scoring window

}
