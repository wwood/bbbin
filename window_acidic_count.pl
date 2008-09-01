#!/usr/bin/perl -w
#this is a program that counts the no. of basic and acidic aas in a file of seqs
#seq that has the highest no. of bases and lowest basic/acidic aa ratio
#INPUT:       perl window.pl --infile aaatest --winsize 20
#OUTPUT:      SeqID Desc WindowSeq PosOfWindowInSeq(start) No.BasicAAs
#             No.AcidicAAs Basic/AcidicRatio

#use lib "/home/leishcyc/lib/perl5/site_perl/5.8.5/";
#use Bio::Perl;

use Bio::SeqIO;
use Getopt::Long;


my $infile = "";
my $winsize = 10;		#default
GetOptions('infile=s' => \$infile, 'winsize=i' => \$winsize);
die("Usage = window.pl --infile  <fasta file> [--winsize <int>]\n") if(!$infile);

my $seqio = Bio::SeqIO->new('-file' => $infile, '-format' => 'fasta');

print "ID\tBasic\tAcidic\tB/A Ratio\n";

my $sep = "\t";

while (my $seqobj = $seqio->next_seq) {
  print $seqobj->display_id.$sep;

  # Preferably a whole window, otherwise the whole length of the sequence
  # if there isn't enough of it
  my $size = $winsize < $seqobj->length ? $winsize : $seqobj->length;

  my $subseq = $seqobj->subseq(1,$size);

  my $bcount = 0;
  my $acount = 0;

  foreach $i (0..(length $subseq)-1) {
    my $aa = substr $subseq, $i, 1;
    if ($aa=~/[KRH]/ig) {
      $bcount++;
    }				# counts basic aas in window
    if ($aa=~/[DE]/ig) {
      $acount++;
    }				# counts acidic aas in window
  }
  my $ratio = $acount == 0 ? $bcount : $bcount/$acount; #acount==0 means a divide by zero problem.
  print $bcount.$sep.$acount.$sep.$ratio."\n";
}
