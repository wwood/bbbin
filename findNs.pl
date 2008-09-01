#!/usr/bin/perl -w

use Bio::SeqIO;


#$tofind = 'N';

my $fastas;
if ($#ARGV == 0){
  $fastas = Bio::SeqIO->new(-file=> $ARGV[0],
			    '-format' => 'fasta');
} else {
  $fastas = Bio::SeqIO->new(-fh=> \*STDIN,
			    '-format' => 'fasta');
}

while ($seq = $fastas->next_seq()){
  $me = $seq->seq();
  #print $me."\n";;

  print '>'.$seq->display_id()."\n";

  $offset = 0;


  while (length $me != -1){

    if ($me =~ s/([atgcATGC]+)(N+)//){
      
      #print "$1..$2\n";

      $a = $1;
      $b = $2;

      $o = $offset+1+length($a); # +1 because 1-based indices for start but not end
      print $o.' - ';
      $offset += length $a;
      $p = $offset+length($b);
      print $p;
      $offset += length $b;

      $d = $p-$o+1;
      print " ($d)\n";

      #print "o$offset\n";
      #exit;
    }
    else {last;}
  }
}

print "\n";
