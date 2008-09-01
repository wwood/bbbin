#!/usr/bin/perl -w

use Bio::SeqIO;

# Separate a fasta file into single files, so
# each file has one 1 sequence in it.


my $fastas;
if ($#ARGV == 0){
  $fastas = Bio::SeqIO->new(-file=> $ARGV[0],
			    '-format' => 'fasta');
} else {
  $fastas = Bio::SeqIO->new(-fh=> \*STDIN,
			    '-format' => 'fasta');
}

while ($seq = $fastas->next_seq()){

  my $id = $seq->display_id;

  my $out = Bio::SeqIO->new(-file=>">$id.fa", -format=>'fasta');
  $out->write_seq($seq);

}
