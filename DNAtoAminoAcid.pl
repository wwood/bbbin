#!/usr/bin/perl -w

use Bio::SeqIO;

$original = Bio::SeqIO->new(-fh => \*STDIN,
			    '-format' => 'fasta');

$original->alphabet('protein');

#$line = <>;
#chomp $line;
$seq   = $original->next_seq();

#Bio::Seq->new('-seq'=>"$line",
#		      '-moltype' => 'dna');


$aa = $seq->translate()->seq();
print "$aa\n";
