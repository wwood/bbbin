#!/usr/bin/perl -w

use Bio::SeqIO;

$original = Bio::SeqIO->new(-fh => \*STDIN,
			    '-format' => 'fasta');

$original->alphabet('dna');

#$line = <>;
#chomp $line;
while (defined($seq   = $original->next_seq())){

#Bio::Seq->new('-seq'=>"$line",
#		      '-moltype' => 'dna');


$aa = $seq->translate()->seq();
print ">".$seq->display_id."\n";
print "$aa\n";
}
