#!/usr/bin/perl -w

use Bio::Seq;
use Bio::SeqIO;
 
# create a sequence object of some DNA
my $seq = Bio::Seq->new(-id => 'testseq', -seq => 'CATGTAGATAG');
 
# print out some details about it
print "seq is ", $seq->length, " bases long\n";
print "revcom seq is ", $seq->revcom->seq, "\n";
 
# write it to a file in Fasta format 
my $out = Bio::SeqIO->new(-file => '>testseq.fsa', -format => 'Fasta');
$out->write_seq($seq);
