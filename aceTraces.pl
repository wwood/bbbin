#!/usr/bin/perl -w

# Get all the sequences from the ace file, including singlets,
# and print them out as a FASTA file.


# Building an input stream
use Bio::Assembly::IO;

my $usage = "usage: $0 <infile.ace>\n";
my $acefile = shift or die $usage;






my $io = new Bio::Assembly::IO(-file=>$acefile, -format=>'ace');


my $assembly = $io->next_assembly();

# Get all the assemblies
foreach $contig ($assembly->all_contigs) {
  my $consensus = $contig->get_consensus_sequence;
  print '>Contig'.$contig->id."\n".$consensus->seq."\n";
}


# Get all the singlets
foreach $contig ($assembly->all_singlets) {
  my $consensus = $contig->get_consensus_sequence;
  print '>Singlet'.$contig->id."\n".$consensus->seq."\n";
}


