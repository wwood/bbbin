#!/usr/bin/perl -w

# Takes a FASTA file that has alignments in it and gives back just the first sequence in the alignment.
# An example of these FASTA files are the ones that are output from the pseudogene.org pipeline.
#
# EXAMPLE
#>chr9q_ENSCINP00000028892.1  NoBandData  2574153  2574284  -  ENSCINP00000028892  1  44  45  0.98  0  0  0  0  1e-15  1.000  0  0  1  com(2574153..2574284)  132  .  .  GENE-SINGLE  GENE-SINGLE      
#TPKYI  AYTTKSNLRFLLFAKRLIIASGSITISLCWCILFFNLKV
#TPKYI  AYTTKSNLRFLLFAKRLIIASGSITISLCWCILFFNLKV
#
# Assumes that the order of the file lines are 
# 1. fasta name
# 2. first sequence
# 3. second sequence
#

while (my $description = <STDIN>){
    my $first = <STDIN>;
    my $second = <STDIN>;
    
    #remove the junk from the sequence
    $first =~ s/\-*\s*//g;
    $second =~ s/\-*\s*//g;

    print $description;
    print $first;
    print "\n";
}
