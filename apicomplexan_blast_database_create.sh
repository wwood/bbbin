#!/bin/bash

cd /blastdb

# soft link the necessary files to the /blastdb folder
# still missing a few species from this section
ln -s ~/phd/data/berghei/genome/plasmodb/6.0/PbergheiAllTranscripts_PlasmoDB-6.0.fasta /blastdb
ln -s ~/phd/data/berghei/genome/plasmodb/6.0/PbergheiAnnotatedProteins_PlasmoDB-6.0.fasta /blastdb

ln -s ~/phd/data/vivax/genome/plasmodb/6.0/PvivaxAnnotatedTranscripts_PlasmoDB-6.0.fasta /blastdb
ln -s ~/phd/data/vivax/genome/plasmodb/6.0/PvivaxAnnotatedProteins_PlasmoDB-6.0.fasta /blastdb

ln -s ~/phd/data/falciparum/genome/plasmodb/6.0/PfalciparumGenomic_PlasmoDB-6.0.fasta /blastdb


# concatenate the databases together
cat \
PfalciparumAnnotatedTranscripts_PlasmoDB-6.0.fasta\
 TgondiiAnnotatedTranscripts_ToxoDB-5.2.fasta\
 GeneDB_Etenella_Genes\
 PbergheiAllTranscripts_PlasmoDB-6.0.fasta\
 PvivaxAnnotatedTranscripts_PlasmoDB-6.0.fasta\
 >apicomplexa.nucleotide.fa
cat PfalciparumAnnotatedProteins_PlasmoDB-6.0.fasta TgondiiAnnotatedProteins_ToxoDB-5.2.fasta GeneDB_Etenella_Proteins PbergheiAnnotatedProteins_PlasmoDB-6.0.fasta >apicomplexa.protein.fa 
cat PfalciparumGenomic_PlasmoDB-6.0.fasta >apicomplexa.genome.fa


formatdb -p F -i apicomplexa.nucleotide.fa
formatdb -i apicomplexa.protein.fa
formatdb -p F -i apicomplexa.genome.fa
