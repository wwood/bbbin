#!/bin/bash

cd /blastdb

# soft link the necessary files to the /blastdb folder
# still missing a few species from this section
ln -s ~/phd/data/berghei/genome/plasmodb/6.0/PbergheiAllTranscripts_PlasmoDB-6.0.fasta
ln -s ~/phd/data/berghei/genome/plasmodb/6.0/PbergheiAnnotatedProteins_PlasmoDB-6.0.fasta

ln -s ~/phd/data/vivax/genome/plasmodb/6.0/PvivaxAnnotatedTranscripts_PlasmoDB-6.0.fasta
ln -s ~/phd/data/vivax/genome/plasmodb/6.0/PvivaxAnnotatedProteins_PlasmoDB-6.0.fasta

ln -s ~/phd/data/falciparum/genome/plasmodb/6.0/PfalciparumGenomic_PlasmoDB-6.0.fasta

ln -s ~/phd/data/yoelii/genome/plasmodb/6.0/PyoeliiAllTranscripts_PlasmoDB-6.0.fasta 
ln -s ~/phd/data/yoelii/genome/plasmodb/6.0/PyoeliiAnnotatedProteins_PlasmoDB-6.0.fasta 

ln -s "/home/ben/phd/data/Toxoplasma gondii/ToxoDB/5.2/TgondiiME49Genomic_ToxoDB-5.2.fasta"
ln -s "/home/ben/phd/data/Toxoplasma gondii/ToxoDB/5.2/TgondiiME49AnnotatedTranscripts_ToxoDB-5.2.fasta"
ln -s "/home/ben/phd/data/Toxoplasma gondii/ToxoDB/5.2/TgondiiME49AnnotatedProteins_ToxoDB-5.2.fasta"


# concatenate the databases together
echo "Concatenating the fasta files.."
cat\
 PfalciparumAnnotatedTranscripts_PlasmoDB-6.0.fasta\
 TgondiiAnnotatedTranscripts_ToxoDB-5.2.fasta\
 GeneDB_Etenella_Genes\
 PbergheiAllTranscripts_PlasmoDB-6.0.fasta\
 PvivaxAnnotatedTranscripts_PlasmoDB-6.0.fasta\
 PyoeliiAllTranscripts_PlasmoDB-6.0.fasta\
 >apicomplexa.nucleotide.fa
cat\
 PfalciparumAnnotatedProteins_PlasmoDB-6.0.fasta\
 TgondiiAnnotatedProteins_ToxoDB-5.2.fasta\
 GeneDB_Etenella_Proteins\
 PbergheiAnnotatedProteins_PlasmoDB-6.0.fasta\
 PvivaxAnnotatedProteins_PlasmoDB-6.0.fasta\
 PyoeliiAnnotatedProteins_PlasmoDB-6.0.fasta\
 >apicomplexa.protein.fa 
cat\
 PfalciparumGenomic_PlasmoDB-6.0.fasta\
 TgondiiME49Genomic_ToxoDB-5.2.fasta\
 >apicomplexa.genome.fa


# Create Blast databases
echo "formating apicomplexan-wide databases.."
formatdb -p F -i apicomplexa.nucleotide.fa
formatdb -i apicomplexa.protein.fa
formatdb -p F -i apicomplexa.genome.fa

#species-specific blast databases
echo "formating toxo databases.."
formatdb -i TgondiiME49AnnotatedProteins_ToxoDB-5.2.fasta
formatdb -p F -i TgondiiME49AnnotatedTranscripts_ToxoDB-5.2.fasta
formatdb -p F -i TgondiiME49Genomic_ToxoDB-5.2.fasta

# Create Blat databases
#faToTwoBit apicomplexa.protein.fa apicomplexa.protein.fa.2bit
#faToTwoBit apicomplexa.nucleotide.fa apicomplexa.nucleotide.fa.2bit
#faToTwoBit apicomplexa.genome.fa apicomplexa.genome.fa.2bit

