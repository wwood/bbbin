#!/bin/bash

cd /blastdb

PLASMODB_VERSION='6.4'
TOXODB_VERSION='6.0'
CRYPTODB_VERSION='4.3'

# soft link the necessary files to the /blastdb folder
# still missing a few species from this section
ln -s ~/phd/data/berghei/genome/PlasmoDB/$PLASMODB_VERSION/PbergheiAllTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta
ln -s ~/phd/data/berghei/genome/PlasmoDB/$PLASMODB_VERSION/PbergheiAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta
ln -s ~/phd/data/berghei/genome/PlasmoDB/$PLASMODB_VERSION/PbergheiGenomic_PlasmoDB-$PLASMODB_VERSION.fasta

ln -s ~/phd/data/Plasmodium\ vivax/genome/PlasmoDB/$PLASMODB_VERSION/PvivaxAnnotatedTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta
ln -s ~/phd/data/Plasmodium\ vivax/genome/PlasmoDB/$PLASMODB_VERSION/PvivaxAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta

ln -s ~/phd/data/falciparum/genome/PlasmoDB/$PLASMODB_VERSION/PfalciparumGenomic_PlasmoDB-$PLASMODB_VERSION.fasta
ln -s ~/phd/data/falciparum/genome/PlasmoDB/$PLASMODB_VERSION/PfalciparumAnnotatedTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta
ln -s ~/phd/data/falciparum/genome/PlasmoDB/$PLASMODB_VERSION/PfalciparumAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta

ln -s ~/phd/data/falciparum/genome/fullmal/5onepass/onepass_pf.fas

ln -s ~/phd/data/Plasmodium\ chabaudi/genome/PlasmoDB/$PLASMODB_VERSION/PchabaudiGenomic_PlasmoDB-$PLASMODB_VERSION.fasta
ln -s ~/phd/data/Plasmodium\ chabaudi/genome/PlasmoDB/$PLASMODB_VERSION/PchabaudiAnnotatedTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta
ln -s ~/phd/data/Plasmodium\ chabaudi/genome/PlasmoDB/$PLASMODB_VERSION/PchabaudiAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta

ln -s ~/phd/data/yoelii/genome/PlasmoDB/$PLASMODB_VERSION/PyoeliiAllTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta 
ln -s ~/phd/data/yoelii/genome/PlasmoDB/$PLASMODB_VERSION/PyoeliiAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta 
ln -s ~/phd/data/yoelii/genome/PlasmoDB/$PLASMODB_VERSION/PyoeliiGenomic_PlasmoDB-$PLASMODB_VERSION.fasta 

# ToxoDB entries
ln -s "/home/ben/phd/data/Toxoplasma gondii/ToxoDB/$TOXODB_VERSION/TgondiiME49Genomic_ToxoDB-$TOXODB_VERSION.fasta"
ln -s "/home/ben/phd/data/Toxoplasma gondii/ToxoDB/$TOXODB_VERSION/TgondiiME49AnnotatedTranscripts_ToxoDB-$TOXODB_VERSION.fasta"
ln -s "/home/ben/phd/data/Toxoplasma gondii/ToxoDB/$TOXODB_VERSION/TgondiiME49AnnotatedProteins_ToxoDB-$TOXODB_VERSION.fasta"

ln -s "/home/ben/phd/data/Neospora caninum/genome/ToxoDB/$TOXODB_VERSION/NeosporaCaninumAnnotatedProteins_ToxoDB-$TOXODB_VERSION.fasta"
ln -s "/home/ben/phd/data/Neospora caninum/genome/ToxoDB/$TOXODB_VERSION/NeosporaCaninumAnnotatedTranscripts_ToxoDB-$TOXODB_VERSION.fasta"
ln -s "/home/ben/phd/data/Neospora caninum/genome/ToxoDB/$TOXODB_VERSION/NeosporaCaninumGenomic_ToxoDB-$TOXODB_VERSION.fasta"

# CryptoDB
ln -s ~/phd/data/Cryptosporidium\ parvum/genome/CryptoDB/$CRYPTODB_VERSION/CparvumAnnotatedProteins_CryptoDB-$CRYPTODB_VERSION.fasta
ln -s ~/phd/data/Cryptosporidium\ parvum/genome/CryptoDB/$CRYPTODB_VERSION/CparvumAnnotatedTranscripts_CryptoDB-$CRYPTODB_VERSION.fasta
ln -s ~/phd/data/Cryptosporidium\ parvum/genome/CryptoDB/$CRYPTODB_VERSION/CparvumGenomic_CryptoDB-$CRYPTODB_VERSION.fasta

ln -s ~/phd/data/Cryptosporidium\ homonis/genome/CryptoDB/$CRYPTODB_VERSION/ChomonisAnnotatedProteins_CryptoDB-$CRYPTODB_VERSION.fasta
ln -s ~/phd/data/Cryptosporidium\ homonis/genome/CryptoDB/$CRYPTODB_VERSION/ChomonisAnnotatedTranscripts_CryptoDB-$CRYPTODB_VERSION.fasta
ln -s ~/phd/data/Cryptosporidium\ homonis/genome/CryptoDB/$CRYPTODB_VERSION/ChomonisGenomic_CryptoDB-$CRYPTODB_VERSION.fasta

# concatenate the databases together
echo "Concatenating the fasta files.."
cat\
 PfalciparumAnnotatedTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta\
 TgondiiME49AnnotatedTranscripts_ToxoDB-$TOXODB_VERSION.fasta\
 NeosporaCaninumAnnotatedTranscripts_ToxoDB-$TOXODB_VERSION.fasta\
 GeneDB_Etenella_Genes\
 PbergheiAllTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta\
 PvivaxAnnotatedTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta\
 PyoeliiAllTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta\
 CparvumAnnotatedTranscripts_CryptoDB-$CRYPTODB_VERSION.fasta\
 >apicomplexa.nucleotide.fa
cat\
 PfalciparumAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta\
 TgondiiAnnotatedProteins_ToxoDB-$TOXODB_VERSION.fasta\
 NeosporaCaninumAnnotatedProteins_ToxoDB-$TOXODB_VERSION.fasta\
# GeneDB_Etenella_Proteins\
 PbergheiAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta\
 PvivaxAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta\
 PyoeliiAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta\
 CparvumAnnotatedProteins_CryptoDB-$CRYPTODB_VERSION.fasta\
 BabesiaWGS.fasta_with_names\
 TANN.GeneDB.pep\
 TPA1.pep\
 >apicomplexa.protein.fa 
cat\
 PfalciparumGenomic_PlasmoDB-$PLASMODB_VERSION.fasta\
 TgondiiME49Genomic_ToxoDB-$TOXODB_VERSION.fasta\
 >apicomplexa.genome.fa


# Create Blast databases
echo "formating apicomplexan-wide databases.."
formatdb -p F -i apicomplexa.nucleotide.fa
formatdb -i apicomplexa.protein.fa
formatdb -p F -i apicomplexa.genome.fa

#species-specific blast databases
echo "formating toxo databases.."
formatdb -i TgondiiME49AnnotatedProteins_ToxoDB-$TOXODB_VERSION.fasta
formatdb -p F -i TgondiiME49AnnotatedTranscripts_ToxoDB-$TOXODB_VERSION.fasta
formatdb -p F -i TgondiiME49Genomic_ToxoDB-$TOXODB_VERSION.fasta

echo "formating neospora databases.."
formatdb -i NeosporaCaninumAnnotatedProteins_ToxoDB-$TOXODB_VERSION.fasta
formatdb -p F -i NeosporaCaninumAnnotatedTranscripts_ToxoDB-$TOXODB_VERSION.fasta
formatdb -p F -i NeosporaCaninumGenomic_ToxoDB-$TOXODB_VERSION.fasta

echo "formating crypto databases.."
formatdb -i CparvumAnnotatedProteins_CryptoDB-$CRYPTODB_VERSION.fasta
formatdb -p F -i CparvumAnnotatedTranscripts_CryptoDB-$CRYPTODB_VERSION.fasta
formatdb -p F -i CparvumGenomic_CryptoDB-$CRYPTODB_VERSION.fasta

echo "formatting berghei"
formatdb -i PbergheiAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta
formatdb -p F -i PbergheiAllTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta
formatdb -p F -i PbergheiGenomic_PlasmoDB-$PLASMODB_VERSION.fasta

echo "formatting yoelii"
formatdb -i PyoeliiAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta
formatdb -p F -i PyoeliiAllTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta
formatdb -p F -i PyoeliiGenomic_PlasmoDB-$PLASMODB_VERSION.fasta

echo "formatting falciparum"
formatdb -i PfalciparumAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta
formatdb -p F -i PfalciparumAnnotatedTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta
formatdb -p F -i PfalciparumGenomic_PlasmoDB-$PLASMODB_VERSION.fasta

echo "formatting chabaudi"
formatdb -i PchabaudiAnnotatedProteins_PlasmoDB-$PLASMODB_VERSION.fasta
formatdb -p F -i PchabaudiAnnotatedTranscripts_PlasmoDB-$PLASMODB_VERSION.fasta
formatdb -p F -i PchabaudiGenomic_PlasmoDB-$PLASMODB_VERSION.fasta

echo "creating Theileria databases.."
cat TANN.GeneDB.pep TPA1.pep >theileria.pep
formatdb -i theileria.pep

echo "creating full malaria one pass database"
cat onepass_pf.fas >fullmal
formatdb -i fullmal -p F


# Create Blat databases
#faToTwoBit apicomplexa.protein.fa apicomplexa.protein.fa.2bit
#faToTwoBit apicomplexa.nucleotide.fa apicomplexa.nucleotide.fa.2bit
#faToTwoBit apicomplexa.genome.fa apicomplexa.genome.fa.2bit

