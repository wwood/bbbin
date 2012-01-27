#!/usr/bin/env ruby

pros = %w(
PbergheiAnnotatedProteins_PlasmoDB-8.2.fasta
PvivaxAnnotatedProteins_PlasmoDB-8.2.fasta

TgondiiME49AnnotatedProteins_ToxoDB-7.2.fasta

CmurisAnnotatedProteins_CryptoDB-4.6.fasta
ChominisAnnotatedProteins_CryptoDB-4.6.fasta
CparvumAnnotatedProteins_CryptoDB-4.6.fasta

BbovisT2BoAnnotatedProteins_PiroplasmaDB-1.1.fasta
TannulataAnkaraAnnotatedProteins_PiroplasmaDB-1.1.fasta
TparvaMugugaAnnotatedProteins_PiroplasmaDB-1.1.fasta
)

pros.each do |pro|
  Dir.chdir('/blastdb') do
    puts `makeblastdb -in #{pro}`
  end
end

pros.each do |pro|
  puts `blastp -query '/home/ben/phd/data/Plasmodium falciparum/genome/PlasmoDB/8.2/PfalciparumAnnotatedProteins_PlasmoDB-8.2.fasta' -db /blastdb/#{pro} -evalue 1e-5 -out /home/ben/bin/plasmarithm/input_data/falciparum8.2versus#{pro}.blast.csv -outfmt 6`
end
