#!/usr/bin/ruby

# Given a gff file in the pseudopipe output format, and a proteome fasta file, generate a fasta file of all the protein sequences

require 'bio'
require 'getoptlong'
require 'rdoc/usage'


# == Usage
#
# hello -f proteome_fasta -g pseudogene_gff

opts = GetoptLong.new(
          [ '--help', '-h', GetoptLong::NO_ARGUMENT ],
          [ '--proteome', '-f', GetoptLong::REQUIRED_ARGUMENT ],
          [ '--pseudogenes', '-g', GetoptLong::REQUIRED_ARGUMENT ]
          )
          
        
# Parse arguments
pseudofile = nil
proteomefile = nil
opts.each do |opt, arg|
  case opt
  when '--help'
    RDoc::usage
    exit
  when '--proteome'
    proteomefile = arg
  when '--pseudogenes'
    pseudofile = arg
  end
end

# Make sure both arguments are present
if !proteomefile || !pseudofile
  RDoc::usage
  exit
end



#the real stuff
# read in the gff file, assigning each gene to an array
gffs = Bio::GFF.new(File.open(pseudofile, 'r').read).records
protein_ids = Hash.new
for g in gffs
  protein_ids[g.frame] = true
end


# read in the fasta file, then go through the sequences looking for the hits
fasta = Bio::FlatFile.open(Bio::FastaFormat, File.open(proteomefile))
entry = fasta.next_entry
while entry
  if protein_ids[entry.entry_id]
    puts entry.to_s
  end


  entry = fasta.next_entry
end
