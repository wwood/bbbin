#!/usr/bin/env ruby

# Extract a sequence surrounding a particular gene model in EuPathDB. Requires the GFF file and the genomic sequence file to be handy.
require 'rubygems'
require 'bio'
require 'reubypathdb'

if __FILE__ == $0
  options = {
  :genomic_sequence_fasta_filename => "/home/ben/phd/data/Toxoplasma gondii/genome/ToxoDB/6.3/TgondiiGT1Genomic_ToxoDB-6.3.fasta",
  :gff_filename => "/home/ben/phd/data/Toxoplasma gondii/genome/ToxoDB/6.3/TgondiiGT1_ToxoDB-6.3.gff",
  :organism_name => 'Toxoplasma_gondii_GT1'
  }
  unless ARGV.length == 1
    $stderr.puts "Usage: eupathdb_genomic_sequence_extractor.rb <gene_model_id>"
    exit
  end
  gene_id_of_interest = ARGV[0]
  
  # Sort through the gff file, looking for the correct ID.
  gff = EupathDBGFF.new(options[:gff_filename])
  gene_of_interest = gff.next_gene
  while gene_of_interest.name.gsub('apidb|','') != gene_id_of_interest
    gene_of_interest = gff.next_gene
    break if gene_of_interest.nil?
  end
  raise Exception, "Unable to find gene ID #{gene_id_of_interest}" if gene_of_interest.nil?
  
  start = nil
  stop = nil
  if gene_of_interest.strand == '+'
    start = gene_of_interest.cds[0].from
    stop = gene_of_interest.cds[gene_of_interest.cds.length-1].to
  elsif gene_of_interest.strand == '-'
    start = gene_of_interest.cds[gene_of_interest.cds.length-1].from
    stop = gene_of_interest.cds[0].to
  else
    raise
  end
  
  $stderr.puts "Found a gene on chromosome '#{gene_of_interest.seqname}' called '#{gene_of_interest.name}', from #{start} to #{stop}"
  p gene_of_interest
  
  # Extract the genomic coordinates
  found = false
  seq_name_naked = gene_of_interest.seqname.gsub(/^>[^|]+|/,'')
  Bio::FlatFile.foreach(options[:genomic_sequence_fasta_filename]) do |seq|
    if matches = seq.entry_id.match(/^>[^|]+|([a-zA-Z0-9_]+)$/)
      p [matches[1], seq_name_naked]
      if matches[1] == seq_name_naked #if we found the chromosome
        puts "found!"
        # extract the genomic region
        puts seq.seq.length
      end
    else
      raise Exception, "Couldn't parse scaffold id #{seq.entry_id}"
    end
  end
end