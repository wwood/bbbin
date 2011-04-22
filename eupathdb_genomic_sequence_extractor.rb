#!/usr/bin/env ruby

# Extract a sequence surrounding a particular gene model in EuPathDB. Requires the GFF file and the genomic sequence file to be handy.
require 'rubygems'
require 'bio'
require 'reubypathdb'
require 'optparse'

if __FILE__ == $0
  options = {
  :genomic_sequence_fasta_filename => "/home/ben/phd/data/Toxoplasma gondii/genome/ToxoDB/6.3/TgondiiGT1Genomic_ToxoDB-6.3.fasta",
  :gff_filename => "/home/ben/phd/data/Toxoplasma gondii/genome/ToxoDB/6.3/TgondiiGT1_ToxoDB-6.3.gff",
  :organism_name => 'Toxoplasma_gondii_GT1',
  :upstream_length => 200,
  }
  
  o = OptionParser.new do |opts|
    opts.banner = [
      'Usage: program [-f genomic_sequence_fasta_filename] [-g gff_filename] [-o organism_name]'
    ]
    opts.on("-f", "--genome_fasta_file FASTA_FILENAME", "A fasta file of the genomic sequence") do |f|
      options[:genomic_sequence_fasta_filename] = f
    end
    opts.on("-g", "--genome_gff_file GFF_FILENAME", "A GFF file of the genome") do |f|
      options[:gff_filename] = f
    end
    opts.on("-o", "--organism_name ORGANISM_NAME", "The name of the organism used in the fasta file - e.g. 'Toxoplasma_gondii_GT1'") do |f|
      options[:organism_name] = f
    end
    opts.on("-l", "--length UPSTREAM_LENGTH", "How much sequence to extract 5' of the stop codon (including the stop codon)") do |f|
      options[:upstream_length] = f.to_i
    end
  end
  o.parse!
  
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
  
  # Where in the chromosome fasta sequence should be extracted?
  start = nil
  stop = nil
  if gene_of_interest.strand == '+'
    stop = gene_of_interest.cds[gene_of_interest.cds.length-1].to.to_i
    start = stop-options[:upstream_length]+2
  elsif gene_of_interest.strand == '-'
    start = gene_of_interest.cds[gene_of_interest.cds.length-1].from.to_i+1
    stop = start+options[:upstream_length]-2
  else
    raise
  end
  
  $stderr.puts "Found a gene on chromosome '#{gene_of_interest.seqname}' called '#{gene_of_interest.name}'. Extracting #{start} to #{stop}."
  
  # Extract the genomic coordinates
  found = false
  seq_name_naked = gene_of_interest.seqname.gsub(/^[^|]+\|/,'')
  Bio::FlatFile.foreach(options[:genomic_sequence_fasta_filename]) do |seq|
    if matches = seq.entry_id.match(/^>[^|]+|([a-zA-Z0-9_]+)$/)
      if matches[1] == seq_name_naked #if we found the chromosome
        # extract the genomic region
        s = seq.seq[start-2..stop-1]
        puts ">#{gene_id_of_interest}_#{options[:upstream_length]}_upstream_from_stop_codon"
        if gene_of_interest.strand == '-'
          s = Bio::Sequence::NA.new(s).reverse_complement.to_s.upcase
        end
        puts s
      end
    else
      raise Exception, "Couldn't parse scaffold id #{seq.entry_id}"
    end
  end
end