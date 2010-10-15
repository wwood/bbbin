#!/usr/bin/env ruby

class EuPathDB
  class Domain
    attr_accessor :start, :stop, :protein_accession, :domain_accession, :domain_name, :source, :evalue
  end
  
  class DomainAnnotationFile
    def initialize(io)
      @io = io
    end
    
    # Read each line of @io and return a hash of 
    # accession numbers to arrays of domain objects
    def domain_hash
      hash = {}
      @io.each_line do |line|
        domain = parse_domain_from_annotation_line(line)
        unless domain.nil?
          hash[domain.protein_accession] ||= []
          hash[domain.protein_accession].push domain
        end
      end
      hash
    end
    
    # Return a Domain object from a line such as this:
    # (tab separated, but there is extra whitespace on the e-value)
    # PFI1830c  PFAM  PF03011 PFEMP 642 810   1.4E-49
    def parse_domain_from_annotation_line(annotation_line)
      splits = annotation_line.split("\t")
      d = Domain.new
      d.protein_accession = splits[0].strip
      d.source = splits[1].strip
      d.domain_accession = splits[2].strip
      d.domain_name = splits[3].strip
      d.start = splits[4].strip.to_i
      d.stop = splits[5].strip.to_i
      d.evalue = splits[6].strip.to_f
      d
    end
  end
end



# Acting as a script
if $0 == __FILE__
  require 'rubygems'
  require 'bio'
  
  # Input a fasta file
  unless ARGV.length == 2
    $stderr.puts "Usage: #{$0} <fasta_file> <eupathdb_domain_definition_file>"
    exit
  end
  
  fasta_file = ARGV[0]
  domain_file = ARGV[1]
  
  domain_hash = EuPathDB::DomainAnnotationFile.new(File.open(domain_file)).domain_hash
  Bio::FlatFile.open(fasta_file).entries.each do |s|
    sequence = s.seq
    puts ">#{s.definition}"
    
    if domain_hash[s.entry_id]
      domain_hash[s.entry_id].each do |d|
        (d.start..d.stop).each do |n|
          # -1 is required to convert between 1-based indices of the EuPath file
          # and 1 based indices of ruby string representation 
          sequence[n-1..n-1] = 'X'
        end
      end
    end
    puts sequence
  end
end
