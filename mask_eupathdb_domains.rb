#!/usr/bin/env

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
