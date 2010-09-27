#!/usr/bin/env

class EuPathDB
  class Domain
    attr_accessor :start, :stop, :protein_accession, :domain_accession, :source, :evalue
  end
  
  class DomainAnnotationFile
    def initialize(io)
      @io = io
    end
    
    # Read each line of @io and return a hash of 
    # accession numbers to arrays of domain objects
    def domain_hash
      hash = {}
      io.each_line do |line|
        domain = parse_domain_from_annotation_line(line)
        hash[domain.protein_accession] ||= []
        hash[domain.protein_accession].push domain
      end
      hash
    end
    
    # Return a Domain object from a line such as this:
    # (tab separated, but there is extra whitespace on the e-value)
    # PFI1830c  PFAM  PF03011 PFEMP 642 810   1.4E-49
    def parse_domain_from_annotation_line(annotation_line)
      d = Domain.new
      d.protein_accession = annotation_line[0].strip
      d.source = annotation_line[1].strip
      d.domain_accession = annotation_line[2].strip
      d.start = annotation_line[3].strip.to_i
      d.stop = annotation_line[4].strip.to_i
      d.evalue = annotation_line[5].strip
    end
  end
end
