#!/usr/bin/env ruby

require 'ensembl'

# add the cds_length method
module Ensembl
  module Core
    class Transcript <DBConnection
      # The length of coding sequence that this transcript encodes
      def cds_length
        protein_seq.length*3
        #coding_region_cdna_end - coding_region_cdna_start
      end
      
      # 'ENSP' + id, instead of ENST + id. Should work
      # for non-human species.
      def stable_protein_id
        if matches = stable_id.match(/(.*)T(.*)/)
          return "#{matches[1]}P#{matches[2]}"
        else
          raise ParseException, "Unable to parse stable transcript id #{stable_id}"
        end
      end
    end
  end
end


# Script
if __FILE__ == $0
  # Is OrthoMCL using v53 like it says on the data sources page?
  include Ensembl::Core
  DBConnection.connect('homo_sapiens',60)
  
  $stdin.each_line do |line|
    whitespaced = line.split(/[,\s]+/)
    whitespaced.each do |line|
      gene_id = line.strip
      
      g = Gene.find_by_stable_id gene_id
      if g.nil?
        $stderr.puts "Couldn't find gene #{gene_id}"
      else
        # ok, good, now, what's the longest protein that this gene encodes?
        transcripts = g.transcripts
        
        transcripts.each do |t|
          puts [t.stable_id, t.cds_length, t.protein_seq].join(' ')
        end
        
        max_length_transcript = transcripts.max {|a,b|
          a.cds_length <=> b.cds_length
        }.cds_length
        
        # Some genes have multiple max_lengths that are the same
        max_transcripts = transcripts.select do |t|
          t.cds_length == max_length_transcript
        end
        p max_transcripts
        
        puts [
        gene_id,
        max_transcripts.collect{|t| t.stable_protein_id}.join(",")
        ].join "\t"
      end
    end
  end
end