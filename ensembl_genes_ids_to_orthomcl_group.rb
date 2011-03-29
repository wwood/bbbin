#!/usr/bin/env ruby


# A script that takes a list of ENSEMBL gene identifiers
# and maps them to OrthoMCL groups. Requires the orthomcl 
# groups file to be downloadable locally.
#
# Only tested on Ruby 1.9.2!

require 'ensembl'

# add the cds_length method
module Ensembl
  module Core
    class Transcript <DBConnection
      # The length of coding sequence that this transcript encodes
      def cds_length
        #protein_seq.length*3
        return 0 unless biotype == 'protein_coding'
        coding_region_cdna_end - coding_region_cdna_start
      end
    end
  end
end

module Bio
  class OrthoMCL
    class Group
      attr_accessor :group_id
      attr_accessor :genes
      
      def self.create_from_groups_file_line(groups_file_line)
        group = self.new
        if matches = groups_file_line.match(/(OG[\d_]+):(.+)/)
          group.group_id = matches[1]
          group.genes = matches[2].strip.split(' ') 
        else
          raise Exception, "Failed to parse OrthoMCL line #{groups_file_line}"
        end
        return group
      end
    end
  end
end


# Script
if __FILE__ == $0
  orthomcl_groups_gzip_filename = ARGV[0]
  orthomcl_groups_gzip_filename ||= '/home/ben/phd/data/orthomcl/v4/groups_OrthoMCL-4.txt.gz' #hack
  
  # Is OrthoMCL using v53 like it says on the data sources page?
  include Ensembl::Core
  DBConnection.connect('homo_sapiens',60)
  
  $stdin.each_line do |line|
    whitespaced = line.split(/[,\s]+/)
    whitespaced.each do |line|
      gene_id = line.strip
      $stderr.puts "Processing #{gene_id} .."
      
      g = Gene.find_by_stable_id gene_id
      if g.nil?
        $stderr.puts "Couldn't find gene #{gene_id}"
      else
        # ok, good, now, what's the longest protein that this gene encodes?
        transcripts = g.transcripts
        
        # Find the maximum length of the CDS
        max_length_transcript = transcripts.max {|a,b|
          a.cds_length <=> b.cds_length
        }.cds_length
        
        # Some genes have multiple max_lengths that are the same, I'm unsure what OrthoMCL's policy on this is.
        max_transcripts = transcripts.select do |t|
          t.cds_length == max_length_transcript
        end
        
        # Convert to Protein IDs, since that is what OrthoMCL uses
        max_translation_ids = max_transcripts.collect{|t| t.translation.stable_id}
        
        # Extract the group IDs for each Ensembl Protein ID
        groups = []
        max_translation_ids.each do |t|
          liner = `zcat #{orthomcl_groups_gzip_filename} |grep #{t}`.strip
          if liner.match(/\n/)
            $stderr.puts "Multiple groups found for protein ID #{t}! Whack."
          else
            unless liner.nil?
              groups.push Bio::OrthoMCL::Group.create_from_groups_file_line liner
            end
          end
        end
        
        # error checking
        if groups.empty?
          $stderr.puts "No OrthoMCL groups found for gene id #{gene_id}, using protein ID(s)#{max_translation_ids.join(',')}"
        elsif groups.length > 1
          $stderr.puts "Too many OrthoMCL groups found for gene id #{gene_id}, using protein ID(s)#{max_translation_ids.join(',')}"
        else
          puts [
          gene_id,
          max_translation_ids.join(','),
          groups[0].group_id
          ].join "\t"
        end
      end
    end
  end
end
