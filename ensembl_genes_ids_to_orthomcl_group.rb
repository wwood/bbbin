#!/usr/bin/env ruby


# A script that takes a list of ENSEMBL gene identifiers
# and maps them to OrthoMCL groups. Requires the orthomcl 
# groups file to be downloadable locally.
#
# Only tested on Ruby 1.9.2!

require 'rubygems'
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
  require 'optparse'
  
  # Default options
  options = {
    :orthomcl_groups_filename => '/home/ben/phd/data/orthomcl/v4/groups_OrthoMCL-4.txt.gz', #ben is the creator of this script, so gets dibs.
    :species => 'homo_sapiens', #because humans are closer to God.
    :ensembl_version => 56, # v3 and v4 of OrthoMCL use this (personal communication)
    :find_ensembl_names => false,
  }
  # Parse command line arguments
  o = OptionParser.new do |opts|
    opts.banner = [
      'Usage: ensembl_genes_ids_to_orthomcl_group.rb -o <orthomcl_gzip_groups_filename> [fasta_filename]',
      "fasta file can also be piped in on STDIN. Requires zcat, gzip",
    ]
    
    opts.on('-o','--orthomcl-gzip-groups-filename GZIP_FILENAME','Path to the OrthoMCL groups file (gzipped), downloadable from orthomcl.org') do |filename|
      options[:orthomcl_groups_filename] = filename
    end
    opts.on('-s','--species SPECIES','Connect to a this species\' Ensembl database. (Default \'homo_sapiens\')') do |s|
      options[:species] = s
    end
    opts.on('-v','--version VERSION','Connect to a specific version number of the Ensembl database. (Default 56 since that is what OrthoMCL v3 and v4 use)') do |s|
      options[:ensembl_version] = s
    end
    opts.on('-n','--find_names','Find genes by their name, rather than the default searching by Ensembl gene id (e.g. MFN2 instead of ENSG00000116688') do
      options[:find_ensembl_names] = true
    end
  end
  o.parse!
  
  include Ensembl::Core
  DBConnection.connect(options[:species],options[:ensembl_version].to_i)
  
  ARGF.each_line do |line|
    whitespaced = line.split(/[,\s]+/)
    whitespaced.each do |line|
      gene_id = line.strip
      $stderr.puts "Processing #{gene_id} .."
      
      g = nil
      if options[:find_ensembl_names]
        geneses = Gene.find_all_by_name gene_id
        if geneses.length > 1
          $stderr.puts "Found #{geneses.length} matches to '#{gene_id}', ignoring"
        else
          g = geneses[0] #the case of none found is handled later.
        end
      else
        g = Gene.find_all_by_stable_id gene_id
      end
      
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
        
        $stderr.puts "Found protein IDs with maximal length of CDS #{max_length_transcript}: #{max_translation_ids.join(', ')}"
        
        # Extract the group IDs for each Ensembl Protein ID
        hits = [] #the protein IDs that hit
        groups = [] #the orthomcl groups that are retrieved
        max_translation_ids.each do |t|
          $stderr.puts "Attempting to find translated protein ID #{t} from the OrthoMCL file"
          liner = `zcat '#{options[:orthomcl_groups_filename]}' |grep #{t}`.strip
          if liner.match(/\n/)
            $stderr.puts "Multiple groups found for protein ID #{t}! Whack."
          else
            unless liner == ''
              hits.push t
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
          hits.join(', '),
          groups[0].group_id
          ].join "\t"
        end
      end
    end
  end
end
