#!/usr/bin/env ruby

$:.unshift File.join(File.dirname(__FILE__))
require 'orthomcl_species_jumper'
require 'rubygems'
require 'reach'

# Return all those groups that have exactly one gene in a list of species

if __FILE__ == $0
  options = {
  :orthomcl_groups_filename => '/home/ben/phd/data/orthomcl/v5/groups_OrthoMCL-5.txt.gz', #ben is the creator of this script, so gets dibs.
  :input_species_codes => nil,
  }
  
  # parse options
  o = OptionParser.new do |opts|
    opts.banner = [
      'Usage: orthomcl_species_jumper.rb -g <orthomcl_groups_filename> -i <input_species_orthomcl_species_code> -o <output_species_orthomcl_species_code>',
      "\nA list of input IDs is piped in via STDIN. Requires grep, and zcat if the orthomcl file is gzipped\n",
    ]
    
    opts.on('-g','--orthomcl-gzip-groups-filename GZIP_FILENAME','Path to the OrthoMCL groups file (either gzipped or not - that is autodetected), downloadable from orthomcl.org') do |filename|
      options[:orthomcl_groups_filename] = filename
    end
    opts.on('-i','--input-species-codes SPECIES_CODES','REQUIRED. OrthoMCL species codes to inspect separated by commas e.g. \'hsap,pfal\'.') do |s|
      options[:input_species_codes] = s.split(',')
    end
  end
  o.parse!
  
  unless options[:input_species_codes].length > 0
    raise Exception, "Need to define the species codes that are being filtered for.., see usage"
  end
  
  # First, grep out those groups that have the first species in them to narrow down the search space
  groups = Bio::OrthoMCL::Group.groups_by_grep(
    options[:orthomcl_groups_filename],
    options[:input_species_codes][0]
  )
  
  # Reject those groups that don't fit the criteria
  groups.reject! do |g|
    reject = false
    options[:input_species_codes].each do |code|
      unless g.genes_with_species_code(code).length == 1
        reject = true
        break
      end
    end
    reject
  end
  
  puts groups.reach.group_id.join("\n")
end
