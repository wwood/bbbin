#!/usr/bin/env ruby

require 'optparse'

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
      
      # return all genes with the species code e.g. hsap
      def genes_with_species_code(species_code)
        @genes.select do |g|
          g.match(/^#{species_code}\|/)
        end
      end
      
      def genes_without_species_codes
        @genes.collect do |g|
          split_species_and_id(g)[1]
        end
      end
      
      # Split up the species code and the gene ID
      #   pfal|PF10_0178 => ['pfal','PF10_0178']
      def split_species_and_id(gene_id)
        if matches = gene_id.match(/^([a-z]{3,4})\|(.+)$/)
          matches[1..2]
        else
          raise ParseException, "Couldn't parse OrthoMCL gene ID `#{gene_id}'"
        end
      end
      
      # Return an array of groups that are found by grepping the groups file
      # If the groups file path ends in .gz, then zcat is used to uncompress
      # the groups file before the grep
      def self.groups_by_grep(groups_file, grep_string)
        use_zcat = groups_file.match(/gz$/) #is this a gz file, or just a regular text file?
        if use_zcat
          cmd = "zcat '#{groups_file}' |egrep '#{grep_string}'"
        else
          cmd = "grep '#{grep_string}' '#{groups_file}'"
        end
        
        groups = []
        `#{cmd}`.strip.split(/\n/).each do |line|
          groups.push create_from_groups_file_line(line)
        end
        return groups
      end
    end
  end
end

if __FILE__ == $0
  options = {
  :orthomcl_groups_filename => '/home/ben/phd/data/orthomcl/v5/groups_OrthoMCL-5.txt.gz', #ben is the creator of this script, so gets dibs.
  :input_species_code => nil,
  :output_species_codes => nil,
  :inverse => false, #as per grep -v
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
    opts.on('-i','--input-species-code SPECIES_CODE','OrthoMCL species code e.g. hsap. Default nil, meaning orthomcl group identifiers are to be fed in. If dash, then species inputs are not checked and mapping is based purely on the orthomcl gene IDs') do |s|
      options[:input_species_code] = s
    end
    opts.on('-o','--output-species-codes SPECIES_CODE','output OrthoMCL species code(s), comma-separated. Default nil, meaning inputs are only mapped to OrthoMCL group IDs') do |s|
      options[:output_species_codes] = s.split(',')
    end
    opts.on('-v','--inverse','output OrthoMCL genes NOT matching the input species codes') do
      options[:inverse] = true
    end
  end
  o.parse!
  
  add_species_code = lambda do |species_code, gene_id|
    "#{species_code}|#{gene_id}"
  end
  
  # split on line breaks and whitespace
  matchers = ARGF.read.split(/\s+/)
  
  # Narrow the list of groups to iterate through for performance reasons
  grep_term = nil
  if options[:input_species_code] == '-'
    # Don't bother matching on species code, since we don't know it
    grep_term = matchers.join('|')
  elsif options[:input_species_code].nil?
    # matching on OrthoMCL group name so just match all of those
    grep_term = matchers.collect{|g| "#{g}:"}.join('|')
  else
    # match on species ID
    grep_term = '.' #could be faster likely, but can't be bothered testing it well enough
  end
  
  # Cache all the groups, grepping on lines with stuff in them (grep not really needed)
  groups = Bio::OrthoMCL::Group.groups_by_grep(
    options[:orthomcl_groups_filename], grep_term
  )
  $stderr.puts "Found #{groups.length} groups to further filter. Filtering now.."
  
  matchers.each do |gene_id|
      # if searching by genes and not by groups
      # multiple genes can be found, because the some gene names are the beginnings of others, e.g. MAL13P1.15 and MAL13P1.150
      # remove those genes that are longer
      found_group = false
      this_groups = groups.select do |g|
        selected = nil
        if options[:input_species_code] == '-'
          # Don't bother matching on species code, since we don't know it
          selected = g.genes_without_species_codes.include? gene_id
        elsif options[:input_species_code].nil?
          # matching on OrthoMCL group name, so we are already happy
          selected = (g.group_id == gene_id)
        else
          # match on species ID
          selected = g.genes.include? add_species_code.call(options[:input_species_code],gene_id)
        end
        found_group = g if selected
        selected
      end
      
      # Error checking
      if !found_group
        $stderr.puts "No groups found for input #{gene_id}, skipping"
        next
      elsif this_groups.length > 1
        $stderr.puts "More than expected (#{groups.length}) OrthoMCL groups found for input #{gene_id}, skipping"
        next
      end
      
      # output
      group = found_group
      to_output = []
      if options[:input_species_code]
        to_output.push gene_id
      end
      to_output.push group.group_id
      unless options[:output_species_codes].nil?
        if options[:inverse]
          to_output.push group.genes.reject{|d|
            options[:output_species_codes].include?(group.split_species_and_id(d)[0])
          }.join(',')
        else
          options[:output_species_codes].each do |code|
            to_output.push group.genes_with_species_code(code).join(',')
          end
        end
      end
      puts to_output.join("\t")
  end
end
