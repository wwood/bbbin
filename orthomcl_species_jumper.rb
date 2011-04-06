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
    end
  end
end

if __FILE__ == $0
  options = {
  :orthomcl_groups_filename => '/home/ben/phd/data/orthomcl/v4/groups_OrthoMCL-4.txt.gz', #ben is the creator of this script, so gets dibs.
  :input_species_code => nil,
  :output_species_codes => nil
  }
  
  # parse options
  o = OptionParser.new do |opts|
    opts.banner = [
      'Usage: orthomcl_species_jumper.rb -g <orthomcl_groups_filename> -i <input_species_orthomcl_species_code> -o <output_species_orthomcl_species_code>',
      "A list of input IDs is piped in via STDIN. Requires grep, and zcat if the orthomcl file is gzipped",
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
  end
  o.parse!
  
  add_species_code = lambda do |species_code, gene_id|
    "#{species_code}|#{gene_id}"
  end
  
  # split on line breaks and whitespace
  ARGF.each_line do |line|
    line.split(/\s+/).each do |gene_id|
      lines = []
      use_zcat = options[:orthomcl_groups_filename].match(/gz$/) #is this a gz file, or just a regular text file?
      if options[:input_species_code]
        # Are we grepping for species?
        to_grep = gene_id
        to_grep = add_species_code.call(options[:input_species_code],gene_id) unless options[:input_species_code] == '-'
        
        if use_zcat
          cmd = "zcat '#{options[:orthomcl_groups_filename]}' |grep '#{to_grep}'"
        else
          cmd = "grep '#{to_grep}' '#{options[:orthomcl_groups_filename]}'"
        end
        #$stderr.puts cmd
        lines = `#{cmd}`.strip.split(/\n/)
        #$stderr.puts lines
      else
        if use_zcat
          cmd = `zcat '#{options[:orthomcl_groups_filename]}' |grep '^#{gene_id}:'`
        else
          cmd = "grep '^#{gene_id}:' #{options[:orthomcl_groups_filename]}"
        end
        lines = `#{cmd}`.strip.split(/\n/)
      end
      
      # convert to parsed OrthoMCL groups 
      groups = lines.collect do |l|
        Bio::OrthoMCL::Group.create_from_groups_file_line l
      end
      
      # if searching by genes and not by groups
      # multiple genes can be found, because the some gene names are the beginnings of others, e.g. MAL13P1.15 and MAL13P1.150
      # remove those genes that are longer
      if options[:input_species_code]
        groups = groups.select do |g|
          if options[:input_species_code] == '-'
            # Don't bother matching on species code, since we don't know it
            g.genes_without_species_codes.include? gene_id
          else
            # match on species ID
            g.genes.include? add_species_code.call(options[:input_species_code],gene_id)
          end
        end
      end
      
      # Error checking
      if groups.length == 0
        $stderr.puts "No groups found for input #{gene_id}, skipping"
        next
      elsif groups.length > 1
        $stderr.puts "More than expected (#{groups.length}) OrthoMCL groups found for input #{gene_id}, skipping"
        next
      end
      
      # output
      group = groups[0]
      to_output = []
      if options[:input_species_code]
        to_output.push gene_id
      end
      to_output.push group.group_id
      unless options[:output_species_codes].nil?
        options[:output_species_codes].each do |code|
          to_output.push group.genes_with_species_code(code).join(',')
        end
      end
      puts to_output.join("\t")
    end
  end
end