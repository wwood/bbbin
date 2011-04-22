#!/usr/bin/env ruby

# Interrogate a EuPathDB gene information table file, extracting information
# about a particular genes (or perhaps in the future, a set of genes)

require 'optparse'
require 'rubygems'
require 'eupathdb_gene_information_table' #from the reubypathdb rubygem

if $0 == __FILE__
  options = {
    :quit_after_first => true,
    :gene_information_table_path => '/home/ben/phd/data/Toxoplasma gondii/genome/ToxoDB/6.0/TgondiiME49Gene_ToxoDB-6.0.txt',
    :grep_hack => false,
  }
  o = OptionParser.new do |opts|
    opts.banner = [
      'Usage: eupathdb_gene_information_file_extractor.rb <eupathdb_gene_identifier> <key>'
    ]
    opts.on("-g", "--grep-hack", "Invoke the grep hack to make gene information file extraction faster") do
      options[:grep_hack] = true
    end
    if ARGV.length != 2
      $stderr.puts opts.banner
      exit
    end
  end
  o.parse!
  
  wanted_gene_id = ARGV[0]
  info = ARGV[1]
  
  gene = EuPathDBGeneInformationFileExtractor.new(options[:gene_information_table_path]).extract_gene_info(wanted_gene_id, 1000)
  
  if gene.get_info(info).nil?
    table = gene.get_table(info)
    if table.nil?
      $stderr.puts "No attribute `#{info} found!`"
    else
      p table
    end
  else
    puts gene.get_table(info)
  end
end
