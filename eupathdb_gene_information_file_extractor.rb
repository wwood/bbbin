#!/usr/bin/env ruby

# Interrogate a EuPathDB gene information table file, extracting information
# about a particular genes (or perhaps in the future, a set of genes)

require 'optparse'
require 'rubygems'
require 'eupathdb_gene_information_table' #from the reubypathdb rubygem

# A class for extracting gene info from a gene info file
class EuPathDBGeneInformationFileExtractor
  def initialize(filename = nil)
    @filename = filename
  end

  # Returns a EuPathDBGeneInformation object corresponding to the wanted key. If
  # there are multiple in the file, only the first is returned.
  def extract_gene_info(wanted_gene_id)
    EuPathDBGeneInformationTable.new(File.open(@filename)).each do |gene|
      return gene if wanted_gene_id == gene.get_info('Gene Id')
    end
    return nil
  end
end

if $0 == __FILE__
  options = {
    :quit_after_first => true,
    :gene_information_table_path => '/home/ben/phd/data/Toxoplasma gondii/genome/ToxoDB/6.0/TgondiiME49Gene_ToxoDB-6.0.txt'
  }
  o = OptionParser.new do |opts|
    opts.banner = [
      'Usage: eupathdb_gene_information_file_extractor.rb <eupathdb_gene_identifier> <key>'
    ]
    if ARGV.length != 2
      $stderr.puts opts.banner
      exit
    end
  end
  o.parse!

  wanted_gene_id = ARGV[0]
  info = ARGV[1]

  gene = EuPathDBGeneInformationFileExtractor.new(options[:gene_information_table_path]).extract_gene_info(wanted_gene_id)

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
