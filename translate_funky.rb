#!/usr/bin/env ruby

# Take in a fasta file of nucleotide sequences corresponding to a collection of proteins, and then translate them. However, the translation table is modifyable on the command line

require 'bio'
require 'optparse'
require 'bio-logger'

translation_table = Bio::CodonTable.copy(1)
table_modifications = {}

options = {
  :logger => 'stderr',
  :trace => 'info',
}

o = OptionParser.new do |opts|
  opts.banner = "translate_funky.rb [-c <codon_modifications>] <transcript_fasta_file>"

  opts.on('-c', "--codon-modifications CODON_MODS", "Modify the codon table e.g. 'ATG:L' to change ATG from M to L") do |v|
    splits = v.split(',')
    splits.each do |split|
      if matches = split.match(/^([ATGC]{3}):([A-Za-z])$/)
        table_modifications[matches[1].downcase] = matches[2]
      else
        raise Exception, "Unexpected table codon modification format '#{split}', expected something like 'ATG:L'"
      end
    end
  end
  
  opts.on(nil, "--debug", "Run very verbosely") do |v|
    options[:trace] = 'funky-translate:debug'
  end
end.parse!


Bio::Log::CLI.trace(options[:trace])
Bio::Log::CLI.logger(options[:logger]) #defaults to STDERR not STDOUT
log = Bio::Log::LoggerPlus.new('funky-translate')
Bio::Log::CLI.configure('funky-translate')


table_modifications.each do |codon, aa|
  translation_table[codon] = aa
end


log.debug "Using CodonTable #{translation_table.inspect}"

Bio::FlatFile.foreach(ARGF) do |entry|
  puts ">#{entry.definition}"
  puts Bio::Sequence::NA.new(entry.seq).translate(1,translation_table)
end
