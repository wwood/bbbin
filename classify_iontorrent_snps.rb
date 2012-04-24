#!/usr/bin/env ruby

# Takes in a reference fasta file and a file of SNPs from dnadiff output e.g. out.snps

# Outputs information about the snps

require 'bio'
require 'csv'
require 'optparse'

class SNP
  attr_reader :ref, :aligned
  attr_accessor :position, :ref_name
  
  def ref=(ref)
    if ref == '.'
      @ref = nil
    else
      @ref = ref
    end
  end
  
  def aligned=(aligned)
    if aligned == '.'
      @aligned = nil
    else
      @aligned = aligned
    end
  end
  
  def type
    return 'insertion' if @ref.nil? and !@aligned.nil?
    return 'deletion' if deletion?
    return 'substitution' if !@ref.nil? and !@aligned.nil?
    p inspect
    raise
  end
  
  def deletion?
    !@ref.nil? and @aligned.nil?
  end
  
  def homopolymer_length(seq)
    count = 1 #count itself
    
    # How many out in front?
    i = 1
    next_base = seq[@position+i]
    while next_base == ref
      count += 1
      i += 1
      next_base = seq[@position+i]
    end

    # How many out back?
    i = -1
    next_base = seq[@position+i]
    while next_base == ref
      count += 1
      i -= 1
      next_base = seq[@position+i]
    end
    
    return count
  end
end






# Parse options
OVERALL = 'overall'
DELETION = 'deletion'
options = {
  :operation => OVERALL,
}
o = OptionParser.new do |opts|
  opts.banner = [
    'Usage: classify_iontorrent_snps.rb <reference fasta file> <dnadiff out.snps>',
    "\tTakes in a reference fasta file and a file of SNPs from dnadiff output e.g. out.snps\n"+
    "\tOutputs information about the snps\n"
  ]
  opts.on("-o", "--overall", "Overall classification of SNPs IN/DEL/SUB breakdown") do |f|
    options[:operation] = OVERALL
  end
  opts.on('-d', "--deletion", "Investigate the deletions specifically") do |e|
    options[:operation] = DELETION
  end
end
o.parse!
if ARGV.length != 2
  $stderr.puts o.banner
  exit
end



# Read in the SNPs
snps = []
CSV.foreach(ARGV[1],:col_sep => "\t") do |row|
  snp = SNP.new
  snp.ref = row[1]
  snp.aligned = row[2]
  snp.position = row[0].to_i
  snp.ref_name = row[10]
  snps.push snp
end

types = {}
snps.each do |snp|
  types[snp.type] ||= 0
  types[snp.type] += 1
end



# Classify each of SNPs
if options[:operation] == OVERALL
  types.each do |type, count|
    puts [type,count].join("\t")
  end



# Print out more specifics on the deletions
elsif options[:operation] == DELETION
  # Read in the fasta sequences
  seqs = {}
  Bio::FlatFile.foreach(ARGV[0]){|s| seqs[s.definition.split(' ')[0]] = s}
  
  #Print out the surrounding sequences around the deletions
  snps.each do |snp|
    next unless snp.deletion?
    
    ref = seqs[snp.ref_name]
    surround4 = ref.seq[snp.position-2..snp.position+1]
    homopolymer_length = snp.homopolymer_length(ref.seq) 
    
    puts [
      snp.ref_name,
      snp.position,
      snp.ref,
      surround4,
      homopolymer_length,
    ].join("\t")
  end
end






