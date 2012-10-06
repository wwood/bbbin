#!/usr/bin/env ruby

# Takes in a reference fasta file and a file of SNPs from dnadiff output e.g. out.snps

# Outputs information about the snps

require 'bio'
require 'csv'
require 'optparse'
require 'bio-assembly'
require 'bio-logger'

class SNP
  attr_reader :ref, :aligned
  attr_accessor :position, :ref_name, :contig_name, :contig_position
  
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




SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
# Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
Bio::Log::CLI.logger('stderr'); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


# Parse options
OVERALL = 'overall'
DELETION = 'deletion'
DELETION_STRANDEDNESS = 'deletion_strandedness'

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
  opts.on("--deletion-strandedness", "Investigate the strand bias at deletion positions specifically. Requires --pileup") do |e|
    options[:operation] = DELETION_STRANDEDNESS
  end
  # opts.on("--ace ACE_FILE", "ACE file representing the assembly") do |e|
    # options[:ace_file] = e
  # end
  opts.on("--pileup PILEUP_FILE", "pileup file representing the assembly") do |e|
    options[:pileup_file] = e
  end
end
o.parse!
if ARGV.length != 2
  $stderr.puts o.banner
  exit
end



# Read in the SNPs
log.info "Reading SNPs"
snps = []
CSV.foreach(ARGV[1],:col_sep => "\t") do |row|
  snp = SNP.new
  
  # [P1]  [SUB] [SUB] [P2]  [BUFF]  [DIST]  [LEN R] [LEN Q] [FRM] [TAGS]
  snp.ref = row[1]
  snp.aligned = row[2]
  snp.position = row[0].to_i
  snp.contig_position = row[3].to_i
  snp.ref_name = row[10]
  snp.contig_name = row[11]
  snps.push snp
end

types = {}
snps.each do |snp|
  types[snp.type] ||= 0
  types[snp.type] += 1
end

log.info "Read in #{snps.length} SNPs"


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
  
  
elsif options[:operation] = DELETION_STRANDEDNESS
  asm = Bio::Assembly.open(options[:ace_file], :ace)
  
  # Get a hash of the contig names
  log.info "Reading in assemblies for each contig"
  contigs = {}
  asm.each_contig do |contig|
    raise unless contigs[contig.name].nil?
    contigs[contig.name] = contig
  end
  log.info "Read in assembly info for #{contigs.length} different contigs"
  
  # Read in the fasta sequences
  seqs = {}
  Bio::FlatFile.foreach(ARGV[0]){|s| seqs[s.definition.split(' ')[0]] = s}
  
  gag_type_contexts = %w(GAAG CTTC  AGGC GCCT  GCCG CGGC  GCCA TGGC)
  
  #Print out the surrounding sequences around the deletions
  snps.each do |snp|
    next unless snp.deletion?
    
    ref = seqs[snp.ref_name]
    surround4 = ref.seq[snp.position-2..snp.position+1]
    
    # Disregard non-gag types
    next unless gag_type_contexts.include?(surround4)
    
    # How many reads are 1) in each direction and 2) have a deletion at that position?
    reads_here = contigs[snp.contig_name].find_reads_in_range(snp.contig_position, snp.contig_position+1)
    #p reads_here.length
    #reads_here.each do |r| p r.name; p r.orientation; end
    forwards = reads_here.select{|read| read.orientation == 'C'}.length
    backwards = reads_here.select{|read| read.orientation == 'U'}.length
    
    
    puts [
      snp.ref_name,
      snp.contig_name,
      snp.position,
      snp.contig_position,
      snp.ref,
      surround4,
      forwards,
      backwards,
      reads_here.length
    ].join("\t")
  end
end






