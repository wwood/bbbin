#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'pp'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Takes a tab separated file containing a (possibly modified) output from a DTAselect run, and use some algorithm to divy up the spectra that match multiple peptides.\n\n"

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length > 1
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

class SelectedProtein
  attr_accessor :identifier

  attr_accessor :sequence_count, :spectrum_count, :sequence_coverage, :length, :molwt, :pi, :descriptive_name

  attr_accessor :peptides

  def initialize
    @peptides = []
  end

  def unique_spectra
    return 0 if @peptides.nil? or @peptides.empty?
    num = @peptides.select{|pep| pep.parent_proteins.length == 1}.collect{|pep| pep.number}.reduce(:+)
    num ||= 0
    return num
  end

  def estimated_spectral_count
    # How many unique spectra are there for each protein that shares a peptide with the current peptide
    return 0 if @peptides.nil? or @peptides.empty?
    peptide_shares = []
    peptides.each do |peptide|
      log.debug "Tallying peptide #{peptide.identifier}, which is has #{peptide.number} spectra shared among #{peptide.parent_proteins.length} proteins"
      log.debug "These proteins have #{peptide.parent_proteins.collect{|pro| pro.unique_spectra}.inspect} unique spectra each"
      total_linked_unique_spectra = peptide.parent_proteins.collect{|pro| pro.unique_spectra}.reduce(:+)
      peptide_shares.push unique_spectra.to_f/total_linked_unique_spectra*peptide.number
    end
    return peptide_shares.reduce(:+)
  end

  def log
    Bio::Log::LoggerPlus[LOG_NAME]
  end
end

class Peptide
  attr_accessor :identifier

  attr_accessor :xcorr, :deltcn, :obs_mono_mz, :cal_mono_mz, :ppm, :delta_amu, :ion_percent, :number, :sequence

  attr_accessor :parent_proteins
  def initialize
    @parent_proteins = []
  end
end

# Hashes of identifiers to objects
proteins = {}
hits = {}

# Read in the tab separated file
reading_header = true
current_protein = nil

ARGF.each_line do |line|
  splits = line.chomp.split("\t")
  log.debug "Parsing line `#{line.chomp}'"

  if reading_header
    log.debug "reading header"
    if splits[1] == 'Locus'
      reading_header = false
    end
    next
  end

  # OK, now we are reading the actual table, not the header
  if splits[1] == 'U'
    log.debug "New protein now being parsed"
    # start a new protein
    current_protein = SelectedProtein.new
    ident = splits[2]
    current_protein.identifier = ident

    i = 3
    current_protein.sequence_count = splits[i].to_i; i+=1
    current_protein.spectrum_count = splits[i].to_i; i+=1
    current_protein.sequence_coverage = splits[i].to_f; i+=1
    current_protein.length = splits[i].to_i; i+=1
    current_protein.molwt = splits[i].to_f; i+=1
    current_protein.pi = splits[i].to_f; i+=1
    current_protein.descriptive_name = splits[i]

    if proteins[ident]
      raise "Unexpectedly found the same protein identifier twice: #{ident}"
    end
    proteins[ident] = current_protein



  elsif !splits[2].nil? and splits[2].match(/^Similarities: /)
    # do nothing, this information is encoded elsewhere already



#  elsif splits[3] == 'Proteins'
#    # Done processing, except for the bits down the bottom which aren't included
#    finish_protein.call(current_protein)
#    break

  elsif splits[2] == 'Filename'
    log.debug "Current line is all headers"
    next #Current line is all headers


  # Have to test for both columns 2 and 3 here because sometime the similarity bits go over multiple lines
  elsif !splits[2].nil? and !splits[3].nil?
    log.debug "New spectra now being parsed"
    # Record a spectra
    ident = splits[2]
    raise "Unexpected hits name `#{ident}', from line `#{line.chomp}'" unless ident.length > 10

    pep = hits[ident]
    if pep.nil?
      pep = Peptide.new
      pep.identifier = ident

      i = 3
      pep.xcorr = splits[i].to_f; i+= 1
      pep.deltcn = splits[i].to_f; i+= 1
      pep.obs_mono_mz = splits[i].to_f; i+= 1
      pep.cal_mono_mz = splits[i].to_f; i+= 1
      pep.ppm = splits[i].to_f; i+= 1
      pep.delta_amu = splits[i].to_f; i+= 1
      pep.ion_percent = splits[i].to_f; i+= 1
      pep.number = splits[i].to_i; i+= 1
      pep.sequence = splits[i]; i+= 1
      #what is the last column? I'm asking EK now

      hits[ident] = pep
    end

    pep.parent_proteins.push current_protein
    current_protein.peptides.push pep
  end
end

total_spectra = hits.collect{|i,pep| pep.number}.reduce(:+)
log.info "Parsed in #{proteins.length} proteins and #{hits.length} peptides, and #{total_spectra} spectra"



# OK, finished parsing the file. Now output the score for each protein
puts [
  'ID',
  'Unique spectra',
  'Estimated total spectra',
  'Normalised spectral count',
  'Description',
].join "\t"
proteins.each do |protein_id, protein|
  log.debug "Now printing protein #{protein_id}, which has #{protein.peptides.length} associated peptides"
  puts [
    protein_id,
    protein.unique_spectra,
    protein.estimated_spectral_count,
    protein.estimated_spectral_count.to_f / total_spectra,
    protein.descriptive_name,
  ].join "\t"
end






