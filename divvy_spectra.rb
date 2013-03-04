#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'pp'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :contaminant_prefix => /^CNTM:/,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Takes a tab separated file containing a (possibly modified) output from a DTAselect run, and use some algorithm to divy up the spectra that match multiple peptides.\n\n"

  opts.on("--merge-proteins LIST_OF_IDENTIFIERS", "Provide a space/tab separated file where the identifiers on each row should be treated as one protein") do |file|
    options[:merge_proteins_file] = file
  end

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

  attr_accessor :sequence_count, :spectrum_count, :sequence_coverage, :length, :molwt, :pi, :validation_status, :descriptive_name

  attr_accessor :peptides

  def initialize
    @peptides = []
  end

  def unique_spectra
    return 0 if @peptides.nil? or @peptides.empty?
    num = @peptides.select{|pep| pep.parent_proteins.length == 1}.collect{|pep| pep.redundancy}.reduce(:+)
    num ||= 0
    return num
  end

  def estimated_spectral_count
    # How many unique spectra are there for each protein that shares a peptide with the current peptide
    return 0 if @peptides.nil? or @peptides.empty?
    peptide_shares = []
    peptides.each do |peptide|
      log.debug "Tallying peptide #{peptide.identifier}, which is has #{peptide.redundancy} spectra shared among #{peptide.parent_proteins.length} proteins"
      log.debug "These proteins have #{peptide.parent_proteins.collect{|pro| pro.unique_spectra}.inspect} unique spectra each"
      total_linked_unique_spectra = peptide.parent_proteins.collect{|pro| pro.unique_spectra}.reduce(:+)
      peptide_shares.push unique_spectra.to_f/total_linked_unique_spectra*peptide.redundancy
    end
    return peptide_shares.reduce(:+)
  end

  def log
    Bio::Log::LoggerPlus[LOG_NAME]
  end
end

class Peptide
  attr_accessor :identifier

  attr_accessor :xcorr, :deltcn, :obs_mono_mz, :cal_mono_mz, :total_intensity, :sp_rank, :sp_score, :ion_proportion, :redundancy, :sequence

  attr_accessor :unique

  attr_accessor :parent_proteins
  def initialize
    @parent_proteins = []
  end
end

# Read in merges, if required
mergers = {}
if options[:merge_proteins_file]
  File.open(options[:merge_proteins_file]).each_line do |line|
    splits = line.strip.split(/\s+/)
    primary_id = splits[0]
    splits.each_with_index do |s, i|
      next if i==0
      raise "This script can only handle two-way merging at the moment, sorry" if splits.length > 2
      raise "ID supposedly matches to multple identifiers: #{splits[1]}" if mergers[s] and mergers[s] != primary_id
      mergers[s] = primary_id
    end
  end

  log.info "Merging of identifiers setup for #{mergers.length} different instances, e.g. #{mergers.to_a[0][0]} => #{mergers.to_a[0][1]}"
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
    if splits[0] == 'Unique'
      reading_header = false
    end
    next
  end

  # OK, now we are reading the actual table, not the header
  if splits[0] != '' and splits[11].nil?
    log.debug "New protein now being parsed"
    # start a new protein
    current_protein = SelectedProtein.new
    ident = splits[0]
    current_protein.identifier = ident

    i = 1
    current_protein.sequence_count = splits[i].to_i; i+=1
    current_protein.spectrum_count = splits[i].to_i; i+=1
    current_protein.sequence_coverage = splits[i].to_f; i+=1
    current_protein.length = splits[i].to_i; i+=1
    current_protein.molwt = splits[i].to_f; i+=1
    current_protein.pi = splits[i].to_f; i+=1
    current_protein.validation_status = splits[i].to_f; i+=1
    current_protein.descriptive_name = splits[i]

    if proteins[ident]
      raise "Unexpectedly found the same protein identifier twice: #{ident}, from line #{line.chomp}"
    end
    proteins[ident] = current_protein




  elsif splits[1] == 'Proteins'
    # Done processing, except for the bits down the bottom which aren't parsed (yet)
    break



  # Now there is a
  else
    log.debug "New spectra now being parsed"
    # Record a spectra
    ident = splits[1]
    raise "Unexpected hits name `#{ident}', from line `#{line.chomp}'" unless ident.length > 10

    pep = hits[ident]
    if pep.nil?
      pep = Peptide.new
      pep.identifier = ident

      i = 2
      pep.xcorr = splits[i].to_f; i+= 1
      pep.deltcn = splits[i].to_f; i+= 1
      pep.obs_mono_mz = splits[i].to_f; i+= 1
      pep.cal_mono_mz = splits[i].to_f; i+= 1
      pep.total_intensity = splits[i].to_f; i+= 1
      pep.sp_rank = splits[i].to_f; i+= 1
      pep.sp_score = splits[i].to_f; i+= 1
      pep.ion_proportion = splits[i].to_f; i+= 1
      pep.redundancy = splits[i].to_i; i+= 1
      pep.sequence = splits[i]

      hits[ident] = pep
    end

    pep.parent_proteins.push current_protein
    current_protein.peptides.push pep
  end
end


# Merge proteins that are known duplicates if need be
mergers.each do |secondary_id, primary_id|
  log.debug "Merging proteins #{primary_id} and #{secondary_id}"
  if proteins[primary_id] and proteins[secondary_id]
    # Do the merge
    log.debug "Both are defined, so doing the complicated merge"

    # Invalidate some things about the primary ID because they are no longer valid
    current_protein = proteins[primary_id]
    current_protein.sequence_count = nil
    current_protein.sequence_coverage = nil
    current_protein.length = nil
    current_protein.molwt = nil
    current_protein.pi = nil
    current_protein.validation_status = nil
    # Keep the primary proteins' description, I reckon

    # When there is spectra that are in the secondary but not the primary, add them to the primary's repertoire.
    primary = proteins[primary_id]
    secondary = proteins[secondary_id]
    primary_peptide_names = primary.peptides.collect{|pep| pep.identifier}
    log.debug "Before transfer of the second protein's peptides, the primary proteins has #{primary.peptides.length} different peptides"
    log.debug "Parent protein IDs of primary peptides: #{primary.peptides.collect{|pep| pep.parent_proteins.collect{|pro| pro.identifier}}.inspect}"
    secondary.peptides.each do |sec_pep|
      unless primary_peptide_names.include?(sec_pep.identifier)
        primary.peptides.push sec_pep
        sec_pep.parent_proteins.push primary
      end
    end
    log.debug "After transfer of the second protein's peptides, the primary proteins has #{primary.peptides.length} different peptides"
    log.debug "Parent protein IDs of primary peptides: #{primary.peptides.collect{|pep| pep.parent_proteins.collect{|pro| pro.identifier}}.inspect}"
    # Remove references second protein from the second peptides
    secondary.peptides.each do |pep|
      pep.parent_proteins.reject!{|pro| pro==secondary}
    end
    log.debug "Parent protein IDs of primary peptides: #{primary.peptides.collect{|pep| pep.parent_proteins.collect{|pro| pro.identifier}}.inspect}"
    # Remove the secondary peptide from the list of peptides
    proteins.delete secondary_id


  elsif proteins[secondary_id]
    raise "You've reached a place in the code that is implemented but untested"
    # Rename the secondary as the primary
    sec = proteins[secondary_id]
    proteins[primary_id] = sec
    proteins.delete secondary_id
    sec.identifier = primary_id
    # The peptide objects should have the correct parent proteins because it is all references

  end #The other two cases do not require any intervention,
end


# Total spectra shouldn't count contaminants, but shared spectra should still be divvied up with
total_contaminating_spectra = proteins.select{|ident, protein| ident.match(options[:contaminant_prefix])}.collect{|i, pro| pro.estimated_spectral_count}.reduce(:+)
total_contaminating_spectra ||= 0

total_spectra = hits.collect{|i,pep| pep.redundancy}.reduce(:+) - total_contaminating_spectra
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
  next if protein_id.match(options[:contaminant_prefix]) #Don't print contaminants

  log.debug "Now printing protein #{protein_id}, which has #{protein.peptides.length} associated peptides"
  puts [
    protein_id,
    protein.unique_spectra,
    protein.estimated_spectral_count,
    protein.estimated_spectral_count.to_f / total_spectra,
    protein.descriptive_name,
  ].join "\t"
end






