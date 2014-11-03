#!/usr/bin/env ruby

require 'rubygems'
require 'bio'
require 'optparse'
require 'csv'
require 'bio-logger'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
options = {
  :email => nil,
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "Usage: entrez_proteins_from_nucleotide_ids.rb --email/-e <email_for_ncbi> <nucleotide_accessions_file>\n\nScript to extract a FASTA file of protein sequences associated with input GenBank nucleotide accession numbers. The input file can either be a one column input, or 2 column tab separated with the second column specifying a lable to prepend each translated sequence with."

  opts.on('-e', "--email EMAIL", "Email address to provide to NCBI") do |v|
    options[:email] = v
  end

  opts.on('-g', "--genbank GENBANK_FILE", "Skip entrez, instead provide a GenBank file to extract the CDS from directly. No input file of accessions is used either.") do |v|
    options[:genbank] = v
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end.parse!

unless options[:email]
  $stderr.puts "ERROR: email address required\n\n"
  $stderr.puts o.banner
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


process_genbank = lambda do |gb_string, prepending, encompassing_name|
  gb = Bio::GenBank.new(gb_string)
  gb.each_cds do |cds|
    ha = cds.to_hash

    # skip pseudogenes
    if ha['protein_id'].nil?
      log.warn "Unable to find any linked protein identifiers for a CDS of #{encompassing_name}, skipping!"
      log.warn "CDS was: #{cds.inspect}"
      next
    end
    if ha['product'].nil?
      log.warn "Unable to find any linked product identifiers for a CDS of #{encompassing_name}, skipping!"
      log.warn "CDS was: #{cds.inspect}"
      next
    end

    # define the name as >protein_id description
    name = ""
    unless prepending.nil?
      name = "#{prepending}|"
    end
    name = "#{name}#{ha['protein_id'][0]} #{ha['product'][0]}"

    # Get the protein translation
    translates = ha['translation']
    if translates.length == 1
      log.debug "Printing #{encompassing_name} CDS #{name}"
      puts ">#{name}"
      puts translates[0]
    else
      $stderr.puts "INFO Skipping #{name}, since #{translates.length} 'translation's were found, expected 1"
    end
  end
end



if options[:genbank]
  process_genbank.call(File.open(options[:genbank]).read, 'abc', 'def')

else
  CSV.foreach(ARGV[0],:col_sep => "\t") do |row|
    line = row[0].strip
    next if line.length == 0
    prepend = row[1]

    $stderr.puts "Fetching #{prepend}|#{line} from GenBank.."
    gb = Bio::NCBI::REST.efetch(line, {"db"=>"nuccore","rettype"=>"gb", "email"=>options[:email]})

    process_genbank.call(gb, prepend, line)
  end
end

