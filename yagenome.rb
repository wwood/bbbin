#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
#require '/home/ben/git/bioruby-hmmer3_report/lib/bio/appl/hmmer/hmmer3/tabular_report'
require '/srv/whitlam/home/users/uqbwoodc/git/bioruby-hmmer3_report/lib/bio/appl/hmmer/hmmer3/tabular_report'
#require 'bio-hmmer3_report'
require 'tmpdir'
require 'bio-hmmer_model'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} -i <proteome.faa> --hmm <hmm>
      
      Takes a fasta file representing a single genome's protein complement and a HMM file. 
      
      Outputs a fasta file of the alignment of the best matching protein from that proteome
      against the HMM, or gaps if no hit was found. \n\n"
      
    opts.on("-i", "--input-fasta PROTEOME", "input proteome [required]") do |arg|
      options[:input_fasta_file] = arg
    end
    opts.on("--hmm HMMFILE", "input HMM [required]") do |arg|
      options[:input_hmm_file] = arg
    end

    # logger options
    opts.separator "\n\tVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length != 0 or options[:input_fasta_file].nil? or options[:input_hmm_file].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  
  fasta_path = File.absolute_path options[:input_fasta_file]
  hmm_path = File.absolute_path options[:input_hmm_file]
  
  hmm_length = Bio::HMMER::Model.parse(File.open(hmm_path)).leng
  
  Dir.mktmpdir do |dir|
    Dir.chdir dir
    
    # Run hmmsearch to find the single-copy gene
    `hmmsearch --tblout hmmsearch.tblout '#{hmm_path}' #{fasta_path} >/dev/null`
    #`hmmsearch --tblout hmmsearch.tblout --incT 999999 '#{hmm_path}' #{fasta_path} -o /dev/null`
    
    # Parse the output file, retrieving the best hit ID (if one exists)
    if File.exist?('hmmsearch.tblout')
      result = Bio::HMMER::HMMER3::TabularReport.new(File.open('hmmsearch.tblout').read)
      best_hit = result.hits[0]
      if best_hit.nil?
        raise "No best hit.."
      end
      best_hit_id = best_hit.target_name
      #log.debug "Found the best hit called #{best_hit_id}"
      
      # Extract the best hit from the input fasta file
      `extract_from_fasta.pl -F '#{best_hit_id}' #{fasta_path} >best_hit.faa`
    
      # Run HMM align against the genome
      `hmmalign --allcol --trim -o best_hitVhmm.sto #{hmm_path} best_hit.faa`
      
      # Convert the stockholm output file to fasta, puts to stdout
      `seqmagick.py convert best_hitVhmm.sto best_hitVhmm.fa`
      puts `cat best_hitVhmm.fa |sed 's/[a-z]//g'`
    else
      puts '>no_hit'
      puts '-'*hmm_length
    end
  end

end #end if running as a script
