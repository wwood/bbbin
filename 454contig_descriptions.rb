#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'

$:.unshift File.join([ENV['HOME'], %w(git bioruby-newbler_outputs lib)].flatten)
require 'bio-newbler_outputs'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} -a alignment_info_file
      
      Takes an alignment info file, and outputs the contig mean and median coverages, as well as the length.\n\n"
      
    opts.on("-a", "--alignment-info 454ALIGNMENT_INFO.TSV", "From newbler [required]") do |f|
      options[:alignment_info_file] = f
    end
    opts.on("-f", "--contig-names-file CONTIG_LIST.TXT454ALIGNMENT_INFO.TXT", "Text file list of contig names, perhaps from 454UnscaffoldedContigs.rb [required]") do |f|
      options[:contig_list_file] = f
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ARGV.length != 0 or options[:alignment_info_file].nil? or options[:contig_list_file].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  
  
  
  target_contigs = []
  log.info "Caching contig names to output for in #{options[:contig_list_file]}..."
  target_contigs = File.open(options[:contig_list_file]).readlines.collect{|s| s.strip}
  log.info "Cached a list of #{target_contigs.length} contigs e.g. #{target_contigs[0]}"
  
  log.info "Parsing #{options[:alignment_info_file]}..."
  contig_hash = Bio::Newbler::AlignmentInfoFile.contig_hash(options[:alignment_info_file])
  log.info "Cached contig information for #{contig_hash.length} contigs"
  
  puts %w(contig length coverage median_coverage).join("\t")
  target_contigs.each do |contig_name|
    if contig_hash[contig_name]
      puts [
        contig_name,
        contig_hash[contig_name].length,
        contig_hash[contig_name].coverage,
        contig_hash[contig_name].median_coverage,
      ].join("\t")
    else
      log.warn "Unable to find contig #{contig_name} in the AlignmentInfo file, skipping."
    end
  end
  
  
  
  
  
end #end if running as a script