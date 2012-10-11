#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

$:.unshift File.join([ENV['HOME'], %w(git bioruby-newbler_outputs lib)].flatten)
require 'bio-newbler_outputs'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME}
      
      Returns a list of contigs with some minimum coverage, and that are not scaffolded at all by newbler.\n\n"
    # Example option
    opts.on("-s", "--scaffolds-txt 454SCAFFOLDS.TXT", "From newbler [required]") do |f|
      options[:scaffolds_file] = f
    end
    opts.on("-m", "--min-coverage MIN_COVERAGE", "From newbler [required]") do |f|
      options[:min_coverage] = f.to_f
      raise if options[:min_coverage] < 0 
    end
    opts.on("-c", "--contig-graph 454CONTIGGRAPH.TXT", "From newbler [required unless -a/--alignment-info is defined]") do |f|
      options[:contigs_graph_file] = f
    end
    opts.on("-a", "--alignment-info 454ALIGNMENT_INFO.TXT", "From newbler [required unless -c/--contig-graph is defined]") do |f|
      options[:alignment_info_file] = f
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") do |q|
      Bio::Log::CLI.trace('error')
    end
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") do | name |
      options[:logger] = name
    end
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG") do | s |
      Bio::Log::CLI.trace(s)
    end
  end
  o.parse!
  if ARGV.length != 0
    $stderr.puts o
    exit 1
  end
  # Setup logging
  Bio::Log::CLI.logger(options[:logger]) #bio-logger defaults to STDERR not STDOUT, I disagree
  log = Bio::Log::LoggerPlus.new(LOG_NAME)
  Bio::Log::CLI.configure(LOG_NAME)
  
  
  # Read in the scaffolds file, these are the ones that should not be printed out
  scaffolded_contigs = {}
  log.info "Parsing scaffolds file #{options[:scaffolds_file]}"
  CSV.foreach(options[:scaffolds_file], :col_sep => "\t") do |row|
    scaffolded_contigs[row[5]] = true if row[4]=='W'
  end
  log.info "Found #{scaffolded_contigs.length} contigs in scaffolds"
  
  low_coverage_count = 0
  ignored_contigs_count = 0
  printed_contigs_count = 0
  
  process_contig = lambda do |contig_name, coverage|
    if coverage < options[:min_coverage]
      low_coverage_count += 1
    elsif scaffolded_contigs[contig_name]
      ignored_contigs_count += 1
    else
      printed_contigs_count += 1
      puts contig_name
    end
  end
  
  # Either use the contigs graph file or the alignment info file
  if options[:contigs_graph_file]
    log.info "Parsing contig graph file #{options[:contigs_graph_file]}"
    CSV.foreach(options[:contigs_graph_file], :col_sep => "\t") do |row|
      break if row[0] == 'C'
      raise unless row[0].to_i > 0 #indicates unexpected file structure
      
      coverage = row[3].to_f
      contig_name = row[1]
      process_contig.call(contig_name, coverage)
    end
    
  else
    last_contig = nil
    count = 0
    total_length = 0
    log.info "Parsing alignment info file #{options[:alignment_info_file]}"
    hash = Bio::Newbler::AlignmentInfoFile.contig_hash(options[:alignment_info_file])
    hash.each do |name, contig|
      process_contig.call(name, contig.coverage)
    end
  end
  
  log.info "Ignored due to low coverage: #{low_coverage_count}"
  log.info "Ignored due to being in a scaffold: #{ignored_contigs_count}"
  log.info "Printed contigs: #{printed_contigs_count}"

end















