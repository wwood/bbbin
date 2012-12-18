#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <fastq_file1> [<fastq_file2> ..]
      
      Take one or more fastq files that contain interleaved contents i.e. read1_1, read1_2, read2_1, read2_2, .. and parse them into output files.
      
      CHECK before you run this script that it has the expected pairing structure
      e.g. 
      \n\n"

    # logger options
    opts.separator "\nVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length == 0
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  ARGV.each do |fastq_file|
    raise "not fastq file found: #{fastq_file}" unless File.exists?(fastq_file)
    base = File.basename fastq_file
    fq1 = "#{base}_1.fq.gz"
    fq2 = "#{base}_2.fq.gz"
    
    log.info "Creating #{fq1} .."
    `zcat '#{fastq_file}' | grep -A3 ' 1:' |grep -v '^--$' |pigz -p 5 >#{fq1}`
    fq1_size = File.size fq1
    log.info "#{fq1} size #{fq1_size}"
    
    log.info "Creating #{fq2} .."
    `zcat '#{fastq_file}' | grep -A3 ' 2:' |grep -v '^--$' |pigz -p 5 >#{fq2}`
    fq2_size = File.size fq2
    log.info "#{fq2} size #{fq2_size}"
  end
end #end if running as a script