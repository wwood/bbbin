#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-faster'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} -n <sample_name> <fastq1> [<fastq2> ..]
      
      Takes one or more fastq files, and transforms them:
      * all reads maximum 100bp long
      * all reads scores either quality 'T' or 'f', assumes input quality is either X or S.%
      * renames all fastq reads
      
      Prints out reads on STDOUT.
      \n\n"
      
    opts.on("-n", "--sample-name ARG", "Name of the sample [required]") do |arg|
      options[:sample_name] = arg
    end
    opts.on("-d", "--directory DIR", "Directory with fastq files in them (ending in .fastq) [not used by default]") do |arg|
      options[:directory] = arg
    end

    # logger options
    opts.separator "\nVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if options[:sample_name].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  counter = 0
  last_name = nil
  
  fastq_files = nil
  if options[:directory]
    fastq_files = Dir.glob(File.join(options[:directory],'*.fastq')) 
  else
    if ARGV.length > 0
      fastq_files = ARGV
    else
      fastq_files = [:stdin] #Bio::Faster#new accepts this argument
    end
  end
  log.info "Found #{fastq_files.length} fastq files to put together"
  
  fastq_files.each do |fastq_file|
    log.info "Working through #{fastq_file}"
    Bio::Faster.new(fastq_file).each_record(:quality => :raw) do |sequence_header, sequence, quality|
      splits = sequence_header.split(' ')
      splits2 = splits[0].split('/') #split on e.g. 1/2
      
      if last_name != splits2[0]
        counter += 1
      end
      last_name = splits2[0]
      
      
      new_name = counter.to_s+'/'+splits2[1]+' '+options[:sample_name]
      raise "Unexpected quality scores in #{quality}" unless quality.match(/^[XS]+$/)
      
      new_quality = quality.gsub('S','T').gsub('X','f')
      
      if sequence.length > 100
        sequence = sequence[0...100]
        new_quality = new_quality[0...100]
      end
      
      print '@'
      puts new_name
      puts sequence
      puts '+'
      puts new_quality
    end
  end
  
end #end if running as a script
