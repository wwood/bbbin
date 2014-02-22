#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-samtools'
require 'bio'
require 'progressbar'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :window_length => 100,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <arguments>

      Take an (indexed) bam file and and plot the coverages of windows from the base sequence, recording the GC content of the window\n\n"

    opts.separator "\nRequired:\n\n"
    opts.on('-b',"--bam FILENAME", "BAM file to get input from. Must be indexed [required]") do |arg|
      options[:bam_file] = arg
    end
    opts.on('-f',"--fasta FILENAME", "Fasta file of the reference. Must be indexed [required]") do |arg|
      options[:fasta_file] = arg
    end

    opts.separator "\nOptional:\n\n"
    opts.on("--window-size LENGTH", "Length of window to cover [default: #{options[:window_length]}]") do |arg|
      options[:window_length] = arg.to_i
      raise unless options[:window_length]>0
    end

    # logger options
    opts.separator "\nVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length != 0 or options[:bam_file].nil? or options[:fasta_file].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


  # For each window of the reference sequence
  sam = Bio::DB::Sam.new(
    :bam => options[:bam_file],
    :fasta => options[:fasta_file],
  )
  sam.load_index
  sam.load_reference
  sam.open
  log.debug "Loaded SAM object: #{sam.to_s}"

  puts %w(
    sequence
    start
    stop
    gc
    coverage
    ).join "\t"

  Bio::FlatFile.foreach(options[:fasta_file]) do |s|
    identifier = s.definition.split(/\s/)[0]
    length = s.seq.length
    progress = ProgressBar.new(identifier, length/options[:window_length])

    1.step(length, options[:window_length]) do |i|
      progress.inc
      stop = i+options[:window_length]-1
      #p [identifier, i, stop]
      average = sam.average_coverage(identifier, i, stop)
      #p average
      seq_chunk = s.seq.seq[i..stop]

      puts [
        identifier,
        i,
        stop,
        Bio::Sequence::NA.new(seq_chunk.seq).gc_percent,
        average,
      ].join("\t")
    end
    progress.finish
  end
end #end if running as a script





