#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :num_processes => 10,
    :fwd_read_grep => ' 1:',
    :rev_read_grep => ' 2:',
    :file_type => 'fastq',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <fastq_file1> [<fastq_file2> ..]

      Take one or more fastq files that contain interleaved contents i.e. read1_1, read1_2, read2_1, read2_2, .. and parse them into output files.

      CHECK before you run this script that it has the expected pairing structure
      e.g.
      \n\n"

    opts.on('--slash','read headers are formatted like "@FCC0WM1ACXX:2:1101:2167:2180#GTCCAGAA/1", so split on the /1 and /2') do
      options[:fwd_read_grep] = '\/1'
      options[:rev_read_grep] = '\/2'
    end
    # this fails because the @ is not replaced with a @
    #opts.on('--fasta','file is fasta, not fastq [default: false]') do
    #  options[:file_type] = 'fasta'
    #end
    opts.on('--threads NUM',Integer,'file is fasta, not fastq [default: #{options[:num_processes]}]') do |i|
      options[:num_processes] = i
    end

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
    if options[:file_type] == 'fasta'
      fq1 = "#{base}_1.fa.gz"
      fq2 = "#{base}_2.fa.gz"
    end

    [fq1, fq2].each do |fq|
      if File.exists?(fq)
        raise "Cowardly failing because #{fq} already exists!"
      end
    end

    num_following_lines = {'fastq' => '3', 'fasta' => '1'}[options[:file_type]]
    log.debug("Using #{num_following_lines} lines following")

    log.info "Creating #{fq1} .."
    `zcat '#{fastq_file}' | grep -A#{num_following_lines} '#{options[:fwd_read_grep]}' |grep -v '^--$' |pigz -p #{options[:num_processes]} >#{fq1}`
    fq1_size = File.size fq1
    log.info "#{fq1} size #{fq1_size}"

    log.info "Creating #{fq2} .."
    `zcat '#{fastq_file}' | grep -A#{num_following_lines} '#{options[:rev_read_grep]}' |grep -v '^--$' |pigz -p #{options[:num_processes]} >#{fq2}`
    fq2_size = File.size fq2
    log.info "#{fq2} size #{fq2_size}"
  end
end #end if running as a script
