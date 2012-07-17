#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'progressbar'

$LOAD_PATH.unshift('/home/ben/git/bioruby-sra_fastq_dumper/lib')
require 'bio-sra_fastq_dumper'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stdout',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>
    
    Take a directory of .lite.sra files, and correct them using acacia, outputing the results in a second directory filled with directories (one for each sra lite file) \n\n"

    opts.on("-i", "--input-directory DIRECTORY", "a directory full of SRA .lite.sra files [required]") do |arg|
      options[:in_directory] = File.absolute_path(arg)
    end
    opts.on("-o", "--output-directory DIRECTORY", "a directory to put output files in [required]") do |arg|
      options[:out_directory] = File.absolute_path(arg)
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
if ARGV.length != 0 or options[:in_directory].nil? or options[:out_directory].nil?
  $stderr.puts o
  exit 1
end
# Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


# Collect a list of the input sra lite files to process, advise
log.debug "Inspecting #{options[:in_directory]}"
sralites = Dir.entries(options[:in_directory]).select{|s| s.match(/\.lite\.sra$/)}
log.info "Found #{sralites.length} SRA lite files in #{options[:in_directory]}"

# Setup progress bar
progress = ProgressBar.new('acacia', sralites.length)

# For each input file,
sralites.each do |sralite|
  # convert to a fastq file in a temporary directory
  Bio::FastqDumper.dump_in_tmpdir(File.join options[:in_directory], sralite) do |fastq|
    # mkdir and cd to an aptly named direct in the output folder
    outdir = File.join(options[:out_directory], sralite)
    if File.exist?(outdir)
      log.info "Skipping SRA file #{sralite}, since it already appears to be processed (output directory exists)"
      next
    end
    log.debug "Creating directory #{options[:out_directory]}"
    Dir.mkdir outdir
    # Dir.chdir options[:in_directory] #acacia cannot handle the fastq file not being in the current directory, I think
    
    # Run acacia and output the files into the working directory
    f = 
    acacia_command = [
      'java',
      '-jar',
      #'/srv/whitlam/home/projects/DATA_MODELLING/acacia-1.50.b04.jar'
      '/home/ben/bioinfo/acacia/acacia-1.50.b05.jar',
      %w(-DFASTA=FALSE -DFASTQ=TRUE),
      "-DOUTPUT_DIR=#{outdir}",
      %w(-DOUTPUT_PREFIX=acacia -DMID_FILE=null),
      "-DFASTQ_LOCATION=#{File.absolute_path fastq}"
    ].flatten
    log.debug "Running acacia with '#{acacia_command}'"
    Bio::Command.call_command_open3(acacia_command) do |stdin, stdout, stderr|
      err = stderr.read
      raise err if err != ''
      p stdout
    end
    puts `ls`
    
    progress.inc
  end
end
progress.finish





