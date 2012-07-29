#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'progressbar'
require 'parallel'
require 'peach'

$LOAD_PATH.unshift(File.join([ENV['HOME'],%w(git bioruby-sra_fastq_dumper lib)].flatten))
require 'bio-sra_fastq_dumper'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stdout',
  :acacia_jar => '/srv/whitlam/bio/apps/sw/acacia/1.51-beta/acacia-1.51.b02.jar', #'/srv/whitlam/bio/apps/sw/acacia/1.51-beta/acacia-1.51.b0.jar', #/srv/whitlam/home/projects/DATA_MODELLING/acacia-1.50.b05.jar',
  :num_threads => 1,
  :max_file_size => 1024*1024*10, #default to 10MB maximum
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>
    
    Take a directory of .lite.sra files, and correct them using acacia, outputing the results in a second directory filled with directories (one for each sra lite file) \n\n"

  opts.on("-i", "--input-directory DIRECTORY", "a directory full of SRA .lite.sra files [required (or -f)]") do |arg|
    options[:in_directory] = File.absolute_path(arg)
  end
  opts.on("-f", "--sralite-list-file FILE", "a directory to put output files in [required (or -i)]") do |arg|
    options[:sralites_list] = File.absolute_path(arg)
  end
  opts.on("-o", "--output-directory DIRECTORY", "a directory to put output files in [required]") do |arg|
    options[:out_directory] = File.absolute_path(arg)
  end
  
  opts.on("-p", "--parallel-threads NUM_THREADS", "How many threads to split things into (1 thread per acacia run) [default #{options[:num_threads]}]") do |arg|
    options[:num_threads] = arg.to_i
  end
  opts.on("-m", "--max-file-size NUM_BYTES", "Ignore files larger than this size [default #{options[:max_file_size]}]") do |arg|
    options[:max_file_size] = arg.to_i
  end
  opts.on("-a", "--acacia-jar ACACIA_JAR_PATH", "Use a non-standard acacia version [default #{options[:acacia_jar]}]") do |arg|
    options[:acacia_jar] = arg
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
if ARGV.length != 0 or options[:out_directory].nil? or (options[:in_directory].nil? and options[:sralites_list].nil?)
  $stderr.puts o
  exit 1
end
# Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


# Collect a list of the input sra lite files to process, advise
if options[:sralites_list]
  log.debug "Reading SRA lite files from #{options[:sralites_list]}"
  sralites = File.open(options[:sralites_list]).readlines.collect{|l| l.strip}
  log.info "Found #{sralites.length} SRA lite files in #{options[:sralites_list]}"
elsif options[:in_directory]
  log.debug "Inspecting #{options[:in_directory]}"
  sralites = Dir.entries(options[:in_directory]).select{|s| s.match(/\.lite\.sra$/)}
  log.info "Found #{sralites.length} SRA lite files in #{options[:in_directory]}"
end

# Setup progress bar
progress = ProgressBar.new('sralite-acacia', sralites.length)

# For each input file,
#Parallel.each(sralites, :in_threads => options[:num_threads]) do |sralite|
#sralites.peach(options[:num_threads]) do |sralite|
sralites.each do |sralite|  # convert to a fastq file in a temporary directory
  lite_path = sralite
  if options[:in_directory]
    lite_path = File.join options[:in_directory], sralite
  end
  
  # Ensure that the file actually exists
  unless File.exist?(sralite)
    log.error "Unable to find file #{sralite}, skipping"
    next
  end
  
  # Ensure that the file is not larger than the maximum allowable size
  size = File.size(sralite)
  if options[:max_file_size] < size
    log.info "Skipping SRA file #{sralite} because it is too big (#{size} vs. max #{options[:max_file_size]})"
    next
  end
  
    # mkdir and cd to an aptly named direct in the output folder
  outdir = File.join(options[:out_directory], File.basename(sralite))
  if File.exist?(outdir)
    log.info "Skipping SRA file #{sralite}, since it already appears to be processed (output directory exists)"
    next
  end
  log.debug "Creating directory #{outdir}"
  Dir.mkdir outdir
  
  # convert to a fastq file in a temporary directory
  begin
    Bio::FastqDumper.dump_in_tmpdir(lite_path) do |fastq|
      # Dir.chdir options[:in_directory] #acacia cannot handle the fastq file not being in the current directory, I think
      log.debug "In temporary directory '#{Dir.getwd}'"
      
      # Run acacia and output the files into the working directory
      acacia_command = [
        'java',
        '-XX:+UseConcMarkSweepGC', #use concurrent garbage collection, as suggested by Lauren.
        '-jar',
        options[:acacia_jar],
        %w(-DFASTA=FALSE -DFASTQ=TRUE),
        "-DOUTPUT_DIR=#{outdir}",
        %w(-DOUTPUT_PREFIX=acacia -DMID_FILE=null),
        "-DFASTQ_LOCATION=#{fastq}"
      ].flatten
      log.debug "Running acacia with '#{acacia_command}'"
      Bio::Command.call_command_open3(acacia_command) do |stdin, stdout, stderr|
        #stdin.close
        out = stdout.readlines
        err = stderr.read
        raise err if err != ''
      end
      
      progress.inc
    end
  rescue Exception => e #if fastq-dump fails
    log.error "Something went wrong while executing acacia on #{sralite}, skipping: #{e}"
    next
  end
end
progress.finish





