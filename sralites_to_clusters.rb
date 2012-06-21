#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'

$LOAD_PATH.unshift('/home/ben/git/bioruby-sra_fastq_dumper/lib')
require 'bio-sra_fastq_dumper'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :sra_run_folder => '.',#"#{raise}",
  :target_folder => '.',#"#{raise}"
}
o = OptionParser.new do |opts|
  #TODO Fill in usage, description and option parsing below
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>
    
    Description of what this program does...\n"
  # Example option
  opts.on("-e", "--eg", "description [default: #{options[:eg]}]") do |f|
    options[:operation] = OVERALL
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





# Iterate over each run in the run folder
Dir.foreach(options[:sra_run_folder]) do |entry|
  p entry
  next if %w(. ..).include?(entry)
  
  if entry.match(/.sra$/)
    log.debug "Processing #{entry}"
  else
    log.debug "Ignoring non sra-file #{entry}"
    next
  end

  # Create a new temporary folder using the ruby plugin
  Bio::FastqDumper.dump_in_tmpdir(entry) do |fastq|  
    # run acacia
    command = [
      'java',
      '-jar',
      #'/srv/whitlam/home/projects/DATA_MODELLING/acacia-1.50.b04.jar'
      '/home/ben/bioinfo/acacia/acacia-1.50.b05.jar',
      %w(-DFASTA=FALSE -DFASTQ=TRUE),
      "-DOUTPUT_DIR=#{Dir.getwd}",
      %w(-DOUTPUT_PREFIX=acacia -DMID_FILE=null),
      "-DFASTQ_LOCATION=#{fastq}"
    ].flatten
    log.debug "Running acacia with '#{command}'"
    puts `ls`
    Bio::Command.call_command_open3(command) do |stdin, stdout, stderr|
      err = stderr.read
      raise err if err != ''
      p stdout
    end
    puts `ls`
    
    # run uclust, saving the results to a new file in the target folder.
    # sort the sequences using uclust
    command = [
      'usearch',
      '-sort',
      "acacia_#{fastq}.fa",
      '-output',
      'seqs.sorted.fasta',
    ]
    Bio::Command.call_command_open3(command) do |stdin, stdout, stderr|
      err = stderr.read
      raise err if err != ''
    end
    puts `ls`
        
    # cluster the sequences with uclust
    command = [
      'usearch',
      '-cluster',
      'seqs.sorted.fa',
      '-output',
      '/tmp/ta.uc'
    ]
    Bio::Command.call_command_open3(command) do |stdin, stdout, stderr|
      err = stderr.read
      raise err if err != ''
    end
    puts `ls`
  end
end





