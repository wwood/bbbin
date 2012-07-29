#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'progressbar'
require 'bio'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stdout',
    :db => '/srv/whitlam/bio/db/merged_gg_silva/merged_gg_silva.fna',
    :taxonomy => '/srv/whitlam/bio/db/merged_gg_silva/merged_gg_silva_taxo.txt',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} -f fasta_file_list.txt
      
      Takes a list of fasta files, and runs QIIME's assign_taxonomy.py on them\n\n"
      
    opts.on("-f", "--fasta-list FILE_OF_FASTA_FILES", "One path to a fasta file line per line [required]") do |arg|
      options[:fasta_files_list_file] = arg
    end
    opts.on("-d", "--db SEQUENCE_DATABASE", "Path to the basename of a bwa formatted database [default: #{options[:db]}]") do |arg|
      options[:db] = arg
    end
    opts.on("-t", "--taxonomy TAXONOMY_FILE", "Path to a taxonomy file defining the taxonomy in the -d/--db [default: #{options[:taxonomy]}]") do |arg|
      options[:taxonomy] = arg
    end

    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length != 0
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  files_to_interrogate = File.open(options[:fasta_files_list_file]).readlines.collect{|s| s.strip}.reject{|s| s==''}
  log.info "Found #{files_to_interrogate.length} fasta files to assign taxonomy to. Here we go.."

  progress = ProgressBar.new('assign_taxonomy',files_to_interrogate.length)
  files_to_interrogate.each do |fasta|
    assign_taxonomy_command = [
      'assign_taxonomy.py',
      '-m',
      'bwasw',
      '-t',
      options[:taxonomy],
      '-i',
      fasta,
      '-d',
      options[:db],
      '--threads',
      '24',
      '-o',
      File.dirname(fasta),
    ].flatten
    log.debug "Running assign_taxonomy with '#{assign_taxonomy_command}'"
    Bio::Command.call_command_open3(assign_taxonomy_command) do |stdin, stdout, stderr|
      err = stderr.read
      raise err if err != ''
    end
    progress.inc
  end
  progress.finish
end #end if running as a script




