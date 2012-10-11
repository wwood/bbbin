#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'tmpdir'
require 'tempfile'
require 'date'
require 'fileutils'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
}
o = OptionParser.new do |opts|
  #TODO Fill in usage, description and option parsing below
  opts.banner = "
    Usage: #{SCRIPT_NAME}
    
    Takes a folder full of sra lite files, and converts it into a blast database, of course.\n"
  
  opts.on("-f", "--folder FOLDER", "A folder full of sra-lite files [required]") do |f|
    options[:sra_folder] = f
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
if ARGV.length != 0 or options[:sra_folder].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]) #bio-logger defaults to STDERR not STDOUT, I disagree
log = Bio::Log::LoggerPlus.new(LOG_NAME)
Bio::Log::CLI.configure(LOG_NAME)


options[:sra_folder] =  File.expand_path options[:sra_folder]






# Create a fasta file with today's date as the name
date = DateTime.now
fasta_filename = "#{Date.today.strftime('%Y%m%d')}.fasta"
if File.exist?(fasta_filename) or File.exist?("#{fasta_filename}.nin")
  raise "A fasta file called #{fasta_filename} already exists, so I'm cowardly falling on my sword so I don't accidentally fall on my sword."
end
fasta = File.expand_path(fasta_filename)
log.debug "Creating a temporary FASTA file #{fasta} to hold all the sequences"

# For each sra file, create a temporary directory
original_directory = Dir.getwd
Dir.open(options[:sra_folder]).each do |file|
  next if %w(. ..).include?(file)
  
  matches = file.match(/^([A-Z\d]+)\.lite\.sra$/)
  unless matches
    log.info "Skipping file #{file} as it wasn't judged to be a lite sra file"
    next
  end
  
  Dir.mktmpdir do |dir|
    Dir.chdir dir
    
    log.debug "Inside temporary directory #{Dir.getwd}"
    # Then soft link the sra lite file into the temporary directory
    `ln -s #{File.join(options[:sra_folder], file)} .`
    log.debug "Soft-linking complete, now the temporary directory contains #{Dir.entries('.').join(', ')}"
    
    # un-sra-lite it with the ncbi tools
    `fastq-dump #{file}`
    log.debug "fastq-dump is done, now the temporary directory contains #{Dir.entries('.').join(', ')}"
    
    # convert to fasta, piping the output onto the big fasta file
    fastq_name = "#{matches[1]}"+'.fastq'
    unless File.exist?(fastq_name)
      log.error "Unexpected lack of fastq file #{fastq_name}, fail."
      next
    end
    command = "awk '{print \">\" substr(\$0,2);getline;print;getline;getline}' '#{fastq_name}' >>#{fasta}"
    log.debug "Finished fastq-dump, now converting to fasta format with #{command}" 
    `#{command}`
  
    # Clean up and remove the temporary directory
    log.debug "Finished with #{file}, moving on."
    Dir.chdir original_directory
  end
end


# Format the database with the fasta file
command = "makeblastdb -in #{fasta} -dbtype nucl"
log.info "Fasta file creation complete, now creating the blast database with '#{command}'"
`#{command}`

# delete the big fasta file
log.info "BLAST database built, now removing the fasta file"
FileUtils.rm(fasta_filename)


