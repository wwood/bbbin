#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :no_headers => false,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <arguments>
      
      Take one or more fasta files corresponding to transcriptomes, and output a table of codon usage counts. \n\n"
      
    opts.on("--no-headers", "Don't print the header row [default #{options[:no_headers]}]]") do |f|
      options[:no_headers] = true
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
  
  nucleotides = %w(A T G C)
  codons = nucleotides
  codons = codons.collect{|c| nucleotides.collect{|n| "#{c}#{n}"}}.flatten
  codons = codons.collect{|c| nucleotides.collect{|n| "#{c}#{n}"}}.flatten
  
  # Print headers
  puts [
    'transcriptome',
    codons.collect{|codon| "#{codon}/#{Bio::Sequence::NA.new(codon).translate}"}
  ].flatten.join("\t") unless options[:no_headers]
  
  ARGV.each do |fasta_file|
    if !File.exists?(fasta_file)
      log.error "Unable to find fasta file #{fasta_file}, skipping"
      # Don't exit because otherwise using it with xargs doesn't work
      next
    end
    
    codon_counts = {}
    num_seqs = 0
    Bio::FlatFile.foreach(fasta_file) do |s|
      s.seq.window_search(3) do |codon|
        codon_counts[codon] ||= 0
        codon_counts[codon] += 1
      end
      num_seqs += 1
    end
    log.info "Processed #{codon_counts.collect{|c,num|num}.reduce(:+)} codons from #{num_seqs} sequences in #{fasta_file}"
    
    # Print the final table
    puts [
      fasta_file,
      codons.collect{|codon|
        c = codon_counts[codon]
        c ||= 0
        c
      }
    ].flatten.join("\t")
  end
  
  
          
end #end if running as a script