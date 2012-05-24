#!/usr/bin/env ruby

require 'optparse'
require 'open3'
require 'bio'
require 'bio-logger'
require 'progressbar'

$:.unshift File.join(ENV['HOME'],'git','bioruby-hmmer_model','lib')
require 'bio-hmmer_model'

$:.unshift File.join(ENV['HOME'],'git','bioruby-hmmer3_report','lib')
require 'bio-hmmer3_report'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Description of what this program does...
"
  opts.on("-t", "--cds CDS_FILE", "A fasta file of the transcripts corresponding to each protein [required]") do |f|
    options[:cds_filename] = f
  end
  opts.on("-h", "--hmms HMM_DATABASE", "A path to where the database of HMMs used to search against the protein database is. This program will be much faster if you have this file index with 'hmmfetch --index' [required]") do |f|
    # Remote the .ssi at the end if it is given, since this isn't passed to the hmmfetch call later on
    options[:hmms_filename] = f.gsub(/.ssi$/,'')
  end
  opts.on("-r", "--hmmsearch-result RESULT", "A path to the result of the hmmsearch that was run previously [required]") do |f|
    options[:results_filename] = f
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

# Setup logging
Bio::Log::CLI.logger(options[:logger]) #bio-logger defaults to STDERR not STDOUT, I disagree
log = Bio::Log::LoggerPlus.new(LOG_NAME)
Bio::Log::CLI.configure(LOG_NAME)

# Cache the transcript sequences in a hash
transcripts = {}
Bio::FlatFile.foreach(options[:cds_filename]) do |seq|
  name = seq.definition.split(' ')[0]
  raise unless transcripts[name].nil?
  transcripts[name] = seq.seq
end
log.debug "Cached #{transcripts.length} sequences e.g. #{transcripts.keys[0]} => #{transcripts[transcripts.keys[0]]}" if log.debug?




# Initialise codon probabilities
probabilities = {}
ordered_possibilities = %w(A T C G)
keys = ordered_possibilities
2.times do
  keys = keys.collect{|k| ordered_possibilities.collect{|n| "#{k}#{n}"}.flatten}.flatten
end
keys.each do |key|
  probabilities[key] = {}
end
AA = %w(A C D E F G H I K L M N P Q R S T V W Y)
AA.each do |aa|
  probabilities.keys.each do |codon|
    probabilities[codon][aa] = 0.0
  end
end


progress = ProgressBar.new('hmm_codons', `grep '^//' '#{options[:results_filename]}' |wc -l`.to_i)

# For each query HMM in the alignment
report_count = 0
Bio::HMMER::HMMER3.reports(File.open(options[:results_filename])) do |report|
  hits = report.hits
  parsed_hmm = nil
  
  # For each domain in each alignment
  hits.each do |hit|
    hit.hsps.each do |hsp|
      # Parse the HMM model in this report
      if parsed_hmm.nil? #Only need to parse it once for each report since the HMM doesn't change within one report
        hmm_id = report.query_accession
        # Extract the name of the HMM model that is being used in this alignment
        # Extract the HMM model from hmm index and feed it into the biogem parser
        
        command = "hmmfetch '#{options[:hmms_filename]}' '#{hmm_id}'"
        log.debug "Running command `#{command}`" if log.debug?

        # execute command and capture both stdout, and stderr
        Open3.popen3(command) do |stdin, stdout, stderr|
          stdin.close
          
          # log.debug 'before stderrror'
          # stderror = stderr.readlines
          # log.debug 'after stderrror'
          out = stdout.read
          
          hmms = Bio::HMMER::Model.models(out)
          unless hmms.length == 1
            raise "Unexpected number of HMMs found in the HMM database. Output of the command\n'#{command}'\nwas\n#{stdout}"
          end
          
          parsed_hmm = hmms[0]
        end
      end
      
      # Extract the part of the transcript sequence that is matching the HMM
      query_name = hit.sequence_name.split(' ')[0]
      from_to = (hsp.alifrom-1)*3 .. (hsp.ali_to)*3
      cds_match = transcripts[query_name][from_to]
      log.debug "CDS: #{cds_match.scan(/(...)/).join(' ')}" if log.debug?
      
      # Check that translating the cds gives the sequence that is reported in the HSP for this protein (minus the gaps in the alignment)
      no_gaps = hsp.flatseq.gsub('-','').upcase
      translated = Bio::Sequence::NA.new(cds_match).translate.to_s.upcase
      # In most cases the conceptually translated and protein sequences will match.
      # However, in rare cases in bacteria/archaea the start codon isn't ATG.
      unless translated == no_gaps or no_gaps.gsub(/^./,'M') == no_gaps
        log.error "Expected protein sequence #{no_gaps}, but found #{translated} instead. Skipping as I hope this is a rare case."
        next
      end
      
      # Simultaneously iterate through the alignment, the HMM model and the transcript
      position = 0
      hmmposition = hsp.hmmfrom
      seqposition = 0
      while position < hsp.hmmseq.length
        # For each position where there is a an amino acid matching a HMM position
        # Record the probability of this position in the HMM matched to the codon
        if hsp.hmmseq[position] == '.' or hsp.flatseq[position] == '-'
          log.debug "Skipping position #{position} in alignment because there was a gap" if log.debug?
        else
          codon = cds_match[seqposition..(seqposition+2)]
          
          unless hsp.flatseq[position] == Bio::Sequence::NA.new(codon).translate.to_s or (position == 0 and hsp.flatseq[position] == 'M')
            raise "Unexpected translation found for codon #{codon}"
          end
          log.debug "Found match for codon #{cds_match[seqposition..(seqposition+2)]} #{hsp.flatseq[position]}, which matched #{hsp.hmmseq[position]}" if log.debug?
          
          cur_probabilities = parsed_hmm.match_emissions[hmmposition-1]
          parsed_hmm.alphabet.each_with_index do |aa, i|
            probabilities[codon][aa] += cur_probabilities[i] #parsed_hmm.match_probability(hmmposition, aa)
          end
        end
        hmmposition += 1 unless hsp.hmmseq[position] == '.'
        seqposition += 3 unless hsp.flatseq[position] == '-'
        position += 1
      end
    end
  end
  progress.inc
  #report_count += 1
  #break if report_count > 100
end
progress.finish


# Output
print "\t"
puts [
  AA,
  'max',
  'expected',
].join("\t")
probabilities.each do |codon, aa_probs|
  print codon
  AA.each do |aa|
    print "\t#{aa_probs[aa]}"
  end
  print "\t"
  print aa_probs.max{|ac1,ac2| ac1[1]<=>ac2[1]}[0]
  
  print "\t"
  expected = Bio::Sequence::NA.new(codon).translate
  print expected
  puts
end
