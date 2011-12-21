#!/usr/bin/env ruby

require 'rubygems'
require 'bio'

module Bio
  class DCNLS
    def predictions(bio_msa_object)
      # Iterate through the columns, taking 5 columns at a time
      acceptable_indices = []
       (0..bio_msa_object[bio_msa_object.keys[0]].length-1).each do |start|
        # Test if the current 5 columns are acceptable as NLSs
        num_basics = []
        #puts
        #puts start
        bio_msa_object.alignment_collect do |seq|
          num_basic = 0
          subseq = seq[start..start+4]
          subseq.each_char do |c|
            num_basic += 1 if %w(R K).include?(c)
          end
          # p num_basic
          num_basics.push num_basic
        end
        if num_basics.select{|n| n>3}.length > num_basics.length/2
          acceptable_indices.push start
        end
      end
      return acceptable_indices
    end
  end
end



if __FILE__ == $0
  require 'optparse'
  
  USAGE = "Usage: dcnls.rb [-q] [-s] <multiple_sequence_alignment_file>"
  options = {
    :verbose => true,
    :single_sequence => false,
  }
  o = OptionParser.new do |opts|
    opts.banner = USAGE
    
    opts.on("-q", "--quiet", "Don't print anything except the results") do
      options[:verbose] = false
    end
    
    opts.on("-s", "--single", "Compute results on each sequence in the provided fasta file individually (default: treat the file as a multiple sequence alignment)") do
      options[:single_sequence] = true
    end
  end
  o.parse!
  msa_path = ARGV[0]
  
  print_results = lambda do |acceptable_indices, aln|
    # Print results
    if acceptable_indices.empty?
      puts "no NLSs found"
    else
      puts "Found NLSs at positions: #{acceptable_indices.collect{|s| s+1}.join(',')}"
    end
    
    if options[:verbose]
      acceptable_indices.each do |i|
        aln.alignment_slice(i..i+4).alignment_collect do |seq|
          subseq = seq.seq
          puts subseq
        end
        puts
      end
    end
  end
  
  if options[:single_sequence]
    # Predict each sequence individually, not using a MSA as input
    Bio::FlatFile.open(msa_path).each do |s|
      aln = Bio::Alignment.new([s.seq])
      print_results.call(Bio::DCNLS.new.predictions(aln), aln)
    end
  else
    # Predict a multiple sequence alignment
    # open the MSA file (this part of the code taken from Christian Zmasek's tutorial at http://code.google.com/p/forester/wiki/PhyloBioRuby)
    seq_ary = Array.new
    ff = Bio::FlatFile.auto(msa_path)
    ff.each_entry do |entry|
      seq_ary.push(entry)
    end
    
    # Creates a multiple sequence alignment (possibly unaligned) named
    # 'seqs' from array 'seq_ary'.
    aln = Bio::Alignment.new(seq_ary)
    print_results.call(Bio::DCNLS.new.predictions(aln), aln)
  end  
end
