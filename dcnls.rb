#!/usr/bin/env ruby

require 'rubygems'
require 'bio'

module Bio
  class DCNLS
    # default window size to attribute an NLS to
    DEFAULT_NLS_LENGTH = 5
    
    DEFAULT_NUMBER_BASIC_RESIDUES_IN_NLS = 4
    
    # Given a Bio::Alignment object, 
    # return a list of (bioinformatics style 1-based) indexes
    # of the alignment that correspond to predicted NLSs.
    # 
    # bio_msa_object: the Bio::Alignment object filled with one
    # or more (aligned) sequences
    #
    # options: tweaks to the algorithm, given as a Hash. Acceptable Hash keys:
    # * :required_number_of_basic_residues - default 4
    #
    def predictions(bio_msa_object, options={})
      options[:required_number_of_basic_residues] ||= DEFAULT_NUMBER_BASIC_RESIDUES_IN_NLS
      options[:nls_length] ||= DEFAULT_NLS_LENGTH

      # Iterate through the columns, taking 5 columns at a time
      acceptable_indices = []
       (0..bio_msa_object[bio_msa_object.keys[0]].length-1).each do |start|
        # Test if the current 5 columns are acceptable as NLSs
        num_basics = []
        #puts
        #puts start
        bio_msa_object.alignment_collect do |seq|
          num_basic = 0
          subseq = seq[start..start+options[:nls_length]-1]
          subseq.each_char do |c|
            num_basic += 1 if %w(R K).include?(c)
          end
          num_basics.push num_basic
        end
        if num_basics.select{
            |n| n >= options[:required_number_of_basic_residues]
          }.length > num_basics.length/2
          acceptable_indices.push start+1
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

    opts.on(nil, "--nls-length LENGTH", "Require NLS predictions to be LENGTH residues long") do |length|
      length = length.to_i
      if length < 1
        raise Exception, "Unexpected --nls-length parameter given"
      end
      options[:nls_length] = length
    end

    opts.on("-b", '--basics BASICS', 'Required number of basic residues within a sub-alignment (or subsequence) of the alignment (or seqeunce) to qualify as an NLS. Deafult: 4') do |basics|
      i = basics.to_i
      if i > 5 or i < 1
        $stderr.puts "Unexpected number of basic residues specified: `#{i}' - has to be between 1 and 5"
        exit
      end
      options[:required_number_of_basic_residues] = i
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
        aln.alignment_slice(i-1..i+3).alignment_collect do |seq|
          subseq = seq.seq
          puts subseq
        end
        puts
      end
    end
  end
  
  # Setup options to be used regardless of if it a single sequence or MSA is to be used.
  tweaks = {}
  unless options[:required_number_of_basic_residues].nil?
    tweaks[:required_number_of_basic_residues] = options[:required_number_of_basic_residues]
  end
  unless options[:nls_length].nil?
    tweaks[:nls_length] = options[:nls_length]
  end
    
  if options[:single_sequence]
    # Predict each sequence individually, not using a MSA as input
    Bio::FlatFile.open(msa_path).each do |s|
      aln = Bio::Alignment.new([s.seq])
      print_results.call(Bio::DCNLS.new.predictions(aln, tweaks), aln)
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
    print_results.call(Bio::DCNLS.new.predictions(aln, tweaks), aln)
  end
end
