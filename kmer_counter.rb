#!/usr/bin/env ruby

require 'bio'
require 'pp'
require 'progressbar'
#gem 'bio-faster', :version => '=0.2.2' #later version drop support for fastA parsing, concentrating on FastQ
#require 'bio-faster'
#require 'parallel'

# Initialise the hash of the different
module Bio
  class Sequence
    class NA
      # Return the current object or its reverse complement, whichever
      # has the sequence that comes first in lexigraphical (alphabetical)
      # order
      def lowest_lexigraphical_form
        rev = self.reverse_complement
        to_s < rev.to_s ? self : rev
      end
    end

    class Kmer
      def self.empty_full_kmer_hash(k=4)
        return @empty_full_hash.dup unless @empty_full_hash.nil?
        
        counts = {}
                
        ordered_possibilities = %w(A T C G)
        keys = ordered_possibilities
        (k-1).times do
          keys = keys.collect{|k| ordered_possibilities.collect{|n| "#{k}#{n}"}.flatten}.flatten
        end
        
        keys.each do |key|
          counts[key] = 0
        end
        counts
      end
      
      def self.empty_kmer_hash(k=4)
        return @empty_hash.dup unless @empty_hash.nil?
        
        counts = {}

        # construct an array of all possible kmers
        ordered_possibilities = %w(A T C G)
        keys = ordered_possibilities
        (k-1).times do
          keys = keys.collect{|k| ordered_possibilities.collect{|n| "#{k}#{n}"}.flatten}.flatten
        end

        # remove keys that are not in their lowest lexigraphical form.
        # Because all possible kmers are already present, the reverse complement of these must already be present in the array
        keys.select! do |key|
          Bio::Sequence::NA.new(key).lowest_lexigraphical_form.to_s.upcase == key
        end
        
        keys.each do |key|
          counts[key] = 0
        end
        
        @empty_hash = counts
        return counts
      end
      
      def self.merge_down_to_lowest_lexigraphical_form(hash)
        new_hash = empty_kmer_hash
        hash.each do |kmer, count|
          key = Bio::Sequence::NA.new(kmer).lowest_lexigraphical_form.to_s.upcase
          new_hash[key] += count
        end
        return new_hash
      end
    end
  end
end

#if __FILE__ == $0
  require 'optparse'

  # Parse cmd line options
  USAGE = "Usage: kmer_counter.rb [-w window_size] [-W window_offset] [-m minimum_window_size] [--window-length] [-k kmer_length] [--contig-name] <fasta_filename>"
  options = {
    :window_size => 5000,
    :minimum_window_size => 2000,
    :window_offset => 5000,
    :kmer => 4,
    :contig_name => false,
    :sequence_length => false
  }

  OptionParser.new do |opts|
    opts.banner = USAGE

    opts.on("-w", "--window-size SIZE", "Length of the window to be used [default #{options[:window_size]}]") do |v|
      window = v.to_i
      unless window > 0
        raise Exception, "Unexpected window size specified: #{v} - it must be greater than 0 residues long!"
      end
      options[:window_size] = window
      options[:window_offset] = window
    end

    opts.on("-W", "--window-offset SIZE", "Length of the offset between windows [default #{options[:window_offset]}]") do |v|
      offset = v.to_i
      unless offset > 0
        offset = options[:window_isze]
      end
      options[:window_offset] = offset
    end

    opts.on("-m", "--minimum-window-size SIZE", "Length of the minimum window to be used [default #{options[:minimum_window_size]}]") do |v|
      window = v.to_i
      unless window > 0
        raise Exception, "Unexpected minimum window size specified: #{v} - it must be greater than 0 residues long!"
      end
      options[:minimum_window_size] = window
    end

    opts.on("-k", "--kmer-length SIZE", "Length of the kmer to be used [default #{options[:kmer]}]") do |v|
      window = v.to_i
      unless window > 0
        raise Exception, "Unexpected minimum window size specified: #{v} - it must be greater than 0 residues long!"
      end
      options[:kmer] = window
    end

    opts.on("-n", "--contig-name", "Output the contig name, on top of the default contig chunk name [default: #{options[:contig_name]}]") do |v|
      options[:contig_name] = true
    end

    opts.on("-l", "--window-length", "print the length of the window in the output [default #{options[:sequence_length]}]") do |v|
      options[:sequence_length] = true
    end
    
    opts.on("-t", "--threads NUM_THREADS", "Use this many threads [default #{options[:threads]}]") do |v|
      options[:threads] = v.to_i
      if options[:threads] < 1
        raise Exception, "Unexpected number of threads specified (after converting to integer) - '#{options[:threads]}'"
      end
    end
  end.parse!

  print "ID\t"

  print Bio::Sequence::Kmer.empty_kmer_hash(options[:kmer]).keys.join("\t")
  print "\tWindowLength" if options[:sequence_length]
  print "\tcontig" if options[:contig_name]
  puts

  orig = Bio::Sequence::Kmer.empty_full_kmer_hash(options[:kmer])
  process_window = lambda do |window,kmer,sequence_name,contig_name|
    counts = orig.dup
    num_kmers_counted = 0
    
    window.window_search(options[:kmer],1) do |tetranucleotide|
      str = tetranucleotide.to_s
      next unless str.gsub(/[ATGC]+/,'') == ''
      num_kmers_counted += 1
      counts[str]+=1
      #counts[Bio::Sequence::NA.new(tetranucleotide).lowest_lexigraphical_form.to_s.upcase] += 1
    end
    
    # Merge everything into lowest lexigraphical form
    new_counts = Bio::Sequence::Kmer.merge_down_to_lowest_lexigraphical_form counts
    
    print "#{sequence_name}"
    counts.keys.each do |tetramer|
      print "\t#{counts[tetramer].to_f/num_kmers_counted}"
    end
    print "\t#{window.length}" if options[:sequence_length]
    print "\t#{contig_name}" if options[:contig_name]
    puts
  end

  fasta_filename = ARGV[0]
  progress = ProgressBar.new('kmer_counter', `grep -c '>' '#{fasta_filename}'`.to_i)
  ff = Bio::FlatFile.open(fasta_filename) 
  #Parallel.each(ff, :in_processes => options[:threads]) do |sequence|
  ff.each do |sequence|
    window_counter = 0
    sequence.seq.window_search(options[:window_size],options[:window_offset]) do |window|
      process_window.call(window, options[:kmer], "#{sequence.definition}_#{window_counter}",sequence.definition)
      window_counter += 1
    end
    leftover_length = sequence.seq.length % options[:window_size]
    if leftover_length >= options[:minimum_window_size]
      process_window.call(
      sequence.seq[sequence.seq.length-leftover_length..sequence.seq.length],
      options[:kmer], "#{sequence.definition}_leftover_#{window_counter}",sequence.definition)
    end
    progress.inc
  end
  progress.finish
#end
