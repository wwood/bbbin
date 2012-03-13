#!/usr/bin/env ruby

require 'bio'
require 'pp'

# Initialise the hash of the different
module Bio
  class Sequence
    class Kmer
      def self.empty_kmer_hash(k=4)
        counts = {}

        keys = %w(A T G C)
        (k-1).times do
          keys = keys.collect{|k| %w(A T G C).collect{|n| "#{k}#{n}"}.flatten}.flatten
        end

        new_keys = []
        keys.each do |key|
          unless new_keys.include?(Bio::Sequence::NA.new(key).reverse_complement.to_s.upcase)
          new_keys.push key
          end
        end

        new_keys.each do |key|
          counts[key] = 0
        end
        return counts
      end
    end
  end
end

if __FILE__ == $0
  require 'optparse'

  # Parse cmd line options
  USAGE = "Usage: kmer_counter.rb [-w window_size] [-W window_offset] [-m minimum_window_size] <fasta_filename>"
  options = {
    :window_size => 5000,
    :minimum_window_size => 2000,
    :window_offset => 5000,
    :kmer => 4,
    :contig_name => false,
  }

  OptionParser.new do |opts|
    opts.banner = USAGE

    opts.on("-w", "--window-size SIZE", "Length of the window to be used") do |v|
      window = v.to_i
      unless window > 0
        raise Exception, "Unexpected window size specified: #{v} - it must be greater than 0 residues long!"
      end
      options[:window_size] = window
    end
    
    opts.on("-W", "--window-offset SIZE", "Length of the offset between windows") do |v|
      offset = v.to_i
      unless offset > 0
        offset = options[:window_isze]
      end
      options[:window_offset] = offset
    end
           
    opts.on("-m", "--minimum-window-size SIZE", "Length of the minimum window to be used") do |v|
      window = v.to_i
      unless window > 0
        raise Exception, "Unexpected minimum window size specified: #{v} - it must be greater than 0 residues long!"
      end
      options[:minimum_window_size] = window
    end

    opts.on("-k", "--kmer-length SIZE", "Length of the kmer to be used") do |v|
      window = v.to_i
      unless window > 0
        raise Exception, "Unexpected minimum window size specified: #{v} - it must be greater than 0 residues long!"
      end
      options[:kmer] = window
    end
    
    opts.on("-n", "--contig-name", "Output the contig name, on top of the default contig chunk name [default: false]") do |v|
      options[:contig_name] = true
    end
  end.parse!

  print "ID\t"

  print Bio::Sequence::Kmer.empty_kmer_hash(options[:kmer]).keys.join("\t")
  print "\tWindowLength"
  print "\tcontig" if options[:contig_name]
  puts

  process_window = lambda do |window,kmer,sequence_name,contig_name|
    counts = Bio::Sequence::Kmer.empty_kmer_hash(kmer)
    num_kmers_counted = 0
    window.window_search(options[:kmer],1) do |tetranucleotide|
      next unless tetranucleotide.upcase.gsub(/[ATGC]+/,'') == ''
      num_kmers_counted += 1
      if counts.key?(tetranucleotide.upcase.to_s)
      counts[tetranucleotide.upcase.to_s] += 1
      else
      counts[Bio::Sequence::NA.new(tetranucleotide).reverse_complement.upcase.to_s] += 1
      end
    end
    print "#{sequence_name}"
    counts.keys.each do |tetramer|
      print "\t#{counts[tetramer].to_f/num_kmers_counted}"
    end
    print "\t#{window.length}"
    print "\t#{contig_name}" if options[:contig_name]
    puts
  end

  Bio::FlatFile.open(ARGF).each do |sequence|
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
  end

end
