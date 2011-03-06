#!/usr/bin/env ruby

# This script takes a fasta file of nucleotide sequences and then returns a fasta file that includes a 6 frame translation, 
# except that the 6 frames are returned as nucleotides, not amino acid sequences as is common. The IDs in the returned fasta
# file are as per the EMBOSS program transeq, so adding _1, _2 etc. after each ID (just the first word, not the whole name of the
# sequence.

require 'rubygems'
require 'bio'

if __FILE__ == $0
  Bio::FlatFile.foreach(ARGV[0]) do |seq|
    # translate forwards
    seq_fwd = "#{seq.seq}"
     (0..2).each do |offset|
      puts ">#{seq.entry_id}_#{offset+1}"
      puts seq_fwd #[0..seq_fwd.length-seq_fwd.length%3-1] #previously, I was removing the hanging nucleotides, but given transeq puts an X on the end when this happens, I won't bother
      seq_fwd.gsub!(/^./,'')
    end
    
    # translate backwards
    # Seems the transeq order of 6 frames isn't quite what I was expected. I was thinking
    ## fwd
    ## fwd minus first nucleotide
    ## fwd minus first two nucleotides
    ## reverse_complement
    ## reverse_complement minus first (was the last) nucleotide
    ## reverse_complement minus first two (was the last two) nucleotides
    #But they seem to have mixed up the last two of these, which is somewhat strange, and annoying to my script. Oh well. Such is life.

    seq_back = Bio::Sequence::NA.new(seq.seq).reverse_complement.upcase
    back_seqs = [
    seq_back,
    seq_back[2..seq_back.length],
    seq_back[1..seq_back.length]
    ]
    (0..2).each do |i|
      puts ">#{seq.entry_id}_#{i+4}"
      s = back_seqs[i]
      puts s #[0..s.length-s.length%3-1] #previously, I was removing the hanging nucleotides, but given transeq puts an X on the end when this happens, I won't bother
    end
  end
end