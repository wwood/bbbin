#!/usr/bin/ruby




module Orf
  START_AMINO_ACID_CHAR = 'M'
  STOP_AMINO_ACID_CHAR = '*'
  
  class OrfFinder
    # return an array of filled out OrfThread objects. Each OrfThread is itself an array of orfs.
    # maximal orfs are returned. That is, 'MMYStop' will return MMYStop, not MYStop.
    # incomplete ORFs may also be included, if they fall out of the analysis
    def generate_longest_orfs(nucleotide_sequence)
      thread = OrfThread.new
    
      # convert to more readable format
      trans = [
        Bio::Sequence::NA.new(nucleotide_sequence).translate(1),
        Bio::Sequence::NA.new(nucleotide_sequence).translate(2),
        Bio::Sequence::NA.new(nucleotide_sequence).translate(3)
      ]
    
      phase_offset = 0
      thread_array = []
      trans.each do |seq|
        thread = OrfThread.new
      
        # retrieve the first bit that doesn't start with an M character
        if first = seq.match(/^([^M\*]*\*)/)
          o = Orf.new
          o.start = phase_offset + 3*first.offset(0)[0]
          o.stop = phase_offset + 3*first.offset(0)[1]-1
          o.aa_sequence = first[1]
          thread.push o
        end
      
        # retrieve the full orfs
        seq.scan(/(M.*?\*)/){ #do |match| doesn't work here - tests fail.
          o = Orf.new
          o.start = phase_offset + 3*$~.offset(0)[0]
          o.stop = phase_offset + 3*$~.offset(0)[1]-1
          o.aa_sequence = $~.to_s
          thread.push o
        }
      
        # retrieve the partial orf at the end if it exists
        if last = seq.match(/(M[^\*]*)$/)
          o = Orf.new
          o.start = phase_offset + 3*last.offset(0)[0]
          o.stop = phase_offset + 3*last.offset(0)[1]-1
          o.aa_sequence = last[1]
          thread.push o
        end
      
        phase_offset += 1
        thread_array.push thread
      end
    
      return thread_array
    end
    
    # Return the Orf object representing the longest orf. It is possible the orf is a fragment,
    # if it encounters the end (or start) of the sequence while in an orf.
    def longest_orf(nucleotide_sequence)
      thread_array = generate_longest_orfs(nucleotide_sequence)
      return nil if thread_array.empty?
      return thread_array.flatten.max
    end
    
    # Return the Orf object representing the longest 'full' orf. That is, the orf that contains both a start
    # and a stop codon
    def longest_full_orf(nucleotide_sequence)
      thread_array = generate_longest_orfs(nucleotide_sequence)
      return nil if thread_array.empty?
      best = thread_array.flatten.max{|a,b| 
        if a.fragment? and b.fragment? #if both fragments, meh. Keep some order though, otherwise other functions get confused maybe?
          a<=>b
        elsif a.fragment?#a is fragment but not b, therefore b wins
          -1
        elsif b.fragment?#b is fragment but not a, therefore a wins
          1
        else
          a<=> b #neither is fragments, so just the normal comparison is fine.
        end
      }
      if !best or best.fragment?
        return nil
      else
        return best
      end
    end
    
    # Longest ORF with a start codon, but not necessarilly a stop codon.
    def longest_m_orf(nucleotide_sequence)
      thread_array = generate_longest_orfs(nucleotide_sequence)
      return nil if thread_array.empty?
      best = thread_array.flatten.max{|a,b| 
        if a.start? and b.start? # both are start. good. please continue
          a<=>b
        elsif a.start?#a is fragment but not b, therefore b wins
          1
        elsif b.start?#b is fragment but not a, therefore a wins
          -1
        else
          a<=> b #if both fragments, meh. Keep some order though, otherwise other functions get confused maybe?
        end
      }
      if !best or !(best.start?)
        return nil
      else
        return best
      end
    end
  end

  class OrfThread <Array
    attr_accessor :orfs
  end

  class Orf
    attr_accessor :start, :stop, :aa_sequence
    
    def fragment?
      mystop = @aa_sequence.length-1
      !(
        (
          @aa_sequence[mystop..mystop] === STOP_AMINO_ACID_CHAR and
            @aa_sequence[0..0] === START_AMINO_ACID_CHAR
        )
      )
    end
    
    # This is a fragment from the end of an ORF?
    def end_fragment?
      mystop = @aa_sequence.length-1
      @aa_sequence[mystop..mystop] === STOP_AMINO_ACID_CHAR and
        @aa_sequence[0..0] != START_AMINO_ACID_CHAR
    end
    
    def start?
      @aa_sequence[0..0] === START_AMINO_ACID_CHAR
    end
    
    def <=>(orf_b)
      self.length <=> orf_b.length
    end

    def length
      @aa_sequence.length
    end
  end
end


# A simple script interface to iterate over a fasta file and return results
if $0 == __FILE__
  require 'bio'
  require 'optparse'
  
  
  def fasta_output(bioseq)
    finder = Orf::OrfFinder.new
    bioseq.each do |seq|
      orf = finder.longest_orf(seq.seq)
      #skip empty sequences
      next if !orf
      
      puts ">#{seq.entry_id}"
      if orf
        puts orf.aa_sequence
      else
        puts
      end
    end
  end
  
  
  options = ARGV.getopts("fh") #f for fasta, no arguments required. h is for help
  if ARGV.length > 1 or options[:h]
    $stderr.puts "Usage: orf_finder.rb [-f] <my.fasta>"
    $stderr.puts "Where my.fasta is the name of the fasta file you want to analyse. The output is the length of the longest ORF found in each sequence."
    $stderr.puts "-f: fasta. output a fasta file of the longest orfs found."
    return
  end
  
  # first argument or failing that, stdin to input the fasta file
  input = ARGV.length == 1 ? ARGV[0] : $stdin
  
  # if fasta is wanted as output, do that
  if options['f']
    fasta_output(Bio::FlatFile.auto(input))
  else

  
    if !options[:f]
      puts [
        'Name',
        'Longest ORF Length',
        'Longest ORF Sequence',
        'Longest Full ORF Length',
        'Longest Full ORF Sequence',
        'Longest Orf with Start Codon Length',
        'Longest Orf with Start Codon Sequence'
      ].join("\t")
    end
    finder = Orf::OrfFinder.new
    Bio::FlatFile.auto(input).each do |seq|
      orf = finder.longest_orf(seq.seq)
      to_print = []
    
      # any old orf
      if !orf
        to_print = [
          seq.entry_id,
          0,
          ''
        ]
      else
        to_print = [
          seq.entry_id,
          orf.length,
          orf.aa_sequence
        ]
      end
    
      #must be a full orf
      orf = finder.longest_full_orf(seq.seq)
      if !orf
        to_print.push nil
        to_print.push nil
      else
        to_print.push orf.length
        to_print.push orf.aa_sequence
      end
      
            #must be a full orf
      orf = finder.longest_m_orf(seq.seq)
      if !orf
        to_print.push nil
        to_print.push nil
      else
        to_print.push orf.length
        to_print.push orf.aa_sequence
      end
    
      puts to_print.join("\t")
    end
  end
  
  

  
  
else #included as module
  require 'bio'
end 