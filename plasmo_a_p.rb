#!/usr/bin/ruby

require 'bio'
require 'signalp'

# PlasmoAP is a program for predicting apicoplast targetting sequences
# in Plasmodium falciparum (malaria).
# 
module Bio
  class PlasmoAP
    
    # Calculate the PlasmoAP score for a sequence (a string of amino acids)
    def calculate_score(sequence)
      # to_s means the sequence can be amino acid string or proper Bio::Sequence::AA object
      signalp = SignalSequence::SignalPWrapper.new.calculate(sequence.to_s)
      
      # Whether it passes one or both of the set tests
      set1 = set2 = false
      
      signal = signalp.classical_signal_sequence?
      return PlasmoAPResult.new(0) if  !signal#Both set of rules need a signal peptide
      cleaved = Bio::Sequence::AA.new(signalp.cleave(sequence))
      
      set1 = set1?(cleaved)
      set2 = set2?(cleaved)
      additional = additional?(cleaved)
      
      
      points = 0
      points += 2 if set1
      points += 2 if set2
      points += 1 if additional
      return PlasmoAPResult.new(points)
    end
    
    private
    
    #      Set 1: a sequence is considered ‘positive’ if i) it starts with a signal peptide, ii) the 15
    #amino acids following the predicted signal peptide do not contain more than 2 acidic
    #residues, iii) the 80 amino acids following the predicted signal peptide cleavage site
    #contain a stretch of 40 amino acids with a total of at least 9 asparagines and/or lysines,
    #and iv) the asparagine/lysine-enriched region has a ratio of basic to acidic residues of at
    #least 5 to 3.
    def set1?(cleaved_amino_acids)
      set1 = false
      return false if cleaved_amino_acids.length <= 15
      aa = Bio::Sequence::AA.new(cleaved_amino_acids[0..14])
      if aa.acidic_count <= 2 # ii)
        aa = Bio::Sequence::AA.new(cleaved_amino_acids[15..15+80-1])
        containing = nil
        # iii) contain a stretch of 40 amino acids with a total of at least 9 asparagines and/or lysines
        aa.window_search(40) do |window|
          if !containing #only the first one is needed
            comp = window.composition
            if comp['N']+comp['K'] >= 9
              containing = window
            end
          end
        end
        
        if containing
          if containing.basic_count.to_f / containing.acidic_count.to_f >= 5.0/3.0 # iv)
            set1 = true
          end
        end
      end
      
      return set1
    end
    
    
    #        Set 2: a sequence is considered ‘positive’ if it i) starts with a signal peptide, ii) if the 22
    #amino acids following the predicted signal peptide cleavage site exhibit a ratio of basic to
    #acidic residues of at least 10 to 7, iii) if the 80 amino acids following the predicted signal
    #peptide cleavage site contain a stretch of 40 amino acids with a total of at least 9
    #asparagines and/or lysines, and iv) if the asparagine/lysine-enriched region has a ratio of
    #basic to acidic residues of at least 10 to 9.
    def set2?(cleaved_amino_acids)
      set2 = false
      return false if cleaved_amino_acids.length <= 21
      aa = Bio::Sequence::AA.new(cleaved_amino_acids[0..21])
      if aa.basic_count.to_f / aa.acidic_count.to_f >= 10.0/7.0 #ii)
        
        # iii) if the 80 amino acids following the predicted signal
        #peptide cleavage site contain a stretch of 40 amino acids with a total of at least 9
        #asparagines and/or lysines
        aa = Bio::Sequence::AA.new(cleaved_amino_acids[22..22+80-1])
        containing = nil
        aa.window_search(40) do |window|
          if !containing #only the first one is needed
            comp = window.composition
            if comp['N']+comp['K'] >= 9
              containing = window
            end
          end
        end
        
        if containing
          if containing.basic_count.to_f / containing.acidic_count.to_f >= 10.0/9.0 # iv)
            set2 = true
          end
        end
      end
      return set2
    end
    
    
    # Additional point
    # basic => nil for none, true for basic, false for acidic
    def additional?(cleaved_amino_acids)
      cleaved_amino_acids.window_search(1) do |aa|
        if aa.basic_count == 1
          return true
        elsif aa.acidic_count == 1
          return false
        end
      end
      return nil
    end

  end # End class PlasmoAP
  
  
  
  class PlasmoAPResult
    attr_reader :points
    def initialize(points)
      @points = points
      raise Exception, "Bad PlasmoAP Score points: #{points}" if points < 0 or points > 5
    end
    
    def to_s
      case @points
      when 0..2
        return '-'
      when 3
        return '0'
      when 4
        return '+'
      when 5
        return '++'
      end
    end  
    
    def ==(another)
      @points == another.points
    end
    
    # '+' or '++' scores were taken as apicoplast targeted in the paper
    # does this result pass that test?
    def apicoplast_targeted?
      @points >= 4
    end
  end
end # End module Bio

if $0 == __FILE__
  runner = Bio::PlasmoAP.new
  
  # print out a list of proteins with yes/no answers
  puts [
    'Name',
    'PlasmoAP Score',
    'Apicoplast Targeted'
  ].join("\t")
  Bio::FlatFile.auto(ARGF).each do |seq|
    result = runner.calculate_score(seq.seq)
    print "#{seq.definition}\t#{result.to_s}\t"
    if result.apicoplast_targeted?
      puts 1
    else
      puts 0
    end
  end
end
