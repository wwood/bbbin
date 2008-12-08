
require 'bio'



class AminoAcid
  POLARS = %w(F M L V W P A C I G)
  POSITIVES = %w(K H R)
  NEGATIVES = %(D E)
  
  def initialize(letter)
    @letter = letter
  end
  
  def polar?
    POLARS.include?(@letter)
  end
  
  def positive?
    POSITIVES.include?(@letter)
  end
  
  def negative?
    NEGATIVES.include?(@letter)
  end
  
  def unknown?
    @letter == '*'
  end
end



amino_acids = %w(A R N D C Q E G H I L K M F P S T W Y V B J Z X *)


puts "    #{amino_acids.join('  ')}"

amino_acids.each do |aa1|
  print "#{aa1} "
  aac1 = AminoAcid.new(aa1)
  amino_acids.each do |aa2|
    aac2 = AminoAcid.new(aa2)
    answer = nil
    if aac1.unknown? or aac2.unknown?
      answer = -1
    elsif aac1.positive?
      if aac2.positive?
        answer = 1
      elsif aac2.negative?
        answer = -1
      else
        answer = -1
      end
    elsif aac1.negative? # no negatives are non-polar
      if aac2.positive?
        answer = -1
      elsif aac2.negative?
        answer = 1
      else
        answer = -1
      end
    elsif aac1.polar? #uncharged but polar
      if aac2.positive?
        answer = -1
      elsif aac2.negative?
        answer = -1
      elsif aac2.polar?
        answer = 1
      else
        answer = -1
      end
    else #uncharged and non-polar
      if aac2.positive?
        answer = -1
      elsif aac2.negative?
        answer = -1
      elsif aac2.polar?
        answer = -1
      else
        answer = 1
      end
    end
    
    if answer == -1
      print " -1"
    else
      print "  #{answer}"
    end
  end
  puts
end