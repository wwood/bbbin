#!/usr/bin/env ruby

# A script for de-IUPAC coding DNA sequences

wobbles = {
  'R'=>['A','G'],
  'Y'=>['C','T'],
  'S'=>['G','C'],
  'W'=>['A','T'],
  'K'=>['G','T'],
  'M'=>['A','C'],
  'B'=>['C','G','T'],
  'D'=>['A','G','T'],
  'H'=>['A','C','T'],
  'V'=>['A','C','G'],
  'N'=>%w(A T G C)
}

ARGF.each do |line|
  splits = line.strip.split(/\s/)
  unless splits.length == 2
    raise Exception, "Expected line format: 'primer1_seq primer2_seq'"
  end

  dewobbler = lambda do |wobbled|
    primer1_roots = ['']
    wobbled.each_char do |char|
      char.upcase!
      if %w(A T G C).include?(char)
        primer1_roots = primer1_roots.collect do |root|
          "#{root}#{char}"
        end
      elsif wobbles[char]
        new_primers = []
        primer1_roots.each do |root|
          wobbles[char].each do |wob|
            new_primers.push "#{root}#{wob}"
          end
        end
        primer1_roots = new_primers
      else
        raise Exception, "Unexpected character detected: #{char}"
      end
    end
    primer1_roots
  end
  
  primer1s = dewobbler.call(splits[0])
  primer2s = dewobbler.call(splits[1])
  primer1s.each do |p1|
    primer2s.each do |p2|
      puts "#{p1} #{p2}"
    end
  end
end
