#!/usr/bin/env ruby

# make all the sequence names uniq in a PHYLIP file

NAME_CHARS = 10

number = 0
first = true
ARGF.each do |line|
  if first
    if line.strip.length == 0
      first = false
      puts
      next
    end
    
    number += 1
    if number == 1
      puts line
    else
      spaces = Math::log10(number).floor+1

      (0..(NAME_CHARS-spaces)).each do |i|
        print line[i..i]
      end
      print number-1
      print line[NAME_CHARS..(line.length-1)]
    end
  else
    puts line
  end
end
