#!/usr/bin/env ruby

# make all the sequence names uniq in a PHYLIP file

NAME_CHARS = 10

number = 0
first = true
taken_names = {}

ARGF.each do |line|
  if first
    if line.strip.length == 0
      first = false
      puts
      next
    end
    
    
    if number == 0
      puts line
      number += 1
    else
      name = ''
      begin #until there is really a unique name
        number += 1
        spaces = Math::log10(number).floor+2
        
        (0..(NAME_CHARS-spaces)).each do |i|
          name<< line[i..i]
        end
        name<< "#{number}"
        name<< line[NAME_CHARS..(line.length-1)]
      end while taken_names[name]
      
      puts name
    end
  else
    puts line
  end
end
