#!/usr/bin/env ruby

# make all the sequence names uniq in a PHYLIP file

DEFAULT_NAME_CHARS = 10

require 'optparse'

# ununiqify_ensembl_tree.rb ../ensembl.cbm48.fa uniqued.phylip consense.outtree
USAGE = "Usage: uniqify_phylip [-n] <phylip_file>"
options = {
  :name_chars => DEFAULT_NAME_CHARS
}
OptionParser.new do |opts|
  opts.banner = USAGE

  opts.on("-n", "--name-chars [NUMBER_CHARACTERS]", Integer, "Number of characters to print out per definition chunk") do |v|
    options[:name_chars] = v
  end
end.parse!

unless ARGV.length == 1
  $stderr.puts USAGE
  exit 1
end
name_chars = options[:name_chars]

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
        spaces = Math::log10(number).floor+2+(10-name_chars)
        
        (0..(name_chars-spaces)).each do |i|
          name<< line[i..i]
        end
        name<< "#{number}"
        name<< ' '
        name<< line[10..(line.length-1)]
      end while taken_names[name]
      
      taken_names[name] = true
      puts name
    end
  else
    puts line
  end
end
