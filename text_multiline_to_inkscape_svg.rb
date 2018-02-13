#!/usr/bin/env ruby

require 'optparse'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME}

Takes lines of text on stdin, and prints out them in a form suitable for pasting into an Inkscape .svg file, where each line is in a separate text object\n\n"
end; o.parse!

#TODO what are you looking at me for? This is your script. Do something.
y = 0
puts '<svg>'
ARGF.each do |line|
  puts "<text y=\"#{y}\">#{line.chomp}</text>"
  y += 10
end
puts '</svg>'
