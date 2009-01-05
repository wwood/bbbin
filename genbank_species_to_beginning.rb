#!/usr/bin/env ruby

require 'bio'

Bio::FlatFile.open($stdin).each do |seq|
  # if it contains brackets, assume that is the species name
  if matches = seq.entry.match(/^\>(.*)\[(.+)\](.*)$/)
    puts ">#{matches[2].gsub(' ', '_')}|#{matches[1]}#{matches[3]}"
  else
    $stderr.puts "Couldn't parse: #{seq.entry}"
    puts "#{seq.entry}"
  end
  puts seq.seq
end

