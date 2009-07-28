#!/usr/bin/env ruby

$stdin.each do |line|
if (matches = line.string.match(/(>.*? ).*/)
puts matches[1]
else
puts line.strip
end
