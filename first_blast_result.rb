#!/usr/bin/env ruby

# Given a blast tab format out, pass through the pipes only the first hsp for each query/hit pair

last_query = nil
last_hit = nil
$stdin.each_line do |line|
  splits = line.strip.split("\t")
  if last_query == splits[0] and
      last_hit == splits[1]
    next
  end
  last_query = splits[0]
  last_hit = splits[1]
  puts line
end

