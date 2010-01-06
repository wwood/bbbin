#!/usr/bin/env ruby

# Given a tab-delimited ruby file, only output 10 hits for each query, because
# anything more just takes up too much time.

current_query_name = nil
current_query_count = 0;
current_hit_name = nil;

$stdin.each_line do |line|
  splits = line.strip.split("\t")
  raise Exception, "Unexpected blast tab format!" unless splits.length == 12
  q = splits[0]
  h = splits[1]

  if q != current_query_name
    current_query_name = q;
    current_query_count = 1;
    current_hit_name = h;
  elsif h != current_hit_name
    current_query_count += 1
    current_hit_name = h
  end

  # too many - don't do anything with these
  next if current_query_count > 10

  # let the others just pass through
  print line
end
