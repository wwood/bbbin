#!/usr/bin/env ruby

# Convert a CSV file to a tab separated values file, keeping in
# mind that Gene annotation columns can contain commas. To get
# around this given as an argument is the number of columns expected
# and the number (0-indexed) of the annotation column.



class CsvToTab

def parse_line(line, num_cols, annotation_column)



# For each line on STDIN
require 'csv'
row = line.split(',')

if row.length != num_cols
last_annotation_column = row.length-(num_cols-annotation_column)-1
line = [
row[0..annotation_column-2],
row[annotation_column-1..last_annotation_column].join(","),
row[last_annotation_column+1..row.length-1]
].flatten.join("\t")
else
line = row.join("\t")
end

return line

end
end


# If running as a script, do that
if $0 == __FILE__

if ARGV.length != 2
$stderr.puts "Usage: #{$0} <num_columns> <annotation_column>"
exit
end

#parse args
num_cols = ARGV[0].to_i
annotation_column = ARGV[1].to_i

t = CsvToTab.new
$stdin.each do |line|
puts t.parse_line(line, num_cols, annotation_column)
end

end
