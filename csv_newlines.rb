#!/usr/bin/ruby

# Reads in a file, and removes newlines that are between \" characters.
# Original purpose is to remove newlines from csv files generated using open office from microarray 'database' excel files so sequences can be extracted.


# usage
if !ARGV[0] or ARGV[0].downcase === 'help'
  $stderr.puts "Usage: csv_newlines.rb <CSV>\n"
  $stderr.puts "The CSV file should have \" as the text delimiter. Field delimiter cannot be \" because the script works by finding odd and even pairs of \" characters on a line. No \" characters should be in the underlying data either - you need to check this."
  exit
end




f = File.open(ARGV[0], 'r')

str = f.gets
buffer = ''

while str
  
  #$stderr.print '.'
  str = str.strip
  
  # If there is an odd number of \" characters, then continue to get lines until fixed
  matched = str.scan(/\"/)
  if matched
    num = matched.length
  else
    num = 0
  end
  #$stderr.puts "matches: #{num}"
  
  if buffer === '' #if there was something already there
    #$stderr.puts "new buffer"
    if (num % 2) == 0
      puts buffer+str
      buffer = ''
    else
      buffer += str
    end
  else # if the start of a new record
    #$stderr.puts "old buffer"
    if (num % 2) == 0
      buffer += str
    else
      puts buffer+str
      buffer = ''
    end
  end
  
  
  str = f.gets
end


puts
