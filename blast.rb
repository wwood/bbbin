#!/usr/bin/ruby

# a wrapper script so that blast2 arguments work with blastall, which
# requires creating a database instead of a -j output method
require 'optparse'

options = ARGV.getopts("j:p:e:i:")
if options.values.include?(nil)
  $stderr.puts "blast.rb: Failed to parse j, e, i or p argument, which are required."
  exit
end

# convert the fasta file to an argument
protein = 'F'
if ['blastx','blastp'].include?(options['p'])
  protein = 'T'
end
cmd = "formatdb -i #{options['j']} -p #{protein} -V </dev/null"
system cmd

# run the actual blast
cmd = "blastall  -i #{options['i']} -d #{options['j']} -p #{options['p']} -e #{options['e']}"
system(cmd)
