#!/usr/bin/ruby

# blastbig program queryfile blastdbname
# blastdbname is in /blastdb

# Usage: blastbig.sh <program/-p> <input sequences/-i> <blast database/-d> <output format/-m> <other options>

outputname = 'blastresult'
if ARGV[3].to_i == 7
    outputname += '.xml'
else
    outputname += '.txt'
end


cmd = "nice blastall -v10 -b10 -a 4 -p #{ARGV[0]} -i #{ARGV[1]} -d /blastdb/#{ARGV[2]} -m #{ARGV[3]} #{ARGV[4]} -o #{outputname}"

# use old style XML because then bioruby can handle it
# no need to anymore, now that the bioruby blast parser is fixed.
#if ARGV[3].to_i == 7
#   cmd += ' -V'
#end    

ENV['BLASTMAT'] = '/home/ben/bioinfo/blastdata'
puts cmd
system(cmd)
