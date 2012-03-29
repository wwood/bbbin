#!/usr/bin/env ruby

require 'bio'
require 'bio-samtools'

current_contig_name = nil
current_coverages = []

class Array
def median

return nil if empty?
len = length
sorted = sort
return (sorted[len/2] + sorted[(len+1)/2]) / 2
end

def mean
inject{ |sum, el| sum + el }.to_f / length
end
end


puts %w(SeqID median mean length).join("\t")
print_median = lambda do |contig, coverages|
  puts "#{contig}\t#{coverages.median}\t#{coverages.mean}\t#{coverages.length}"
end

ARGF.each do |line|
  pileup = Bio::DB::Pileup.new(line)
  if pileup.ref_name == current_contig_name
   # if the same just add the coverage
   current_coverages.push pileup.coverage
  else
   # if changed, print the last one (if there was one), setup the new one, add coverage
   print_median.call(current_contig_name, current_coverages) unless current_contig_name.nil?
   current_contig_name = pileup.ref_name
   current_coverages = [pileup.coverage]
  end
end

print_median.call(current_contig_name, current_coverages) unless current_contig_name.nil?
