#!/usr/bin/env ruby

require 'bio'
require 'bio-signalp'

# Looking for NeuroPeptides using basic approaches

puts [
'gene',
'cleaved molecular weight',
'cleaved matches KR?',
'cleaved matches RR?',
'cleaved matches KK?',
].join("\t")

signalp = Bio::SignalP::Wrapper.new
Bio::FlatFile.foreach('Aqu1.pep.fa') do |s|
  seek = s.seq.gsub(/\*/,'').gsub('X','')
  result = signalp.calculate(seek)
  if result.signal?
    cleaved_seq = result.cleave(seek)
    uncleaved_aa = Bio::Sequence::AA.new(seek)
    cleaved_aa = Bio::Sequence::AA.new(cleaved_seq)
    to_write = [
    s.definition
    ]
    
    aa = cleaved_aa
    to_write.push aa.molecular_weight
    to_write.push !(aa.to_s.upcase.match(/KR/).nil?)
    to_write.push !(aa.to_s.upcase.match(/RR/).nil?)
    to_write.push !(aa.to_s.upcase.match(/KK/).nil?)
    
    puts to_write.join("\t")
  end
end
