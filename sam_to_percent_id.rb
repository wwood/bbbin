#!/usr/bin/env ruby

# Take a SAM file, and output a CSV file with:
# ID of the query,position,percent_id,match_length

require 'bio-samtools'

# a SAM line hitting
# F6F53:106:426	16	gi|374477287|gb|JH600414.1|	6371	0	10M	*	0	0	TTTGCCGTAC	*	XT:A:R	NM:i:0	X0:i:5	X1:i:127	XM:i:0	XO:i:0	XG:i:0	MD:Z:10
# a SAM line failing
# F6F53:9:85	4	*	0	0	*	*	0	0	ATCTGGGGACACGGTTCGA	*


ARGF.each_line do |line|
  next if line.match(/^@/) #skip header lines
  splits = line.strip.split("\t")
  if 

