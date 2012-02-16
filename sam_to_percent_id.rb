#!/usr/bin/env ruby

# Take a SAM file, and output a CSV file with:
# ID of the query,position,percent_id,match_length

require 'bio-samtools'

sam = Bio::DB::SAM.open(:tam => ARGF)

