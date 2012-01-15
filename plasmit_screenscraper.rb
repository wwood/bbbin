#!/usr/bin/env ruby

require 'open3'
require 'progressbar'
require 'reubypathdb'

if __FILE__ == $0
  USAGE = 'plasmit_screenscraper.rb <fasta_file>'
  
  inputs = Bio::FlatFile.open(ARGF).entries
  progress = ProgressBar.new('plasmit',inputs.length)
  
  inputs.each do |aa|
    aa.definition
    progress.inc
    
    # only the first 24 amino acids are used, but given that the length
    # output recorded for a 24 amino acid length protein is 23, I'm playing
    # it safe here
    if aa.seq.length < 25
      $stderr.puts "Found sequence #{aa.definition} with fewer than 25 amino acids, ignoring"
    else
      command = "curl -s 'http://gecco.org.chemie.uni-frankfurt.de/cgi-bin/plasmit/runanalysis.cgi?output=simple&sequence=>dummy%0A#{aa.seq[0..24]}'"
      
      Open3.popen3(command) do |stdin, stdout, stderr, wait_thr|
        stdin.close
        
        result = stdout.readlines
        error = stderr.readlines
        
        unless error.length == 0
          raise Exception, "There appears to be a problem while screenscrapting: #{error}"
        end
        
        result.each do |res_line|
          next unless res_line.match(/Lines read with presumably/)
          matches = res_line.match(/>>(.*?)<\/TD><TD>(.*?)</)
          puts [
            aa.definition,
            matches[2]
          ].join("\t")
        end
      end
    end
    
    sleep 1 #don't hammer the server too much
  end
end