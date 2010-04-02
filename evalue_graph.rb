#!/usr/bin/env ruby

# Make a graph from blast results showing the distribution of e-values obtained in the best hit

require 'tempfile'

Tempfile.open('evalue_data_points') do |tempfile|
  # Dump all the evalues without the log scale to a tempfile,
  # which will then be graphed
  last_entry = 'definitely not a gene!'
  $stdin.each do |line|
    splits = line.split "\t"
    raise unless splits.length == 12
    
    # Only record one e-value per query
    next if last_entry == splits[0]
    last_entry = splits[0]
    
    evalue = splits[10].to_f
    evalue = 1.0e-200 if evalue == 0.0 #this appears to be about the limit of BLAST, at the moment
    
    
    tempfile.puts Math.log10(evalue).round
  end    
  tempfile.close
  
  # Write a file that contains the commands to use in R
  pdf_name = '/tmp/evalues.pdf'
  r_commands = [
"data=read.table('#{tempfile.path}')",
"pdf('#{pdf_name}')",
"hist(data[,1], xlab='-log(E-value)', title='histogram of E-values')",
"dev.off()",
  ]
  Tempfile.open('evalue_r_cmds') do |tempfile|
    r_commands.each do |cmd|
      tempfile.puts cmd
    end
    tempfile.close
    
    `R --no-save <#{tempfile.path}`
  end
  
  # Open the pdf for the viewing by the user
  `gnome-open '#{pdf_name}'`
end