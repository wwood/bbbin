#!/usr/bin/env ruby
#
# Takes the output of an exonerate run (that has been run with the 
# --showtargetgff or the other showgff and gives prints out just the GFF part of
# the file
 

START_GFF = '^\# \-\-\- START OF GFF DUMP \-\-\-'
END_GFF = '^\# \-\-\- END OF GFF DUMP \-\-\-'

in_gff = false
ARGF.each do |line|
  case in_gff
  when false
    if line.match(/#{START_GFF}/)
      in_gff = true
    end
    
  when true
    if line.match(/#{END_GFF}/)
      puts
      in_gff = false
    else
      unless line.match(/^\#/)
        print line
      end
    end
  end
end
