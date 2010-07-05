#!/usr/bin/env ruby

# This script is similar to blast. The difference is that it tries to remove 
# problems with one area of the protein hogging all the hits.
#


if __FILE__ == $0
  require 'optparse'
  
  # Parse options
  # should be able to set number deep at each point
  
  # run the initial blast
  # start the blast output file
  
  # while there are proteins still being hit
    # mask out the hit areas
    # redo the blast, this time using the masked sequence as input
    # print the current blast outputs to the big output file
  
  # Print out the big output file to STDOUT
end