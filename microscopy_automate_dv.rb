#!/usr/bin/env ruby

# Given a folder wih multiple stacks present as individual slices, process those stacks so that they are all as one merged stack, and do a maximum projection of each.
# Requires:
# * ImageJ to be accessible on the cmd line.
# * Ruby and the rubygems fastercsv, peach
# * ImageMagick accessible on the cmd line (i.e. 'convert')


require 'rubygems'
require 'tempfile'
require 'optparse'
require 'fastercsv'
require 'peach'

USAGE = "Usage: microscopy_automate [-r <channel_colour_order e.g. 'rgb'>] [<directory>]";
options = {
  :channel_colour_order => 'rgb',
  :num_threads => 2,
}



base_dir = ARGV[0]
base_dir ||= '.' #cwd by default
output_dir = File.join(base_dir,'processed')
max_projection_dir = File.join(output_dir, 'maximum_projections')
scale_bar_max_projection_dir = File.join(output_dir, 'maximum_projections_scale_bar')
three_d_dir = File.join(output_dir, '3D_reconstructions')
Dir.mkdir output_dir unless File.exist?(output_dir)
Dir.mkdir max_projection_dir unless File.exist?(max_projection_dir)
Dir.mkdir scale_bar_max_projection_dir unless File.exist?(scale_bar_max_projection_dir)
Dir.mkdir three_d_dir unless File.exist?(three_d_dir)

stacks = {} #hash of stack id to filenames
stack_colours = %w(Red Green Blue)
stack_metadata = {} #hash of stack_id to hash of metadata key value pairs (as it comes in the lei files, for instance)

# First parse the lei file into a csv file, so that the 

# Process the directory, splitting the filenames up into stacks in the stacks hash
Dir.foreach(base_dir) do |filename|
  next unless filename.match(/.dv$/) #skip '.' and '..' etc.
  
  # convert from dv to tif stack
  filename_for_bfconvert = "#{filename}.stack_%w.tif"
  cmd = "bfconvert '#{filename}' '#{File.join(output_dir, filename_for_bfconvert)}'"
  puts "Running: #{cmd}"
  puts `#{cmd}`
  
  # get the paths to the individual stack files
  channel_paths = []
  Dir.foreach(output_dir) do |f|
    path = File.join(output_dir, f)
    channel_paths.push f if f.match(/#{filename}.stack_\d.tif$/)
  end
  channel_paths.sort!
  p channel_paths
  
  # Do maximal projection
  Tempfile.open(['ij_automate_macro', '.txt']) do |tempfile|# IJ does support commands directly on cmdline, but doesn't exit afterwards, annoyingly, so have to use a tempfile instead
    channel_paths.each do |channel_filename|
      channel_path = File.join(output_dir, channel_filename)
      tempfile.puts "open('#{channel_path}')\n"
      tempfile.puts "run(\"Z Project...\", \"start=1 stop=30 projection=[Max Intensity]\");\n" #HACKK!!
      fn = "#{channel_filename}.mp.png"
      tempfile.puts "saveAs('png','#{File.join(max_projection_dir, channel_filename)}')\n"
      tempfile.puts "close()\n"
    end
    
    tempfile.close
    puts `cat #{tempfile.path}`
    
    # Run ImageJ with the macro
    puts `imagej -b '#{tempfile.path}'`
  end
  
  # Clean up
  channel_paths.each do |ch|
    cmd = "rm '#{File.join output_dir, ch}'"
    puts "Running: #{cmd}"
    puts `#{cmd}`
  end
end
