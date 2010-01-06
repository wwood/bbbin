#!/usr/bin/env ruby

# A script to take a blast xml output and split it into many smaller files,
# without using a rediculous amount of memory, which can happen if
# the whole big xml file is parsed at once.

state_machine_state = 'header' #header, middle

header = ''
ending = [ #all xml file end the same way
  '  </BlastOutput_iterations>',
  '</BlastOutput>'
]
current_filehandle = nil
current_filename_number = 1

start_xml_file = lambda{
  current_filehandle = File.open("#{current_filename_number}.xml", 'w')
  current_filename_number += 1
  current_filehandle.print header
  current_filehandle.print "    <Iteration>\n"
}
finish_xml_file = lambda {
  current_filehandle.puts ending.join("\n")
  current_filehandle.close
}

ARGF.each_line do |line| #for each line of the XML file
  # if we are beginning of a new kind block
  if line == "    <Iteration>\n"
    # if we are doing the first one, then record the header and don't print the currentxml, increment state
    if state_machine_state == 'header'
      state_machine_state = 'middle' #inc state
      start_xml_file.call
    else # if in the middle, print header, currentxml and ending
      finish_xml_file.call # finish the last one
      start_xml_file.call #start another
    end
  elsif line ==  '  </BlastOutput_iterations>'
    # at the end of the file. finish the current, then move state
    finish_xml_file.call # finish the last one
    break
  elsif state_machine_state == 'header'
    # continue recording header
    header += line
  else # just a regular middle line
    raise unless state_machine_state == 'middle' #pointless checking
    current_filehandle.print line
  end
end