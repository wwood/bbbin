#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} fasta_file

      Take a fasta file which has Ns in it (like after a scaffolded assembly), and split
      wherever there is Ns. Print out a fasta file of the splits, encoding the positions of the Ns into the name

      e.g. if the input is
      >scaffold1
      AANNTGT

      Then the output would be
      >scaffold1_1_2
      AA
      >scaffold1_5_8
      TGT
        \n\n"

    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ARGV.length > 1
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


  Bio::FlatFile.foreach(ARGF) do |scaffold|
    raise if scaffold.seq[0] == 'N' #not implemented

    # Remove Ns at the end
    scaffold.seq.gsub!(/N+$/,'')

    not_N_starts = [0]
    not_N_stops = []
    state = 'not_N'

    # find positions of Ns
    # working in zero-based indices here
    position = 0
    scaffold.seq.to_s.each_char do |char|
      if state=='not_N' and char == 'N' #switch from non-N mode to N mode
        state = 'yes_N'
        not_N_stops.push position-1
      elsif state == 'yes_N' and char != 'N' #switch from N mode to non-N mode
        state = 'not_N'
        not_N_starts.push position
      end
      position += 1
    end
    not_N_stops.push position-1 #finish it off

    if not_N_starts.length != not_N_stops.length or state != 'not_N'
      raise "error check fail. I don't expect to end on an N, for instance"
    end

    # print out the fragments
    scaffold_basename = scaffold.entry_id
    last_stop = 1
    not_N_starts.each_with_index do |start, index|
      stop = not_N_stops[index]
      frag = scaffold.seq[start..stop]
      puts ">#{scaffold_basename}_#{index+1}of#{not_N_starts.length}_#{start+1}_#{stop+1}"
      puts frag
    end
  end
end #end if running as a script
























