#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

Given a list of reactions\n\n"

  opts.on("--pairs ARG", "linked list pairs of main compontents [required]") do |arg|
    options[:pairs_file] = arg
  end
  opts.on("--metabolite-masses ARG", "masses file [required]") do |arg|
    options[:metabolites_file] = arg
  end
  opts.on("--output ARG", "output file [required]") do |arg|
    options[:output_file] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME); log.outputters[0].formatter = Log4r::PatternFormatter.new(:pattern => "%5l %c %d: %m", :date_pattern => '%d/%m %T')


log.info "Reading metabolite file"
masses = {}
CSV.foreach(options[:metabolites_file], :col_sep => "\t") do |row|
  raise unless row.length == 4
  compound_id = row[0]
  mass = row[3].to_f
  raise if masses[compound_id]
  masses[compound_id] = mass
end

ins = {}
outs = {}
num_zero_mass_reactions = 0
num_equal_mass_reactions = 0
log.info "Reading linked list"
CSV.foreach(options[:pairs_file], :col_sep => ' ') do |row|
  comp1 = row[0]
  comp2 = row[1]
  mass1 = masses[row[0]]
  mass2 = masses[row[1]]
  if mass1 == 0 or mass2 == 0
    num_zero_mass_reactions += 1
    next
  end
  if mass1 > mass2
    # degradation reaction (assume fwd)
    ins[comp2] ||= []
    ins[comp2].push comp1
    outs[comp1] ||= []
    outs[comp1].push comp2
  elsif mass1 < mass2
    ins[comp1] ||= []
    ins[comp1].push comp2
    outs[comp2] ||= []
    outs[comp2].push comp1
  else
    num_equal_mass_reactions += 1
  end
end
log.info "Found #{num_zero_mass_reactions} reactions that included a zero mass reactant e.g. those with R groups"
log.info "Found #{num_equal_mass_reactions} reactions that didn't change mass"

log.info "Finding reactions where the reactant compound is a reactant only in a single reaction and the product is only a product in a single reaction"
num = 0
File.open(options[:output_file], 'w') do |out|
  ins.each do |comp, inputs|
    if inputs.length == 1
      if outs[inputs[0]].length == 1
        out.puts [
          inputs[0],
          outs[inputs[0]][0],
          ].join("\t")
        num += 1
      end
    end
  end
end
log.info "Found #{num} reactions that fit"
