#!/usr/bin/env ruby
require "pry"
require 'optparse'
require 'bio-logger'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

Takes a multi-hmm file, and creates individual ones. Right now the name of the file is based on the ACC line.\n\n"

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 and ARGV.length != 1
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME); log.outputters[0].formatter = Log4r::PatternFormatter.new(:pattern => "%5l %c %d: %m", :date_pattern => '%d/%m %T')


# HMMER3/b [3.0 | March 2010]
# NAME  KO:K10759 3.1.1.20
# TC    557.6 556.8;
# ACC   HMMsoil99735
# DESC  KO:K10759 3.1.1.20
# LENG  453
# ALPH  amino
# RF    no
# CS    no
# MAP   yes
# DATE  Tue Feb 21 14:05:44 2012
# NSEQ  35
# EFFN  0.798950
# CKSUM 2199757263
# STATS LOCAL MSV      -11.6917  0.69854
# STATS LOCAL VITERBI  -12.6172  0.69854
# STATS LOCAL FORWARD   -6.0051  0.69854
# HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
#             m->m     m->i     m->d     i->m     i->i     d->m     d->d
#   COMPO   2.51273  3.94250  2.89033  2.71593  3.29244  2.64821  3.69688  2.94100  2.79606  2.55776  3.66303  3.06618  3.36211  3.07492  3.06447  2.63679  2.79193  2.69531  4.34974  3.20266
#           2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
#           0.10298  4.07384  2.51522  0.61958  0.77255  0.00000        *
#       1   2.73679  0.71710  4.28335  4.05273  4.01669  3.33190  4.67225  3.18276  3.88082  3.09195  4.21786  4.00623  4.04867  4.22185  3.98238  3.00848  3.21428  2.91851  5.40470  4.30311      1 - -
#           2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
#           0.02763  3.99849  4.72084  0.61958  0.77255  0.56371  0.84186

log.info "Parsing.."

state = :before_acc
previous_lines = []
output_file = nil
num_found = 0

ARGF.each_line do |line|
  if state == :before_acc
    previous_lines << line
    if line[0..2]=='ACC'
      state = :after_acc
      if matches = line.match(/ACC\s+(.+)/)
        output_file = File.open "#{matches[1] }.hmm", 'w'
        output_file.print previous_lines.join
        num_found += 1
      else
        raise "Error parsing ACC line: #{line}"
      end
    end
  elsif state == :after_acc
    output_file.print line
    if line[0..1] == "//"
      output_file.close
      previous_lines = []
      output_file = nil
      state = :before_acc
    end
  end
end
log.info "Found #{num_found} HMMs"


