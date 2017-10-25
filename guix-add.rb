#!/usr/bin/env ruby

require 'bio-commandeer'

require 'optparse'
require 'bio-logger'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
$:.unshift File.join(File.dirname(__FILE__),'..','lib')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :git_commit_arguments => '',
  :coauthor => false,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    helper for writing GNU Guix commit messages\n\n"

  opts.on("-m", "--many", "Use the first pkg when >1 are found [default: croak]") do
    options[:many] = true
  end
  opts.on("-c", "--commit-arguments ARGS", "Pass these arguments to 'git commit' [default: #{options[:git_commit_arguments]}]") do |r|
    options[:git_commit_arguments] = r
  end
  opts.on("--staged", "Only commit staged changes [default: use git commit -a]") do
    options[:staged_only] = true
  end
  opts.on("--coauthor", "Add 'Co-authored by ...' to the commit message [default: #{options[:coauthored]}]") do
    options[:coauthor] = true
  end
  opts.on("--package NAME", "Add this package [default: guess]") do |arg|
    options[:package] = arg
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

diff = Bio::Commandeer.run "git diff #{options[:staged_only] ? '--cached' : ''}"

reg = /^\+\(define-public ([^ ]+)/

news = diff.split("\n").select{|l| l.match reg}.collect{|l| l.match(reg)[1]}

puts "Found #{news.length} new packages: #{news.join(', ')}"
top_pkg = nil
pkg = nil
s = ''
if options[:many] and news.length > 1
  top_pkg = news[0]
  pkg = news.join(', ')
  s = 's'
elsif news.length == 1
  pkg = news[0]
  top_pkg = pkg
elsif options[:package]
  pkg = options[:package]
  top_pkg = pkg
else
  raise "Did not find exactly 1 new package, aborting"
end

lines = Bio::Commandeer.run("git diff #{options[:staged_only] ? '--cached' : ''} --stat").split("\n")
raise lines if lines.length != 2
file = lines[0].strip.split(' ')[0]
msgs = [
  "gnu: Add #{top_pkg}.",
  '',
  "* #{file} (#{pkg}): New variable#{s}.\n"
]

if options[:coauthor]
  msgs.push "\nCo-authored-by: #{Bio::Commandeer.run('git config user.name').strip} <#{Bio::Commandeer.run('git config user.email').strip}>\n"
end

if msgs.select{ |m| m.length > 74 }.length > 1
  log.warn "Long lines detected, fixing manually required."
end
staged_options = '-a'
staged_options = '' if options[:staged_only]
Bio::Commandeer.run "git commit #{staged_options} -F - #{options[:git_commit_arguments]}", :stdin => msgs.join("\n")
