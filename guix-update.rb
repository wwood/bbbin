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
  :staged => false,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Description of what this program does...\n\n"

  opts.on("-s", "--staged", "commit staged changes, not unstaged [default: #{options[:staged]}]") do
    options[:staged] = true
  end
  opts.on("-p", "--package STRING", "Set the package to be committed [default: autodetect]") do |arg|
    options[:package] = arg
  end
  opts.on("-v", "--version STRING", "Set the version to be committed [default: autodetect]") do |arg|
    options[:version] = arg
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


diff = nil
if options[:staged]
  diff = Bio::Commandeer.run "git diff --staged"
else
  diff = Bio::Commandeer.run "git diff"
end
log.debug "diff was:\n#{diff}"

name_reg = /^.*\(define-public ([^ ]+)/
version_reg = /^\+ +\(version \"([^ ]+)\"/

news = diff.split("\n").select{|l| l.match name_reg}.collect{|l| l.match(name_reg)[1]}
version_new = diff.split("\n").select{|l| l.match version_reg}.collect{|l| l.match(version_reg)[1]}

puts "Found #{news.length} packages: #{news.join(', ')}"
puts "Found #{version_new.length} versions: #{version_new.join(', ')}"
defined_pkg = options[:package]
new_version = options[:version]
if (version_new.length == 1 or new_version) and (news.length == 1 or defined_pkg)
  pkg = defined_pkg
  pkg ||= news[0]
  new_version ||= version_new[0]
  lines = Bio::Commandeer.run("git diff --stat").split("\n").select{ |l| l.match(/gnu\/packages\/[^\/]+\.scm/) }
  raise lines if lines.length != 1
  file = lines[0].strip.split(' ')[0]
  cmd = "git commit -F -"
  cmd += " -a" if not options[:staged]
  Bio::Commandeer.run cmd, :stdin => "gnu: #{pkg}: Update to #{new_version}.\n\n* #{file} (#{pkg}): Update to #{new_version}.\n"
else
  raise "Did not find exactly 1 package to update, aborting"
end
