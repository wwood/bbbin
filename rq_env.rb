#!/usr/bin/ruby
#require 'gem'
#require 'rio'

# A Wrapper around RQ so it stops making those annoying errors:
# "ENV['HOME'] is unset!"
# Only do so where there is no home, so normal users can use a different thing.
command = "rq #{ARGV.join(' ')}"
if !ENV['HOME']
  command = "export PATH=$PATH:/usr/local/bin; export HOME='/var/www'; #{command}"
end
system(command)
