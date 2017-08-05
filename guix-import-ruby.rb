#!/usr/bin/env ruby

require 'bio-commandeer'

out = Bio::Commandeer.run "~/git/guix/pre-inst-env guix import gem #{ARGV[0]}"
       
puts ["\n\n(define-public ruby-",
      ARGV[0],
      "\n",
      out.strip,
      ')',
     "\n\n"].join
