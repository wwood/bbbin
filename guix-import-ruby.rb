#!/usr/bin/env ruby

require 'bio-commandeer'

out = Bio::Commandeer.run "guix import gem #{ARGV[0]}"
       
puts ['(define-public ruby-',
      ARGV[0],
      "\n",
      out.strip.gsub("expat","license:expat"),
      ')',
     "\n\n"].join
