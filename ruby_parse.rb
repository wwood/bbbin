#!/usr/bin/ruby

#require 'rubygems'
#gem 'rio'
#require 'rio'

stdin = $stdin.readlines.collect{|l| l.strip}

cmd = ARGV[0]
#puts cmd
eval(cmd)
