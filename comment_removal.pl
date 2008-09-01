#!/usr/bin/perl -w

# Remove all the commented lines from a file

foreach (<>){
  if (!m/^\s*\#/ and !m/^\s*$/){
    print $_;
  }
}
