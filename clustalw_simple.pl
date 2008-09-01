#!/usr/bin/perl -w

#CLUSTALW is annoying to use. This program takes 2 input files specified at the command
# line, concatenates them, and then feeds them into clustalw using the default options
# and the concatenated file

# I tried to use tempfile but it required XXXX's at the end of the name, and i need .fa at the end.
my $tmp_filename = "/tmp/clustalw_simple.tmp.".rand;
my $tmp_out_filename = "/tmp/clustalw_simple.tmp.".rand;


if ($#ARGV=>0){
  $args = join "\" \"", @ARGV;
  print `cat \"$args\" >$tmp_filename`;
} else {
  open TMP, ">$tmp_filename" or die "couldn't open temp file";
  foreach $line (<STDIN>){
    print TMP $line;
  }
  close TMP;
}
print `clustalw -infile=$tmp_filename -outfile=$tmp_out_filename`;
print `cat $tmp_out_filename`;


`rm $tmp_filename $tmp_out_filename`;
