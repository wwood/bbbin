package TestingHelpers;

use strict;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);

require Exporter;
use Test::Simple;

@ISA = qw(Exporter);
# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.
@EXPORT = qw(
	     $cmdline_test
	    );
$VERSION = '0.1';

sub new {
  my $package = shift;
  return bless({}, $package);
}

# Test what happens when certain a certain command line is run:
# Arguments:
# program - string of the program name, which must be in the path or fully specified.
# args = a string of the 
sub cmdline_test {
  my ($pkg, $program, $args, $r1, $r2, $contents) = @_;
  my @rout = @{$r1};
  my @rerr = @{$r2};
  if (!defined($args)){
    $args = '';
  }

  if ($contents){ open IN, ">in" or die "couldn't open temp file";
  print IN $contents; close IN; `$program $args >in.out 2>in.err <in`;
  } else { `$program $args >in.out 2>in.err`; }

  open OUT, "in.out";
  my @actual_out = <OUT>;
  close OUT;
  ok ($#actual_out == $#rout, '$#rout: '."$#rout vs $#actual_out");
  foreach my $i (0..$#rout){
    chomp $actual_out[$i];
    chomp $rout[$i];
    ok ($actual_out[$i] eq $rout[$i], "out line $i: $actual_out[$i] vs $rout[$i]");
  }

  open OUT, "in.err";
  my @actual_err = <OUT>;
  close OUT;
  ok ($#actual_err == $#rerr, '$#rerr: '.($#rerr)." vs $#actual_err");
  foreach my $i (0..$#rerr){
    if ($actual_err[$i]){
    chomp $actual_err[$i];
    chomp $rerr[$i];
    ok ($actual_err[$i] eq $rerr[$i], "err line $i: $actual_err[$i] vs $rerr[$i]");
  } else {
    ok (!defined(chomp $rerr[$i]), "err line $i: vs $rerr[$i] should be undef");
  }
  }
}





1; #stupid return true thing.
