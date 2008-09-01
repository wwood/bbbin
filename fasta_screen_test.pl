#!/usr/bin/perl -w

# Testing routines for fasta_screen.pl. fasta_screen.pl must be in the path

use Test::Simple qw(no_plan);
use TestingHelpers;

# test nothing
@out = ("ATG");
@err = ();
&datest('ATG',\@out,\@err);

# test simple
@out = ("ATG");
@err = ('g');
&datest('ATGg',\@out,\@err);


# 2 lines, both with ok chars
@out = ("ATG","GAC");
@err = ();
&datest("ATG\nGAC\n",\@out,\@err);

# 2 lines, both with bad chars
@out = ("ATG","GAC");
@err = ('c','fs');
&datest("cATG\nfsGAC\n",\@out,\@err);

# test empty line at end
print "Empty line\n";
@out = ("ATG","GAC",'','A');
@err = ('c','fs');
&datest("cATG\nfsGAC\n\nA",\@out,\@err);


# test starting > line
@out = ('>cATG',"GAC",'','A');
@err = ('fs');
&datest(">cATG\nfsGAC\n\nA",\@out,\@err);


print "test starting > line with spaces\n";
# test starting > line
@out = ('>cATG',"GAC",'A');
@err = ('fs','  >k');
&datest(">cATG\nfsGAC\n  >k\nA",\@out,\@err);



# test testing mode
@out = ();
@err = ('>cATG','fs','  >k');
&datest(">cATG\nfsGAC\n  >k\nA",\@out,\@err, '--t');



##test protein ones
#ok
@out = ('ACDEFGHIKLMNPQRSTUVWY');
@err = ();
&datest("ACDEFGHIKLMNPQRSTUVWY",\@out,\@err, '--protein');
#not ok
@out = ('ACDEFGHIKLMNPQRSTUVWY');
@err = ('X');
&datest("ACDEFGHIKLMNPQRSTUVWYX",\@out,\@err, '--protein');

##test protein and nx
@out = ('ACDEFGHIKLMNPQRSTUVWYX');
@err = ();
&datest("ACDEFGHIKLMNPQRSTUVWYX",\@out,\@err, '--protein --nx');



#test other
@out = ('ABC');
@err = ('D');
&datest("ABCD",\@out,\@err, '--other ABC');


sub datest {
  my ($contents,$out,$olderr,$args) = @_;

  # Add the annoying Removed the following characters: bit
  my @arr = ('Removed the following characters:');
  push @arr, @{$olderr};

  my $test = TestingHelpers->new();
  $test->cmdline_test('fasta_screen.pl', $args, $out, \@arr, $contents);
}
