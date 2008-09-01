# Testing procedure for sequenceChop.pl


use Test::Simple qw(no_plan);
use TestingHelpers;


#   1. A fasta file with 10 elements, chopping out 2-5
@out = ('>yey_chopped_2-5','TGCA');
@err = ();
&datest(">yey\nATGCATGCAT\n",2,5,\@out,\@err);

#   2. A fasta file with 10 elements, chopping out 1-5
@out = ('>yey_chopped_1-5','ATGCA');
@err = ();
&datest(">yey\nATGCATGCAT\n",1,5,\@out,\@err);

#   3. A fasta file with 10 elements, chopping out 2-10
@out = ('>yey_chopped_2-10','TGCATGCAT');
@err = ();
&datest(">yey\nATGCATGCAT\n",2,10,\@out,\@err);

#   4. A fasta file with 2 elements and a space in the description line, chopping out 1-1
@out = ('>yey desc_chopped_1-1','A');
@err = ();
&datest(">yey desc\nATGCATGCAT\n",1,1,\@out,\@err);

#   5. Using negative positions
@out = ('>yey desc_chopped_8-9','CA');
@err = ();
&datest(">yey desc\nATGCATGCAT\n",-3,-2,\@out,\@err);

#   6. Using multiple input files
@out = ('>yey desc_chopped_8-9','CA','>yey desc2_chopped_8-9','CC');
@err = ();
&datest(">yey desc\nATGCATGCAT\n>yey desc2\nATGCATGCCTCCT\n",-3,-2,\@out,\@err);

sub datest {
  my ($fasta, $start, $stop, $out, $err) = @_;

  my $test = TestingHelpers->new();
  $test->cmdline_test('sequenceChop.pl', '- '.$start.' '.$stop, $out, $err, $fasta);
}
