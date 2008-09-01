#!/usr/bin/perl -w

# Screen a fasta file to remove all the crap, excluding all lines that start with a '>'


use Getopt::Long;

GetOptions('i=s'=>\$input,
	   'test'=>\$test,
	   'dna'=>\$dna_ok,
	   'protein'=>\$protein_ok,
	   'nx'=>\$nx_ok,
	   'other=s'=>\$other_chars);#,
#	   'underscore'=>\$underscore);

my $dna_chars = 'ATGCRYNVKMSW';
my $protein_chars = 'ACDEFGHIKLMNPQRSTUVWY*BJUOZ';
my $nx_chars = 'NX';



my $ok_chars = '';
$ok_chars .= $dna_chars if $dna_ok;
$ok_chars .= $protein_chars if $protein_ok;
$ok_chars .= $nx_chars if $nx_ok;
$ok_chars .= $other_chars if $other_chars;

#by default, protein and dna are OK
if ($ok_chars eq ''){
  $ok_chars .= $dna_chars;
  $ok_chars .= $protein_chars;
  $ok_chars .= $nx_chars;
}
#print $ok_chars;

my $fh = STDIN;
if (defined($input)){
  open IN, $input or die "could not open file to parse it: $input.";
  $fh = IN;
} elsif ($#ARGV == 0){
  open IN, $ARGV[0] or die "could not open file to parse it: $input.";
  $fh = IN;
}

print STDERR "Removed the following characters:\n";
foreach my $l (<$fh>){
  chomp $l;
  if ($l =~ s/(^>.*)//){
    if (!$test){
      print $1."\n";
    }
  } else {
    # I don't know how to print out each removed character otherwise,
    # so I have to do some looping.
    if ($l =~ m/[^$ok_chars]+/i){
      while ($l =~ s/[^$ok_chars]+//i){
	print STDERR $&;
      }
      print STDERR "\n";
    }
    if (!$test){
      print $l."\n";
    }
  }
}


