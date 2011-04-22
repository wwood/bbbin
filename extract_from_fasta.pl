#!/usr/bin/perl -w

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Std;


our ($opt_r, $opt_f, $opt_g, $opt_v, $opt_F, $opt_s, $opt_i);
getopts('rf:gvFs:i');


# given a fasta file with multiple sequences, extract a single sequence

if ($#ARGV != 0 && $#ARGV != 1)
  {
    print STDERR "usage: extract_from_fasta <seq_name>\n";
    print STDERR "the input fasta file is piped in, or used as the second argument.\n\n";
    print STDERR "Options:\n";
    print STDERR "-r use regular expressions in the names of the sequences\n";
    print STDERR "-f extract multiple sequences. Sequence names/regular expressions are defined in\n";
    print STDERR "   a newline separated list file\n";
    print STDERR "-v reverse. Print out sequences that don't match (cf grep -v)\n";
    print STDERR "-F first only. Match only the first part of the fasta name, i.e. the everything before the first space character or fasta name.\n";
    print STDERR "-s sequence. Match using the sequence of an input fasta file, not the IDs at all. Incompatible with everything else, except -v\n";
    print STDERR "-g input a gff file with the names of the sequences in the first col, \n";
    print STDERR "   and chop out the sequence up so with the coordinates given in the start/stop\n";
    print STDERR "   columns, except with 2kb on each side.\n";
    print STDERR "-i Ignore case when matching (currently only works when using regular expressions (-r))\n";
    print STDERR "\n";
    exit;
  }

# For regular, normal matches given on the command line.
my $seqname = undef;
unless ($opt_f or $opt_s){
  $seqname = shift @ARGV;
}

my @lines = ();
my $seqio = '';
if ($#ARGV==-1){
  $seqio = Bio::SeqIO->new(-fh => \*STDIN, -format => 'Fasta');
} else {
  $seqio = Bio::SeqIO->new(-file => shift @ARGV, -format => 'Fasta');
}



# Check for incompatible arguments
if ($opt_s){
  if ($opt_f or $opt_r or $opt_g or $opt_F){
    print STDERR "-s is incompatible with -f, -r, -F and -g\n";
    exit 1;
  }
}



# If an input file of sequence names is specified, use that instead
# Create a hash of these so that it can be accessed quickly later.
my %seqnames = ();
my @seqnames = ();
if ($opt_f){
  open NAMES, $opt_f or die "couldn't open file with gene names: '$opt_f'";
  foreach(<NAMES>){
    chomp;
    if (s/[^[:print:]]+//g){
      print STDERR "WARNING: non-standard characters found in the name: '$_'. These were removed before extraction.\n"
    }
    if (defined($seqnames{$_})){
      $seqnames{$_} += 1;
    } else {
      $seqnames{$_} = 1;
      if ($opt_r){
	push @seqnames, $_;
      }
    }
  }
}


# This next stanza of code concerns the -s flag
my %seq_to_id_hash = (); #hash of fasta sequences to IDs, so reporting problems when -s is specified is more coherent
if ($opt_s){
  # Read in fasta file of probe sequences
  $input_seqio = Bio::SeqIO->new(-file => $opt_s, -format => 'Fasta');
  while (my $seq = $input_seqio->next_seq){
    # Save the seq as a key in %seq_to_id_hash, where the value is the ID
    my $displayId = $seq->display_id;
    if ($seq->desc){
      $displayId .= ' '.$seq->desc;
    }

    my $new_seq = $seq->seq;
    $new_seq =~ s/\*$//; #remove trailing stop codon if one exists
    $seq_to_id_hash{$new_seq} = $displayId;
    $seqnames{$displayId}++;
  }
}






my $any_matches = 0; # any matches at all?
while (my $seq = $seqio->next_seq){
  # Record whether it matches or not, so we can -v it if necessary
  my $matches = 0;

  # I want the whole of the line, not just the display_id.
  # Annoying how there is no function in bioperl to do this.
  # That is, I want it unless -F (first only) is specified.
  my $displayId = $seq->display_id;
  unless ($opt_F){
    if ($seq->desc){
      $displayId .= ' '.$seq->desc;
    }
  }

  if ($opt_f){ #file of exact matches
    if ($opt_r){
      # Have to go through each of the entries in the array each 
      # time, making this a slow operation. Is it possible to apply
      # many regular expressions on a single string all at once?
      foreach $name (@seqnames){
	my $answer = undef;
	if ($opt_i){
	  $answer = ($displayId =~ m/$name/i);
	} else {
	  $answer = ($displayId =~ m/$name/);
	}
	if ($answer){
	  $matches = 1;
	  $seqnames{$name}++;
	  last; #only print the sequence once.
	}
      }
    }

    # No regular expression - easy.
    else {
      if ($seqnames{$displayId}){
	$matches = 1;
	$seqnames{$displayId}++;
      }
    }
  } elsif ($opt_r){#regular expression match
    my $answer = undef;
    if ($opt_i){
      $answer = ($displayId =~ m/$seqname/i);
    } else {
      $answer = ($displayId =~ m/$seqname/);
    }
    if ($answer){
      $matches = 1;
      $any_matches = 1;
    }
  } elsif ($opt_s){#sequence match
    # remove trailing stop codon
    my $new_seq = $seq->seq;
    $new_seq =~ s/\*$//;
    if ($seq_to_id_hash{$new_seq}){
      $matches = 1;
      $any_matches = 1;
      $seqnames{$seq_to_id_hash{$new_seq}}++;
    }
  } elsif ($displayId eq $seqname){# exact match
    $matches = 1;
    $any_matches = 1;
  }


  # if normal and matches, print the sequence
  # elsif inverse and not matches, print it
  if ((!$opt_v && $matches) || ($opt_v && !$matches)){
    &print_seq($seq);
  }
}



# If there was a file used, check to make sure each of the sequences
# specified were actually found in the file.
if ($opt_f or $opt_s){
  while (my ($key, $value) = each(%seqnames)){
    # $value by default is 1 because otherwise it isn't true in an if
    # statement.
    if ($value == 1){
      print STDERR "WARNING: sequence '$key' not found in the fasta\n";
    } elsif ($value > 2){
      if ($opt_r){
	print STDERR "WARNING: sequence '$key' included ".($value-1)." times in the file. Does it exist twice in the fasta file?, or did it match that many different regular expressions?\n";
      } else {
	print STDERR "WARNING: sequence '$key' included ".($value-1)." times in the file. Does it exist twice in the fasta file?\n";
      }
    }
  }
}
elsif (!$any_matches){
  print STDERR "WARNING: sequence '$seqname' not found in the fasta\n";
}

sub print_seq {
  my $seq = $_[0];
  print '>'.$seq->display_id.' '.$seq->description."\n".$seq->seq."\n";
}

