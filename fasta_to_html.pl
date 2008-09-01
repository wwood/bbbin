#!/usr/bin/perl -w

# Takes a fasta file and creates a directory and tree of individual files so they can be accessed.


use Bio::SeqIO;



if ($#ARGV != 0){
  print STDERR "Usage: $0 <fasta>\n\n";
  exit;
}


#create the fasta file and directory
open INDEX, ">$ARGV[0].html" or die "couldn't open index file";

my $DIR_NAME = "$ARGV[0].individuals";		# name of the directory that holds the individual files
`mkdir $DIR_NAME`;


my $in  = Bio::SeqIO->new(-file => $ARGV[0],
			  -format => 'Fasta');


while ( my $seq = $in->next_seq() ) {
  #use the whole name to be safe
  my $whole = $seq->id.' '.$seq->description;

  my $fa_name = "$DIR_NAME/$whole.html";
  open FA, ">$fa_name" or die "couldn't open '$fa_name'";
  print FA "<pre>\n>".$seq->id.' '.$seq->description."\n".$seq->seq."\n</pre>";
  close FA;

  print INDEX "<a href='$fa_name'>$whole</a><br />\n";
}

close INDEX;
