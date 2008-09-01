#!/usr/bin/perl -w

# Takes in a fasta file from STDIN, and then genscans each of the sequences in that fasta file


use Getopt::Long;
use Bio::SeqIO;

# A default home directory was added when making this an encapsulated program
$PAR_HOME = $ENV{HOME}.'/bioinfo/genscan';
$TEMP_DIR = '/tmp/genscanMultiFasta';

$parfname = "HumanIso.smat";
$postscript = 'unspecified';
GetOptions('pdf:s'=> \$postscript,
	   'parfname=s' => \$parfname,
	   'help'=>\$help,
	   'subopt=s'=>\$subopt,
	   'v'=>\$verbose,
	   'cds'=>\$cds,
	   'fasta=s'=>\$input);
if ($postscript eq ''){
  $postscript = 'default';
}
$parfname = "$PAR_HOME/$parfname";

#print 'm'.$postscript.".\n";
#exit;

if ($help){
  print "Usage: $0 [options]\n".
    "\tFasta file is piped in\n".
    "\tWhere [options] is:\n".
      "\t--help\t\tDisplay this msg\n".
	"\t--pdf[=groupID]\tCreate pdf outputs with genscan's ps output and ps2pdf. groupID is an optional argument inserted into the output filenames.\n".
	  "\t--parfname=file\tUse file as the parameter file.\n".
"\t--fasta=input\t specify the input file directly, rather than through stdin.\n".
	    "\n";
  exit;
}

if ($#ARGV == 0){
    $input = $ARGV[0];
}


if(defined($input)){
  print "input: $input\n";
  $original = Bio::SeqIO->new(-file => $input,
			    '-format' => 'fasta');
} else {
  $original = Bio::SeqIO->new(-fh => \*STDIN,
			    '-format' => 'fasta');
}


if (!(-d($TEMP_DIR))){
  `mkdir $TEMP_DIR`;
}
$tempfile = $TEMP_DIR.'/input-'.rand().".fasta";

`touch $tempfile`;
while ($seq = $original->next_seq()){
  $tmp_out = Bio::SeqIO->new(-file => ">$tempfile",
			     '-format' => 'fasta');
  $tmp_out->write_seq($seq);



  $args = "";
  #print "ps:$postscript.\n";  
  if (!($postscript eq 'unspecified')){
    print "ps:$postscript.\n";
    $id = $seq->display_id();
    $name = "";
    if (!($postscript eq 'default')){
      $name .= $postscript.'_';
    }
    $name .= "$id".'.ps';
    $args = " -ps ".$name;
  }
  if ($verbose){
    $args = "$args -v";
  }
  if ($cds){
    $args = "$args -cds"
  }
  #print "args:".$args.".\n";

  $suboptArg = '';
  if ($subopt){
    $suboptArg = '-subopt '.$subopt;
  }

  #print "genscan $parfname $tempfile $suboptArg $args"."\n";
  print `genscan $parfname $tempfile $suboptArg $args`;
  print '****************************************************************'."\n";

  # Convert ps to pdf if ps was generated
  if (!($postscript eq 'unspecified')){
    `ps2pdf $name`;
    #`rm $name`;
  }
}
