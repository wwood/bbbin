#!/usr/bin/perl -w

# generate the thesis appendix for all the code

open APPENDIX, ">$ENV{HOME}/thesis/doc/thesis/codeAppendix.tex" or die "fuck";


print APPENDIX '\chapter{Code Listing}';

print APPENDIX 'Included here is all the code used in this study. It is instructive to know some general principles used. Firstly, many scripts use unix piping (STDIN, STDERR \& STDOUT), to facilitate the creation of pipelines. Code was for the most part written using bash shell scripting, gnuplot scripting, and Perl. A significant number of the Perl scripts rely on the useful BioPerl modules to interact with biological file formats such as ACE and FASTA. Code listing is alphabetical.'."\n\n\n";




print APPENDIX '\section{Code}'."\n\n";

print APPENDIX '\\begin{tiny}'."\n\n\n";

$bin = "$ENV{HOME}/thesis/bin" or die "binbad";
opendir DIR, "$bin";

chdir $bin;

foreach $file (readdir DIR) {

  if ($file =~ m/~/ || $file=~m/\.\./ || $file=~/^\./ ||$file=~m/^consed$/ ||$file=~m/^testing$/ ||$file=~m/\#/) {
    print "ignored: $file\n";
    next;
  } else {
    $n = $file;
    $n =~ s/_/\\_/g;
    print APPENDIX '\subsection{'.$n.'}';
    print APPENDIX '\verbatiminput{\home/thesis/bin/'.$file."}\n";
    #	print '\verbatiminclude{'.$file."}\n";
  }
}


print APPENDIX '\\end{tiny}'."\n";
