#!/usr/bin/perl -w

#<A HREF="/pub/databases/transfac/doc/factor2.html#AC" target="_new">AC</A>   T00902
#<A HREF="/pub/databases/transfac/doc/ class1.html#AC" target="_new">AC</A>   C0008


opendir DIRECTORY, ".";

while (defined($file = readdir DIRECTORY)) {
  open HTML, $file;

  if ($file =~ m/^\.$/ || $file =~ m/^\.\.$/){
    next;
  }

  #remove everything
  while (defined(pop @matches)) {
    ;
  }


  foreach $line (<HTML>) {
    chomp $line;
    if ($line =~ m/<A HREF=\"\/pub\/databases\/transfac\/doc\/.*.html\#AC\" target=\"_new\">AC<\/A>\s*(\w*)/) {
      push @matches, $1;
    }
  }

  if ($#matches == 0){
    `cp $file ../parsed/$matches[0].html`;
  }
  else {
    print "$#matches: $file\n";
  }
	
}
