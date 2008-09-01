#!/usr/bin/perl -w





#INPUT EXAMPLE
#"#Contig1.1 5799-3300,2804-2042"
#1913
#2016
#"2198-2199"
#2212
#2277
#2280
#2283
#2295
#3363
#5784
#5803
#5928
#6106
#6788
#6789
#6808
#6821
#6823
#6833
#6843
#6850
#6855
#6861
#6923


#"#Contig2.1 4460-1965"
#"between 2015 2016"
#"3173-3715"
#"3716-3186"
#"between 3186 3187"
#"3187-3191"
#3264
#3277
#3500
#3501
#3502
#3606
#3637
#3642
#"3716-3186"
#3854
#4338
#"4341-4344"
#4386
#4401
#"between 4585 4586




$usage = "$0 <ManualSNPsFile>\n\n";
open IN, "$ARGV[0]" or die "$usage";


@all_lines = <IN>;
@all_lines =~ s/\"//g;

$everything = join "\n",@all_lines;
@contigs = split '#',$everything;

foreach $contig (@contigs){
  @lines = split "\n",$contig;
  
  $started = 0;
  foreach $line (@lines){
    #skip blank lines
    if ($line =~ m/^\s$/){next;}

    #if the line starts with a '#' then we are starting a new contig.
    if ($line =~ m/^#/){
      if ($starting==1){die "badly parsed at line: $line";}
      $started = 1;

      $line =~ m/^#(\S+) (\S+)\s*/ or die "badly parsed at line: $line";
      $name = $1;
      $positions = $2;
      
      #parse position; it may contain >1 area
      if ($positions =~ m/(\d+)-(\d+),(\d+)-(\d+)/){
	$start = $1;
	$finish = $4;
      } elsif ($positions =~ m/(\d+)-(\d+)/){
	$start = $1;
	$finish = $2;
      } else {
	die "badly parsed at line: $line";
      }


      open OUT, "SNP_plot_$name.txt";
    }


    # else it is a normal line
    else {
      if ($line =~ m/between (\d+).*?,([atgc\*]+)/)
	push @betweens, 
    }


  }
  close OUT;
}
