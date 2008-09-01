#!/usr/bin/perl -w

# Count how many distinct Match lines are left in the inspected (or not) output file


#V$FKHD/FHXB.01
#     V$FOXJ2_02	126	123	0.939	0.905	Fork head homologous X binds DNA with a dual sequence specificity (FHXA and FHXB)
#
#P$AHBP/HAHB4.01
#     V$HLF_01	135	129	0.938	0.957	Sunflower homeodomain leucine-zipper protein Hahb-4
#


foreach (<STDIN>) {
  if (m/^     (.*?\t.*?\t)/) {
    $hit = $1;
    $hit =~ s/\$//g;
    if (!(grep /$hit/, @hits)){
      push @hits, $hit;
    }
  }
}

$num = $#hits + 1;
print "$num\n";
#print "Distinct Match Records: $num\n";
