#!/usr/bin/perl -w

# Transforms a GFF output from exonerate and transforms
# it so it is suitable for a gff2ps plot.


#A typical exonerate est2genome output:
#NKCluster	exonerate:est2genome	gene	30540	31243	1990	+	.	gene_id 0 ; sequence CL2941546Contig1
# ; gene_orientation -
#NKCluster	exonerate:est2genome	utr5	30540	30712	.	+	.	
#NKCluster	exonerate:est2genome	exon	30540	30712	.	+	.	insertions 0 ; deletions 0
#NKCluster	exonerate:est2genome	splice3	30713	30714	.	+	.	intron_id -1 ; splice_site "CT"
#NKCluster	exonerate:est2genome	splice5	30771	30772	.	+	.	intron_id 1 ; splice_site "AC"
#NKCluster	exonerate:est2genome	utr5	30773	30862	.	+	.	
#NKCluster	exonerate:est2genome	exon	30773	30862	.	+	.	insertions 0 ; deletions 0
#NKCluster	exonerate:est2genome	splice3	30863	30864	.	+	.	intron_id -1 ; splice_site "CT"
#NKCluster	exonerate:est2genome	splice5	30973	30974	.	+	.	intron_id 1 ; splice_site "AC"
#NKCluster	exonerate:est2genome	utr5	30975	31022	.	+	.	
#NKCluster	exonerate:est2genome	exon	30975	31022	.	+	.	insertions 0 ; deletions 0
#NKCluster	exonerate:est2genome	splice3	31023	31024	.	+	.	intron_id -1 ; splice_site "CT"
#NKCluster	exonerate:est2genome	splice5	31075	31076	.	+	.	intron_id 1 ; splice_site "AC"
#NKCluster	exonerate:est2genome	utr5	31077	31167	.	+	.	
#NKCluster	exonerate:est2genome	exon	31077	31167	.	+	.	insertions 0 ; deletions 0
#NKCluster	exonerate:est2genome	splice3	31168	31169	.	+	.	intron_id -1 ; splice_site "CT"
#NKCluster	exonerate:est2genome	splice5	31222	31223	.	+	.	intron_id 1 ; splice_site "AC"
#NKCluster	exonerate:est2genome	exon	31224	31243	.	+	.	insertions 0 ; deletions 0
#NKCluster	exonerate:est2genome	similarity	30540	31243	1990	+	.	alignment_id 0 ; Query CL2941546Contig1



# What I want from this:
#NKCluster	exonerate:est2genome	utr5	30540	30712	.	+	.	sequence CL2941546Contig1
#NKCluster	exonerate:est2genome	exon	30540	30712	.	+	.	sequence CL2941546Contig1
#.. etc. only including utr5 and exon things.



# For some reason bioperl flunks out on these gff's, so I'll just
# write a parser from scratch. I have an old version or something?



my $cur_gene_name = undef;
while (my $feature = <STDIN>){
  @bits = split /\t/, $feature;
  #print STDERR $#bits."\n";
  if ($#bits == 0 or !($bits[1] =~ m/exonerate\:/)){
    #print STDERR "ignoring line: $feature";
    next;
  }

  if ($bits[2] eq 'gene'){
    $bits[8] =~ m/sequence (.*?) /;
    $cur_gene_name = $1;
    $cur_gene_name =~ s/ESTCONTIG\:CLU002\://;
    $cur_gene_name =~ s/\|.*//;
  }

  # I haven't observed utr3, but just in case it is included here.
  elsif ($bits[2] eq 'exon'){# or $bits[2] eq 'utr5' or $bits[2] eq 'utr3'){
    $first = join "\t", @bits[0,1,2,3,4,5,6,7];
    print $first."\t$cur_gene_name\n";
  }

  elsif ($bits[2] eq 'similarity'){# or $bits[2] eq 'utr5' or $bits[2] eq 'utr3'){
    print "$feature";
  }
}




