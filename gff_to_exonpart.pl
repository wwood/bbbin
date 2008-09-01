#!/usr/bin/perl -w

# Convert vulgar into a exonpart for use by AUGUSTUS as a hint

# Input: lines like:
##gff-version 2
##source-version exonerate:protein2genome:local 2.0.0
##date 2008-03-09
##type DNA
#
#
# seqname source feature start end score strand frame attributes
#
# Contig13185exonerate:protein2genome:localgene2750928753757-.gene_id 1 ; sequence jgi|Thaps3|644|fgenesh1_pm.C_chr_9000012 ; gene_orientation +
#  Contig13185exonerate:protein2genome:localcds2794928753.-.
#  Contig13185exonerate:protein2genome:localexon2794928753.-.insertions 6 ; deletions 3 ; frameshifts 1
#  Contig13185exonerate:protein2genome:localsplice52794727948.-.intron_id 1 ; splice_site "GT"
#  Contig13185exonerate:protein2genome:localintron2768627948.-.intron_id 1
#  Contig13185exonerate:protein2genome:localsplice32768627687.-.intron_id 0 ; splice_site "AG"
#  Contig13185exonerate:protein2genome:localcds2750927685.-.
#  Contig13185exonerate:protein2genome:localexon2750927685.-.insertions 0 ; deletions 0
#  Contig13185exonerate:protein2genome:localsimilarity2750928753757-.alignment_id 1 ; Query jgi|Thaps3|644|fgenesh1_pm.C_chr_9000012 ; Align 28754 1 75 ; Align 28678 26 72 ; Align 28603 50 60 ; Align 28543 73 471 ; Align 28069 230 120 ; Align 27686 270 177
# --- END OF GFF DUMP ---



# Output: lines like: (1 to 1 ratio), formatted like GFF
# Contig13185    anchor    exonpart   27949  28753  0   -   .   source=M
# Contig13185    anchor    exonpart   27509  27685  0   -   .   source=M




foreach (<>){
    chomp;
    my @splits = split /\t/, $_;

    if ($#splits != 7 and $#splits != 8){  next; }

    if ($splits[2] eq 'cds'){
      print "$splits[0]\tanchor\texonpart\t$splits[3]\t$splits[4]\t.\t$splits[6]\t.\tsource=M\n";
    } elsif ($splits[2] eq 'splice5'){
      print "$splits[0]\tanchor\tdss\t$splits[3]\t$splits[3]\t.\t$splits[6]\t.\tsource=M\n";
    } elsif ($splits[2] eq 'splice3'){
      print "$splits[0]\tanchor\tass\t$splits[4]\t$splits[4]\t.\t$splits[6]\t.\tsource=M\n";
    }
}
