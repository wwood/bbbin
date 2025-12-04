# motif_rater

A utility for summarizing motif content in genome FASTA files.

## Usage

```bash
cargo run -- \
  --genome-fasta-file genome1.fa --genome-fasta-file genome2.fa \
  --motif ATGC --motif-file motifs.txt
```

Inputs mirror common CoverM FASTA options: provide genome files directly, supply a file with newline-delimited paths via `--genome-fasta-list`, or scan directories with `--genome-fasta-directory`.

Output columns (tab-separated) per genome:

1. `gc_content`: GC fraction across the genome.
2. `length`: total number of bases.
3. `{motif}_count`: occurrences of each motif plus its reverse complement.
4. `{motif}_expected`: expected occurrences given GC content and genome length.
