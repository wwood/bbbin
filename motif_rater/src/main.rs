use anyhow::{anyhow, Context, Result};
use clap::{ArgAction, Parser};
use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
#[command(author, version, about = "Compute motif statistics for genomes")]
struct Args {
    /// Individual FASTA files containing genomes.
    #[arg(short = 'f', long = "genome-fasta-file")]
    genome_fasta_files: Vec<PathBuf>,

    /// Directories to scan for FASTA files (non-recursive).
    #[arg(short = 'd', long = "genome-fasta-directory")]
    genome_fasta_directories: Vec<PathBuf>,

    /// Files listing FASTA paths, one per line.
    #[arg(short = 'l', long = "genome-fasta-list")]
    genome_fasta_lists: Vec<PathBuf>,

    /// Motifs provided on the command line (IUPAC notation).
    #[arg(short = 'm', long = "motif")]
    motifs: Vec<String>,

    /// Files containing motifs, one per line.
    #[arg(long = "motif-file")]
    motif_files: Vec<PathBuf>,

    /// Print header row in output.
    #[arg(long, default_value_t = true, action = ArgAction::Set, value_name = "BOOL")]
    header: bool,
}

#[derive(Clone)]
struct MotifPattern {
    label: String,
    reverse_complement: String,
    allowed_forward: Vec<Vec<u8>>,
    allowed_reverse: Vec<Vec<u8>>,
}

struct GenomeMetrics {
    name: String,
    length: u64,
    gc_count: u64,
    motif_counts: Vec<usize>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let motifs = collect_motifs(&args).context("failed to collect motifs")?;
    if motifs.is_empty() {
        return Err(anyhow!(
            "provide at least one motif using --motif or --motif-file"
        ));
    }

    let genomes = collect_genome_paths(&args).context("failed to collect genome FASTA files")?;
    if genomes.is_empty() {
        return Err(anyhow!(
            "provide at least one genome FASTA via --genome-fasta-file, --genome-fasta-directory, or --genome-fasta-list"
        ));
    }

    if args.header {
        print_header(&motifs);
    }

    for genome in genomes {
        let metrics = analyze_genome(&genome, &motifs)
            .with_context(|| format!("processing genome {:?}", genome))?;
        output_metrics(&metrics, &motifs);
    }

    Ok(())
}

/// Gather motif definitions from CLI flags and files, removing empties and duplicates.
fn collect_motifs(args: &Args) -> Result<Vec<MotifPattern>> {
    let mut motifs: Vec<String> = args
        .motifs
        .iter()
        .map(|m| m.trim().to_string())
        .filter(|m| !m.is_empty())
        .collect();

    for file in &args.motif_files {
        motifs.extend(read_lines(file)?.into_iter().map(|l| l.trim().to_string()))
    }

    motifs.retain(|m| !m.is_empty());

    motifs
        .into_iter()
        .map(|motif| MotifPattern::new(&motif))
        .collect()
}

/// Aggregate genome FASTA paths from explicit files, directory listings, and list files.
fn collect_genome_paths(args: &Args) -> Result<Vec<PathBuf>> {
    let mut paths = BTreeSet::new();

    for file in &args.genome_fasta_files {
        paths.insert(file.clone());
    }

    for list in &args.genome_fasta_lists {
        for line in read_lines(list)? {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            paths.insert(PathBuf::from(trimmed));
        }
    }

    for dir in &args.genome_fasta_directories {
        let entries = fs::read_dir(dir).with_context(|| format!("reading directory {:?}", dir))?;
        for entry in entries {
            let entry = entry?;
            if entry.file_type()?.is_file() {
                let path = entry.path();
                if is_fasta(&path) {
                    paths.insert(path);
                }
            }
        }
    }

    Ok(paths.into_iter().collect())
}

/// Parse a FASTA file and aggregate metrics across all sequences it contains.
fn analyze_genome(path: &Path, motifs: &[MotifPattern]) -> Result<GenomeMetrics> {
    let file = File::open(path).with_context(|| format!("opening {:?}", path))?;
    let reader = BufReader::new(file);

    let mut sequence = String::new();
    let mut length: u64 = 0;
    let mut gc_count: u64 = 0;
    let mut motif_counts = vec![0usize; motifs.len()];

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !sequence.is_empty() {
                process_sequence(
                    &sequence,
                    &mut length,
                    &mut gc_count,
                    motifs,
                    &mut motif_counts,
                );
                sequence.clear();
            }
        } else {
            sequence.push_str(line.trim());
        }
    }

    if !sequence.is_empty() {
        process_sequence(
            &sequence,
            &mut length,
            &mut gc_count,
            motifs,
            &mut motif_counts,
        );
    }

    Ok(GenomeMetrics {
        name: path
            .file_name()
            .and_then(|os| os.to_str())
            .unwrap_or("unknown")
            .to_string(),
        length,
        gc_count,
        motif_counts,
    })
}

/// Update running totals for a single FASTA sequence.
///
/// The function normalizes `U` to `T`, increments total length and GC counts,
/// and tallies matches for each motif on the forward and reverse strands.
fn process_sequence(
    seq: &str,
    length: &mut u64,
    gc_count: &mut u64,
    motifs: &[MotifPattern],
    motif_counts: &mut [usize],
) {
    let normalized: Vec<u8> = seq
        .as_bytes()
        .iter()
        .map(|b| match b.to_ascii_uppercase() {
            b'U' => b'T',
            other => other,
        })
        .collect();

    *length += normalized.len() as u64;
    *gc_count += normalized
        .iter()
        .filter(|b| **b == b'G' || **b == b'C')
        .count() as u64;

    for (motif_index, motif) in motifs.iter().enumerate() {
        let forward = count_matches(&normalized, &motif.allowed_forward);
        let reverse = count_matches(&normalized, &motif.allowed_reverse);

        motif_counts[motif_index] += forward + reverse;
    }
}

/// Count how many windows in `sequence` satisfy the set of allowed bases per position.
fn count_matches(sequence: &[u8], allowed: &[Vec<u8>]) -> usize {
    if allowed.is_empty() || sequence.len() < allowed.len() {
        return 0;
    }

    let motif_len = allowed.len();
    let mut count = 0usize;

    for window_start in 0..=sequence.len() - motif_len {
        let mut matched = true;
        for (offset, allowed_bases) in allowed.iter().enumerate() {
            let base = sequence[window_start + offset];
            if !allowed_bases.iter().any(|candidate| *candidate == base) {
                matched = false;
                break;
            }
        }
        if matched {
            count += 1;
        }
    }

    count
}

/// Emit a tab-separated line containing genome-level metrics and motif statistics.
fn output_metrics(metrics: &GenomeMetrics, motifs: &[MotifPattern]) {
    let gc_fraction = if metrics.length > 0 {
        metrics.gc_count as f64 / metrics.length as f64
    } else {
        0.0
    };

    print!("{}\t{}\t{:.6}", metrics.name, metrics.length, gc_fraction);

    for (idx, motif) in motifs.iter().enumerate() {
        let expected = expected_motif_count(motif, gc_fraction, metrics.length);
        let windows = motif_window_count(metrics.length, motif.allowed_forward.len());
        let observed_rate = if windows > 0 {
            metrics.motif_counts[idx] as f64 / windows as f64
        } else {
            0.0
        };

        print!(
            "\t{}\t{:.3}\t{:.6}",
            metrics.motif_counts[idx], expected, observed_rate
        );
    }

    println!();
}

/// Compute the expected number of motif occurrences given the GC content and sequence length.
///
/// The calculation multiplies the probability of the motif appearing at any position by the
/// number of possible windows. Forward and reverse-complement probabilities are summed for
/// every motif, even when the motif is palindromic, because both orientations are reported
/// independently. For example, the motif `ACG` has probability 0.25³ = 0.015625 on a sequence
/// with GC content 0.5; its reverse complement `CGT` has the same probability. In a 100 bp
/// sequence there are 98 windows of length 3, so the expected total is (0.015625 × 2) × 98 ≈
/// 6.125 occurrences.
fn expected_motif_count(motif: &MotifPattern, gc_fraction: f64, length: u64) -> f64 {
    if motif.allowed_forward.is_empty()
        || length == 0
        || length < motif.allowed_forward.len() as u64
    {
        return 0.0;
    }

    let prob_forward = motif_probability(&motif.allowed_forward, gc_fraction);
    let prob_reverse = motif_probability(&motif.allowed_reverse, gc_fraction);

    let combined_prob = prob_forward + prob_reverse;

    let positions = motif_window_count(length, motif.allowed_forward.len());
    combined_prob * positions as f64
}

/// Probability that a random window with the given GC content matches the motif constraints.
fn motif_probability(allowed: &[Vec<u8>], gc_fraction: f64) -> f64 {
    allowed
        .iter()
        .map(|bases| {
            bases
                .iter()
                .map(|b| base_probability(*b as char, gc_fraction))
                .sum::<f64>()
        })
        .product()
}

/// Probability of observing a base given the GC fraction of the sequence.
fn base_probability(base: char, gc_fraction: f64) -> f64 {
    match base {
        'G' | 'C' => gc_fraction / 2.0,
        'A' | 'T' => (1.0 - gc_fraction) / 2.0,
        _ => 0.0,
    }
}

impl MotifPattern {
    /// Construct the forward and reverse-complement representations of a motif.
    fn new(raw: &str) -> Result<Self> {
        let label = raw.trim().to_ascii_uppercase();
        if label.is_empty() {
            return Err(anyhow!("motif cannot be empty"));
        }

        let allowed_forward = motif_to_allowed(&label)?;
        let reverse_complement = reverse_complement(&label)?;
        let allowed_reverse = motif_to_allowed(&reverse_complement)?;

        Ok(Self {
            label,
            reverse_complement,
            allowed_forward,
            allowed_reverse,
        })
    }
}

/// Convert an IUPAC motif string to allowed bases for each position.
fn motif_to_allowed(motif: &str) -> Result<Vec<Vec<u8>>> {
    motif
        .chars()
        .map(|c| iupac_bases(c).map(|bases| bases.iter().map(|b| *b as u8).collect()))
        .collect::<Option<Vec<Vec<u8>>>>()
        .ok_or_else(|| anyhow!("unsupported IUPAC symbol in motif: {}", motif))
}

fn iupac_bases(symbol: char) -> Option<&'static [char]> {
    match symbol.to_ascii_uppercase() {
        'A' => Some(&['A']),
        'C' => Some(&['C']),
        'G' => Some(&['G']),
        'T' | 'U' => Some(&['T']),
        'R' => Some(&['A', 'G']),
        'Y' => Some(&['C', 'T']),
        'S' => Some(&['G', 'C']),
        'W' => Some(&['A', 'T']),
        'K' => Some(&['G', 'T']),
        'M' => Some(&['A', 'C']),
        'B' => Some(&['C', 'G', 'T']),
        'D' => Some(&['A', 'G', 'T']),
        'H' => Some(&['A', 'C', 'T']),
        'V' => Some(&['A', 'C', 'G']),
        'N' => Some(&['A', 'C', 'G', 'T']),
        _ => None,
    }
}

/// Generate the reverse-complement of an IUPAC motif string.
fn reverse_complement(motif: &str) -> Result<String> {
    motif
        .chars()
        .rev()
        .map(|c| complement_symbol(c).ok_or_else(|| anyhow!("unsupported IUPAC symbol: {}", c)))
        .collect()
}

fn complement_symbol(symbol: char) -> Option<char> {
    match symbol.to_ascii_uppercase() {
        'A' => Some('T'),
        'C' => Some('G'),
        'G' => Some('C'),
        'T' | 'U' => Some('A'),
        'R' => Some('Y'),
        'Y' => Some('R'),
        'S' => Some('S'),
        'W' => Some('W'),
        'K' => Some('M'),
        'M' => Some('K'),
        'B' => Some('V'),
        'D' => Some('H'),
        'H' => Some('D'),
        'V' => Some('B'),
        'N' => Some('N'),
        _ => None,
    }
}

fn print_header(motifs: &[MotifPattern]) {
    print!("genome\tlength\tgc_content");
    for motif in motifs {
        print!(
            "\t{}_count\t{}_expected\t{}_observed_rate",
            motif.label, motif.label, motif.label
        );
    }
    println!();
}

fn motif_window_count(length: u64, motif_len: usize) -> u64 {
    if motif_len == 0 || length < motif_len as u64 {
        0
    } else {
        length - motif_len as u64 + 1
    }
}

fn is_fasta(path: &Path) -> bool {
    match path.extension().and_then(|e| e.to_str()) {
        Some(ext) => matches!(ext.to_ascii_lowercase().as_str(), "fa" | "fasta" | "fna"),
        None => false,
    }
}

fn read_lines(path: &Path) -> Result<Vec<String>> {
    let file = File::open(path).with_context(|| format!("opening {:?}", path))?;
    let reader = BufReader::new(file);
    let mut lines = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.trim_start().starts_with('#') {
            continue;
        }
        lines.push(line);
    }
    Ok(lines)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_motif_and_reverse() {
        let motif = MotifPattern::new("ACG").unwrap();
        let mut length = 0;
        let mut gc = 0;
        let mut counts = vec![0];
        process_sequence(
            "ACGTACG",
            &mut length,
            &mut gc,
            &[motif.clone()],
            &mut counts,
        );
        assert_eq!(counts[0], 3);
    }

    #[test]
    fn expected_counts_double_for_non_palindromic() {
        let motif = MotifPattern::new("ACG").unwrap();
        let gc = 0.5;
        let expected_single = motif_probability(&motif.allowed_forward, gc);
        let expected_total = expected_motif_count(&motif, gc, 1000);
        let positions = 1000 - motif.allowed_forward.len() as u64 + 1;
        assert!((expected_total - expected_single * 2.0 * positions as f64).abs() < 1e-6);
    }

    #[test]
    fn palindromic_expected_counts_are_doubled() {
        let motif = MotifPattern::new("ATAT").unwrap();
        let gc = 0.5;
        let expected_single = motif_probability(&motif.allowed_forward, gc);
        let expected_total = expected_motif_count(&motif, gc, 1000);
        let positions = 1000 - motif.allowed_forward.len() as u64 + 1;
        assert!((expected_total - expected_single * 2.0 * positions as f64).abs() < 1e-6);
    }

    #[test]
    fn iupac_complement_is_inverse() {
        assert_eq!(reverse_complement("ARY").unwrap(), "RYT");
    }
}
