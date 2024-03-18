use std::collections::HashMap;
use std::collections::HashSet;
use std::env;

extern crate clap;
use clap::*;

#[macro_use]
extern crate log;
use log::LevelFilter;
extern crate env_logger;
use env_logger::Builder;

use std::io;
use std::io::prelude::*;
use std::io::BufReader;

extern crate bio;
use bio::io::{fasta, fastq};

extern crate rand;
use rand::{seq::IteratorRandom, thread_rng};

mod translate;

fn main() {
    let mut app = build_cli();
    let matches = app.clone().get_matches();

    match matches.subcommand_name() {
        Some("concatenate") => {
            let m = matches.subcommand_matches("concatenate").unwrap();
            set_log_level(m);

            let mut alignment_files = vec![];
            for v in m.values_of("alignment-files").unwrap() {
                alignment_files.push(v)
            }

            if m.is_present("match-ids") {
                // Expect all FASTA files to have a shared set of IDs

                // Read all files to get a list of the IDs, and save seqs
                let mut marker_to_id_to_aln = HashMap::new();
                let mut marker_to_len = HashMap::new();
                let mut all_ids = HashSet::new();

                for alignment_file in alignment_files {
                    let reader = fasta::Reader::new(
                        std::fs::File::open(alignment_file)
                            .expect(&format!("Failed to open alignment file {}", alignment_file)),
                    );

                    let marker = alignment_file;
                    let mut marker_hash_map = HashMap::new();
                    for record in reader.records() {
                        let r = record.unwrap();
                        let id = r.id().to_string().clone();
                        let seq = std::str::from_utf8(r.seq()).unwrap().to_string();

                        assert!(
                            !marker_hash_map.contains_key(&id),
                            "Found duplicate sequence ID {} in alignment {}",
                            id,
                            alignment_file
                        );
                        marker_hash_map.insert(id.clone(), seq.clone());
                        all_ids.insert(id.clone());
                        if marker_to_len.contains_key(&marker) {
                            assert!(
                                marker_to_len.get(&marker).unwrap() == &seq.len(),
                                "Unexpected length {} for seuqence {} in alignment {}",
                                seq.len(),
                                &id,
                                &alignment_file
                            );
                        } else {
                            marker_to_len.insert(marker, seq.len());
                        }
                    }
                    marker_to_id_to_aln.insert(marker, marker_hash_map);
                }

                let mut id_to_full_seq: HashMap<String, String> = HashMap::new();
                for (marker, id_to_seq) in marker_to_id_to_aln.iter() {
                    for id in &all_ids {
                        let seq = match id_to_seq.get(id) {
                            Some(seq) => (*seq).clone(),
                            None => "-".repeat(*marker_to_len.get(marker).unwrap()),
                        };
                        match id_to_full_seq.get_mut(id) {
                            Some(prev_seq) => *prev_seq = format!("{}{}", prev_seq, seq),
                            None => {
                                id_to_full_seq.insert(id.clone(), seq);
                            }
                        }
                    }
                }

                for (id, seq) in id_to_full_seq {
                    println!(">{}\n{}", id, seq);
                }
            } else {
                // Expect all FASTA files to have the same num seqs, but don't require IDs to match
                let first_reader =
                    fasta::Reader::new(std::fs::File::open(alignment_files[0]).expect(&format!(
                        "Failed to open alignment file {}",
                        alignment_files[0]
                    )));
                let other_readers = alignment_files[1..]
                    .iter()
                    .map(|f| {
                        fasta::Reader::new(
                            std::fs::File::open(f)
                                .expect(&format!("Failed to open alignment file {}", f)),
                        )
                    })
                    .collect::<Vec<_>>();
                info!("Opened {} FASTA readers", other_readers.len() + 1);

                let first_record_iter = first_reader.records();
                let mut other_record_iters = vec![];
                for reader in other_readers {
                    other_record_iters.push(reader.records())
                }

                let mut first = true;
                let mut sequence_lengths = vec![];
                for record1 in first_record_iter {
                    let r1 = record1.unwrap();
                    let mut ids = vec![r1.id().to_string()];
                    let mut descriptions = vec![r1.desc().unwrap_or("").to_string()];
                    let mut seqs = vec![std::str::from_utf8(r1.seq()).unwrap().to_string()];

                    if first {
                        sequence_lengths.push(r1.seq().len());
                    } else {
                        assert!(r1.seq().len() == sequence_lengths[0]);
                    }

                    for (i, ref mut current_record_iter) in
                        other_record_iters.iter_mut().enumerate()
                    {
                        let current_record_res = current_record_iter.next().expect(&format!(
                            "Unexpected number of records in {}",
                            alignment_files[i + 1]
                        ));
                        let current_record = current_record_res.unwrap();

                        if first {
                            sequence_lengths.push(current_record.seq().len());
                        } else {
                            debug!("Sequence lengths: {:?}", sequence_lengths);
                            assert!(
                                current_record.seq().len() == sequence_lengths[i + 1],
                                "Unexpected length in {} in alignment #{}: {}, expected {}",
                                current_record.id(),
                                i + 2,
                                current_record.seq().len(),
                                sequence_lengths[i + 1]
                            );
                        }

                        ids.push(current_record.id().to_string());
                        descriptions.push(current_record.desc().unwrap_or("").to_string());
                        seqs.push(
                            std::str::from_utf8(current_record.seq())
                                .unwrap()
                                .to_string(),
                        );
                    }

                    if first {
                        first = false;
                    }

                    info!("Pairing sequences {}", ids.join(", "));

                    print!(
                        ">{} {}\n{}\n",
                        ids.join(" "),
                        descriptions.join(" "),
                        seqs.join(""),
                    );
                }
            }
        }
        Some("gff_to_fasta") => {
            let m = matches.subcommand_matches("gff_to_fasta").unwrap();
            set_log_level(m);

            let reader = BufReader::new(io::stdin());
            let mut is_passed_header = false;
            for line_res in reader.lines() {
                let line = line_res.expect("Line read fail");
                if is_passed_header {
                    println!("{}", line);
                } else if line == "##FASTA" {
                    is_passed_header = true;
                }
            }
            if !is_passed_header {
                error!("Could not fing '##FASTA' line in input GFF3, so failed");
                std::process::exit(1);
            }
        }
        Some("describe") => {
            let m = matches.subcommand_matches("describe").unwrap();
            set_log_level(m);

            let reader = fasta::Reader::new(io::stdin());
            for record in reader.records() {
                let res = record.unwrap();
                println!("{}\t{}", res.id(), res.seq().len())
            }
        }
        Some("gc") => {
            let m = matches.subcommand_matches("gc").unwrap();
            set_log_level(m);

            let mut num_gc: u64 = 0;
            let mut num_at: u64 = 0;
            let mut num_other: u64 = 0;

            // Poor coding here, duplicating. But annoying to fix
            if m.is_present("fasta") {
                let reader = fasta::Reader::new(io::stdin());

                for record in reader.records() {
                    for c in record.unwrap().seq() {
                        match *c as char {
                            'G' => num_gc += 1,
                            'g' => num_gc += 1,
                            'C' => num_gc += 1,
                            'c' => num_gc += 1,

                            'A' => num_at += 1,
                            'a' => num_at += 1,
                            'T' => num_at += 1,
                            't' => num_at += 1,

                            _ => num_other += 1,
                        }
                    }
                }
            } else {
                let reader = fastq::Reader::new(io::stdin());

                for record in reader.records() {
                    for c in record.unwrap().seq() {
                        match *c as char {
                            'G' => num_gc += 1,
                            'g' => num_gc += 1,
                            'C' => num_gc += 1,
                            'c' => num_gc += 1,

                            'A' => num_at += 1,
                            'a' => num_at += 1,
                            'T' => num_at += 1,
                            't' => num_at += 1,

                            _ => num_other += 1,
                        }
                    }
                }
            };

            if num_at + num_gc == 0 {
                panic!("No A, T, G or Cs found!")
            }
            println!(
                "{} {} {} {}",
                num_gc,
                num_at,
                num_other,
                num_gc as f64 / (num_at + num_gc) as f64
            );
        }
        Some("fasta_to_fastq") => {
            let m = matches.subcommand_matches("fasta_to_fastq").unwrap();
            set_log_level(m);

            let reader = fasta::Reader::new(io::stdin());
            for record_res in reader.records() {
                let record = record_res.expect("Failed to parse FASTA entry");
                println!("@{}", record.id());
                println!("{}", std::str::from_utf8(record.seq()).unwrap());
                println!("+");
                for _ in record.seq().iter().enumerate() {
                    print!("A");
                }
                println!();
            }
        }
        Some("add_genome_name_to_contig") => {
            let m = matches
                .subcommand_matches("add_genome_name_to_contig")
                .unwrap();
            set_log_level(m);

            let genome_fasta_files = m.values_of("genome-fasta-files").unwrap();
            info!(
                "Read in {} genome fasta file paths",
                genome_fasta_files.len()
            );
            let output_dir = m.value_of("output-directory").unwrap();

            for genome in genome_fasta_files {
                let reader = fasta::Reader::new(
                    std::fs::File::open(genome)
                        .expect(&format!("Failed to open genome fasta file {}", genome)),
                );
                let genome_name = std::path::Path::new(&genome)
                    .file_stem()
                    .expect("Failed to parse genome name")
                    .to_str()
                    .expect("Character conversion issue");
                let outpath = format!("{}/{}.renamed_contigs.fna", output_dir, genome_name);
                debug!("Writing output fasta file {} ..", &outpath);
                let mut writer = fasta::Writer::to_file(&outpath)
                    .expect(&format!("Failed to create output file {}", &outpath));
                let mut count: usize = 0;
                for record_res in reader.records() {
                    let record = record_res.expect("Failed to parse FASTA entry");
                    writer
                        .write(
                            &format!("{}~{}", genome_name, record.id()),
                            record.desc(),
                            record.seq(),
                        )
                        .expect("Failed to write FASTA entry");
                    count += 1;
                }
                debug!("Wrote {} sequences", count);
            }
        },
        Some("msa-sample") => {
            let m = matches.subcommand_matches("msa-sample").unwrap();
            set_log_level(m);

            // Read sequences in
            let mut original_order = vec![];
            let mut id_to_seq = HashMap::new();
            let mut sequence_length = None;
            let reader = fasta::Reader::new(
                std::fs::File::open(m.value_of("fasta").unwrap()).expect(&format!(
                    "Failed to open alignment file {}",
                    m.value_of("fasta").unwrap()
                )),
            );
            for record in reader.records() {
                let r = record.unwrap();
                let id = r.id().to_string().clone();
                let seq = std::str::from_utf8(r.seq()).unwrap().to_string();

                assert!(
                    !id_to_seq.contains_key(&id),
                    "Found duplicate sequence ID {} in alignment",
                    id,
                );
                original_order.push(id.clone());
                id_to_seq.insert(id.clone(), seq.clone());
                match sequence_length {
                    Some(l) => {
                        assert!(seq.len() == l, "Found unexpected length in sequence {}", id);
                    }
                    None => {
                        sequence_length = Some(seq.len());
                    }
                }
            }

            // Determine which columns are usable
            let mut useable_columns = vec![];
            let num_seqs = id_to_seq.len();
            let num_seqs_threshold = (num_seqs as f32)
                * m.value_of("min-perc-taxa")
                    .unwrap()
                    .parse::<f32>()
                    .expect("Failed to parse float for min-perc-taxa")
                / 100.;
            info!("Requiring {} positions to be filled", num_seqs_threshold);
            for i in 0..sequence_length.unwrap() {
                let mut num_not_gap = 0;
                for seq in id_to_seq.values() {
                    let ch = seq.chars().nth(i).unwrap();
                    debug!("Found char {}", ch);
                    if ch != '-' && ch != '*' {
                        num_not_gap += 1;
                    }
                }
                if num_not_gap as f32 >= num_seqs_threshold {
                    useable_columns.push(i)
                }
            }
            info!("Found {} columns passing thresholds", useable_columns.len());

            // Sample useable column IDs
            let num_to_choose = m
                .value_of("num-columns")
                .unwrap()
                .parse::<usize>()
                .expect("faild to parse num-columns");
            let mut rng = thread_rng();
            if useable_columns.len() > num_to_choose {
                useable_columns = useable_columns
                    .iter()
                    .choose_multiple(&mut rng, num_to_choose)
                    .into_iter()
                    .map(|r| *r)
                    .collect();
            }

            // Print sampled columns
            info!("Printing {} columns", useable_columns.len());
            for id in original_order {
                let seq = id_to_seq.get(&id).unwrap();
                println!(">{}", id);
                for i in &useable_columns {
                    print!("{}", seq.chars().nth(*i).unwrap());
                }
                println!()
            }
        },
        // Some("contig_dereplicate") => {
        //     panic!();
        //     // For each line
        //         if subject name >= query name
        //             continue
        //         else
        //             Query should be masked in this region
        //             if %id is good enough
        //                 add it to store of masks (just a pipe to sort?)
        //     Sort tab-separated file of masks by sequence name
        //     throw out stretches < 100bp
        //     Add 20bp buffer
        //     print fasta if any unmasked regions
        //     masking algorithm:
        //     foreach contig, make a priority queue of starts/stops with position
        //     iterate queue, changing state if necessary
        //     give last chunk back if open at end
        // }
        Some("backtranslate") => {
            let m = matches.subcommand_matches("backtranslate").unwrap();
            set_log_level(m);

            let reader = fasta::Reader::new(
                std::fs::File::open(m.value_of("input").unwrap()).expect(&format!(
                    "Failed to open input file {}",
                    m.value_of("input").unwrap()
                )),
            );
            let mut writer = fasta::Writer::to_file(m.value_of("output").unwrap())
                .expect(&format!("Failed to open output file {}", m.value_of("output").unwrap()));
            for record_res in reader.records() {
                let record = record_res.expect("Failed to parse FASTA entry");
                let seq = std::str::from_utf8(record.seq()).unwrap();
                let mut new_seq = String::new();
                for aa in seq.chars() {
                    match aa {
                        // Ala, A	GCT, GCC, GCA, GCG	GCN	Ile, I	ATT, ATC, ATA	ATH
                        'A' => new_seq.push_str("GCT"),
                        'I' => new_seq.push_str("ATT"),
                        // Arg, R	CGT, CGC, CGA, CGG; AGA, AGG	CGN, AGR; or
                        // CGY, MGR	Leu, L	CTT, CTC, CTA, CTG; TTA, TTG	CTN, TTR; or
                        // CTY, YTR
                        'R' => new_seq.push_str("CGT"),
                        'L' => new_seq.push_str("CTT"),
                        // Asn, N	AAT, AAC	AAY	Lys, K	AAA, AAG	AAR
                        'N' => new_seq.push_str("AAT"),
                        'K' => new_seq.push_str("AAA"),
                        // Asp, D	GAT, GAC	GAY	Met, M	ATG
                        'D' => new_seq.push_str("GAT"),
                        'M' => new_seq.push_str("ATG"),
                        // Asn or Asp, B	AAT, AAC; GAT, GAC	RAY	Phe, F	TTT, TTC	TTY
                        'B' => new_seq.push_str("AAT"),
                        'F' => new_seq.push_str("TTT"),
                        // Cys, C	TGT, TGC	TGY	Pro, P	CCT, CCC, CCA, CCG	CCN
                        'C' => new_seq.push_str("TGT"),
                        'P' => new_seq.push_str("CCT"),
                        // Gln, Q	CAA, CAG	CAR	Ser, S	TCT, TCC, TCA, TCG; AGT, AGC	TCN, AGY
                        'Q' => new_seq.push_str("CAA"),
                        'S' => new_seq.push_str("TCT"),
                        // Glu, E	GAA, GAG	GAR	Thr, T	ACT, ACC, ACA, ACG	ACN
                        'E' => new_seq.push_str("GAA"),
                        'T' => new_seq.push_str("ACT"),
                        // Gln or Glu, Z	CAA, CAG; GAA, GAG	SAR	Trp, W	TGG
                        'Z' => new_seq.push_str("CAA"),
                        'W' => new_seq.push_str("TGG"),
                        // Gly, G	GGT, GGC, GGA, GGG	GGN	Tyr, Y	TAT, TAC	TAY
                        'G' => new_seq.push_str("GGT"),
                        'Y' => new_seq.push_str("TAT"),
                        // His, H	CAT, CAC	CAY	Val, V	GTT, GTC, GTA, GTG	GTN
                        'H' => new_seq.push_str("CAT"),
                        'V' => new_seq.push_str("GTT"),
                        // START	ATG, CTG, UTG	HTG	STOP	TAA, TGA, TAG	TRA, TAR
                        '*' => new_seq.push_str("TAA"),
                        'X' => new_seq.push_str("NNN"),
                        'U' => new_seq.push_str("NNN"),
                        _ => panic!("Unexpected amino acid {}", aa),
                    }
                }
                writer
                    .write(record.id(), record.desc(), new_seq.as_bytes())
                    .expect("Failed to write FASTA entry");
            }
        },
        Some("translate") => {
            let m = matches.subcommand_matches("translate").unwrap();
            set_log_level(m);

            let reader = fasta::Reader::new(
                std::fs::File::open(m.value_of("input").unwrap()).expect(&format!(
                    "Failed to open input file {}",
                    m.value_of("input").unwrap()
                )),
            );
            let mut writer = fasta::Writer::to_file(m.value_of("output").unwrap())
                .expect(&format!("Failed to open output file {}", m.value_of("output").unwrap()));
            for record_res in reader.records() {
                let record = record_res.expect("Failed to parse FASTA entry");
                let seq = std::str::from_utf8(record.seq()).unwrap();
                let new_seq = translate::translate(seq.as_bytes());
                writer
                    .write(record.id(), record.desc(), new_seq.as_bytes())
                    .expect("Failed to write FASTA entry");
            }
        },
        _ => {
            app.print_help().unwrap();
            println!();
        }
    }
}

fn set_log_level(matches: &clap::ArgMatches) {
    let mut log_level = LevelFilter::Info;
    if matches.is_present("verbose") {
        log_level = LevelFilter::Debug;
    }
    if matches.is_present("quiet") {
        log_level = LevelFilter::Error;
    }
    let mut builder = Builder::new();
    builder.filter(None, log_level);
    if env::var("RUST_LOG").is_ok() {
        builder.parse_filters(&env::var("RUST_LOG").unwrap());
    }
    if builder.try_init().is_err() {
        panic!("Failed to set log level - has it been specified multiple times?")
    }
}

fn build_cli() -> App<'static, 'static> {
    return App::new("bbbin")
        .version(crate_version!())
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Utilities for bioinformtics")
        .args_from_usage("-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
        .subcommand(
            SubCommand::with_name("gc")
                .about("Calculate G+C content of FASTA/Q sequences piped in")
                .arg(Arg::with_name("fasta")
                .help("Input is FASTA not FASTQ")
                .long("fasta")))
        .subcommand(
            SubCommand::with_name("describe")
                .about("Calculate length of each FASTA sequence"))
        .subcommand(
            SubCommand::with_name("fasta_to_fastq")
                .about("Make a FASTQ from a FASTA file, setting all quality values to 'A'"))
        .subcommand(
            SubCommand::with_name("gff_to_fasta")
                .about("Extract the master FASTA record from the bottom of the gff3"))
        .subcommand(
            SubCommand::with_name("concatenate")
                .about("Concatenate fasta files by pasting sequences together, for concatenated gene phylogeny")
                .arg(Arg::with_name("alignment-files")
                    .help("Files to concatenate")
                    .long("alignment-files")
                    .required(true)
                    .takes_value(true)
                    .multiple(true))
                .arg(Arg::with_name("match-ids")
                    .help("Match sequences together based on ID rather than position [default: Match on position]")
                    .long("match-ids")
                ))
        .subcommand(
            SubCommand::with_name("add_genome_name_to_contig")
                .about("Rename the contigs in a genome-wise fasta file")
                .arg(Arg::with_name("genome-fasta-files")
                    .long("genome-fasta-files")
                    .required(true)
                    .takes_value(true)
                    .multiple(true)
                    .help("List of genome fasta files to rename"))
                .arg(Arg::with_name("output-directory")
                    .long("output-directory")
                    .required(true)
                    .takes_value(true)
                    .help("Folder to write new genome fasta files to"))
        )
        .subcommand(
            SubCommand::with_name("msa-sample")
                .about("Sample columns from a multiple sequence alignment")
                .arg(Arg::with_name("fasta")
                    .long("fasta")
                    .required(true)
                    .takes_value(true)
                    .help("Fasta file to sample (a multiple sequence alignment)"))
                .arg(Arg::with_name("num-columns")
                    .long("num-columns")
                    .required(true)
                    .takes_value(true)
                    .help("Number of columns to choose (without replacement)"))
                .arg(Arg::with_name("min-perc-taxa")
                    .long("min-perc-taxa")
                    .takes_value(true)
                    .default_value("50")
                    .help("Minimum percentage of taxa required to retain column (inclusive bound) (default: 50)"))
        )
        .subcommand(
            SubCommand::with_name("backtranslate")
                .about("Backtranslate a protein sequence to a DNA sequence, guessing the codons")
                .arg(Arg::with_name("input")
                    .long("input")
                    .required(true)
                    .takes_value(true)
                    .help("Input file"))
                .arg(Arg::with_name("output")
                    .long("output")
                    .required(true)
                    .takes_value(true)
                    .help("Output file"))
        )
        .subcommand(
            SubCommand::with_name("translate")
                .about("Translate DNA sequences to protein sequences")
                .arg(Arg::with_name("input")
                    .long("input")
                    .required(true)
                    .takes_value(true)
                    .help("Input file"))
                .arg(Arg::with_name("output")
                    .long("output")
                    .required(true)
                    .takes_value(true)
                    .help("Output file"))
        );
}
