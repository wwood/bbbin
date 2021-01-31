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

            let first_reader = fasta::Reader::new(
                std::fs::File::open(alignment_files[0])
                    .expect(&format!("Failed to open alignment file {}", alignment_files[0]))
            );
            let other_readers = alignment_files[1..].iter().map(|f| 
                fasta::Reader::new(
                    std::fs::File::open(f)
                        .expect(&format!("Failed to open alignment file {}", f))
                    )
                )
                .collect::<Vec<_>>();
            info!("Opened {} FASTA readers", other_readers.len()+1);

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

                for (i, ref mut current_record_iter) in other_record_iters.iter_mut().enumerate() {
                    let current_record_res = current_record_iter.next().expect(&format!("Unexpected number of records in {}", alignment_files[i+1]));
                    let current_record = current_record_res.unwrap();

                    if first {
                        sequence_lengths.push(current_record.seq().len());
                    } else {
                        debug!("Sequence lengths: {:?}", sequence_lengths);
                        assert!(
                            current_record.seq().len() == sequence_lengths[i+1], 
                            "Unexpected length in {} in alignment #{}: {}, expected {}",
                            current_record.id(), i+2, current_record.seq().len(), sequence_lengths[i+1]
                        );
                    }

                    ids.push(current_record.id().to_string());
                    descriptions.push(current_record.desc().unwrap_or("").to_string());
                    seqs.push(std::str::from_utf8(current_record.seq()).unwrap().to_string());
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

            let reader = fastq::Reader::new(io::stdin());
            let mut num_gc: u64 = 0;
            let mut num_at: u64 = 0;
            let mut num_other: u64 = 0;

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
        }
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
                     .about("Calculate G+C content of FASTQ sequences piped in"))
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
                    .multiple(true)))
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
        );
}
