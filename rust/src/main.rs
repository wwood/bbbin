use std::env;

extern crate clap;
use clap::*;

#[macro_use]
extern crate log;
use log::LogLevelFilter;
extern crate env_logger;
use env_logger::LogBuilder;

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

            let reader1 = fasta::Reader::new(
                std::fs::File::open(
                    m.value_of("ALIGNMENT1")
                        .expect("Failed to open first FASTA file"),
                )
                .expect("Failed to open first FASTA file"),
            );
            let reader2 = fasta::Reader::new(
                std::fs::File::open(
                    m.value_of("ALIGNMENT2")
                        .expect("Failed to open first FASTA file"),
                )
                .expect("Failed to open first FASTA file"),
            );
            let mut num_residues = vec![0usize; 2];

            let mut first = true;
            for (record1, record2) in reader1.records().zip(reader2.records()) {
                let r1 = record1.unwrap();
                let r2 = record2.unwrap();

                let len1 = r1.seq().len();
                let len2 = r2.seq().len();

                if first {
                    first = false;
                    num_residues[0] = len1;
                    num_residues[1] = len2;
                } else {
                    assert!(num_residues[0] == len1);
                    assert!(num_residues[1] == len2);
                }

                eprintln!("Pairing sequence {} with {}", r1.id(), r2.id());
                print!(
                    ">{} {}\n{}{}\n",
                    r1.id(),
                    r2.id(),
                    std::str::from_utf8(r1.seq()).unwrap(),
                    std::str::from_utf8(r2.seq()).unwrap()
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
    let mut log_level = LogLevelFilter::Info;
    if matches.is_present("verbose") {
        log_level = LogLevelFilter::Debug;
    }
    if matches.is_present("quiet") {
        log_level = LogLevelFilter::Error;
    }
    let mut builder = LogBuilder::new();
    builder.filter(None, log_level);
    if env::var("RUST_LOG").is_ok() {
        builder.parse(&env::var("RUST_LOG").unwrap());
    }
    builder.init().unwrap();
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
                .arg(Arg::with_name("ALIGNMENT1")
                .help("Sets the 1st input file to use")
                .required(true)
                .index(1))
                .arg(Arg::with_name("ALIGNMENT2")
                .help("Sets the 2nd input file to use")
                .required(true)
                .index(2)))
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
