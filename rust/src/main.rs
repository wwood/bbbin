use std::env;

extern crate clap;
use clap::*;

#[macro_use]
extern crate log;
use log::LogLevelFilter;
extern crate env_logger;
use env_logger::LogBuilder;

use std::io;

extern crate bio;
use bio::io::fastq;

fn main() {
    let mut app = build_cli();
    let matches = app.clone().get_matches();

    match matches.subcommand_name() {
        Some("gc") => {
            let m = matches.subcommand_matches("gc").unwrap();
            set_log_level(m);

            let reader = fastq::Reader::new(io::stdin());
            let mut num_gc: u32 = 0;
            let mut num_at: u32 = 0;
            let mut num_other: u32 = 0;

            for record in reader.records() {
                for c in record.unwrap().seq() {
                    match *c as char {
                        'G' => {num_gc += 1},
                        'g' => {num_gc += 1},
                        'C' => {num_gc += 1},
                        'c' => {num_gc += 1},

                        'A' => {num_at += 1},
                        'a' => {num_at += 1},
                        'T' => {num_at += 1},
                        't' => {num_at += 1},

                        _ => {num_other += 1}
                    }
                }
            }
            if num_at + num_gc == 0 {
                panic!("No A, T, G or Cs found!")
            }
            println!("{} {} {} {}",
                     num_gc, num_at, num_other, num_gc as f32 / (num_at+num_gc) as f32);
        },
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
}
