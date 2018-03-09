extern crate bbbin;

use std::env;
use std::str;

extern crate clap;
use clap::*;

extern crate log;
use log::LevelFilter;
extern crate env_logger;
use env_logger::Builder;

fn main(){
    let mut app = build_cli();
    let matches = app.clone().get_matches();

    match matches.subcommand_name() {
        Some("gc") => {
            let m = matches.subcommand_matches("gc").unwrap();
            set_log_level(m);
            let sequence_files: Vec<&str> = m.values_of("sequences").unwrap().collect();
            let interval = value_t!(m.value_of("interval"), usize).unwrap();
            bbbin::gc(&sequence_files, interval, &mut std::io::stdout())
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
        builder.parse(&env::var("RUST_LOG").unwrap());
    }
    builder.init();
}

fn build_cli() -> App<'static, 'static> {
    //-f, --fasta-files=<FILE>...         'Read contig to genome mapping from these fasta files'
    let genome_args: &'static str = "-s, --sequences=<FASTA>...      'Fasta file to calculate G+C for'
                      -i, --interval=<INTEGER>         'Calculate G+C for each interval in the sequence of this size'

                      -v, --verbose       'Print extra debug logging information'
                      -q, --quiet         'Unless there is an error, do not print logging information'";

    return App::new("bbbin")
        .version("0.1.0-pre")
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Ben's bioinformatic bin - various utilities")
        .args_from_usage("-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
        .subcommand(
            SubCommand::with_name("gc")
                .about("Calculate G+C%")
                .args_from_usage(&genome_args));
}
