extern crate assert_cli;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;

    #[test]
    fn test_backtranslate() {
        Assert::main_binary()
            .with_args(&["backtranslate", "--input", "tests/data/backtranslate/random2_3aa.fna", "--output", "/dev/stdout"])
            .succeeds()
            .stdout()
            .is(">random_sequence_length_3_1\n\
                CAATTTATG\n\
                >random_sequence_length_3_2\n\
                TAATAAAAA\n")
            .unwrap();
    }

    #[test]
    fn test_translate() {
        Assert::main_binary()
            .with_args(&["translate", "--input", "tests/data/translate/random2_3aa.fna", "--output", "/dev/stdout"])
            .succeeds()
            .stdout()
            .is(">random_sequence_length_3_1\n\
                MKK\n")
            .unwrap();
    }
}
