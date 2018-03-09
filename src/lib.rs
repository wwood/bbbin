extern crate bio;
extern crate log;

use std::str;

pub fn gc(
    fastas: &Vec<&str>, interval: usize, print_stream: &mut std::io::Write,){

    for fasta in fastas {
        let reader = bio::io::fasta::Reader::from_file(fasta).expect(
            &format!("Unable to open file {}", fasta));
        for record in reader.records() {
            let r = record.expect("Failed to fetch next sequence");
            let id = r.id();
            let seq = r.seq();
            let current_length = seq.len();
            let mut offset: usize = 0;
            for s in seq.chunks(interval) {
                if s.len() == interval {
                    let mut at: u32 = 0;
                    let mut total: u32 = 0;
                    for c in s {
                        if *c == 'A' as u8 ||
                            *c == 'a' as u8 ||
                            *c == 'T' as u8 ||
                            *c == 't' as u8 {
                                at += 1;
                                total += 1;
                            } else if *c == 'C' as u8 ||
                            *c == 'c' as u8 ||
                            *c == 'G' as u8 ||
                            *c == 'g' as u8 {
                                total += 1;
                            }
                    };
                    writeln!(print_stream, "{}\t{}\t{}\t{}\t{}",
                             id, offset+1, offset+interval, at, total)
                        .expect("Failed to output result line");
                    offset += interval;
                }
            }
            while offset+interval < current_length {
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    extern crate tempfile;
    use std::io::prelude::*;

    #[test]
    fn it_works() {
        let mut file = tempfile::NamedTempFile::new()
            .expect("failed to create tempfile");
        writeln!(file, ">s1 10 chars");
        writeln!(file, "ATGCAAAAAA");
        writeln!(file, ">s3 N char");
        writeln!(file, "ATGCAANAAATTT");
        file.flush().unwrap();
        let mut stream = Cursor::new(Vec::new());
        gc(&vec![file.path().to_str().unwrap()], 5, &mut stream);
        assert_eq!(
            "s1\t1\t5\t3\t5\ns1\t6\t10\t5\t5\ns3\t1\t5\t3\t5\ns3\t6\t10\t4\t4\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }
}
