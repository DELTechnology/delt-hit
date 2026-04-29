use std::cmp::Reverse;
use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Instant;

type SelectionIds = Vec<u16>;
type BarcodeIds = Vec<u16>;
type SelectionCounts = HashMap<BarcodeIds, u32>;
type Counts = HashMap<SelectionIds, SelectionCounts>;

fn parse_line(line: &str) -> Option<(SelectionIds, BarcodeIds)> {
    let mut parts = line.trim_end().split('?');
    parts.next()?;

    let mut selection_ids = Vec::new();
    let mut barcode_ids = Vec::new();

    for token in parts {
        let value = token.rsplit('.').next()?.parse::<u16>().ok()?;
        if token.contains('S') {
            selection_ids.push(value);
        } else if token.contains('B') {
            barcode_ids.push(value);
        }
    }

    Some((selection_ids, barcode_ids))
}

fn selection_name(selection_ids: &[u16]) -> String {
    selection_ids
        .iter()
        .map(u16::to_string)
        .collect::<Vec<_>>()
        .join("-")
}

fn id_name(barcodes: &[u16]) -> String {
    barcodes
        .iter()
        .map(u16::to_string)
        .collect::<Vec<_>>()
        .join("_")
}

fn write_counts(output_dir: &Path, counts: Counts) -> io::Result<()> {
    fs::create_dir_all(output_dir)?;

    for (selection_ids, barcode_counts) in counts {
        let mut rows: Vec<(BarcodeIds, u32)> = barcode_counts.into_iter().collect();
        rows.sort_by_key(|(barcodes, count)| (Reverse(*count), barcodes.clone()));

        let file_name = format!("{}_counts.txt", selection_name(&selection_ids));
        let file_path = output_dir.join(file_name);
        let file = File::create(file_path)?;
        let mut writer = BufWriter::new(file);

        let num_codes = rows.first().map(|(barcodes, _)| barcodes.len()).unwrap_or(0);
        for idx in 0..num_codes {
            write!(writer, "code_{}\t", idx)?;
        }
        writeln!(writer, "count\tid")?;

        for (barcodes, count) in rows {
            for value in &barcodes {
                write!(writer, "{}\t", value)?;
            }
            writeln!(writer, "{}\t{}", count, id_name(&barcodes))?;
        }
    }

    Ok(())
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: rust_get_counts <reads_with_adapters.gz> <output_dir>");
        std::process::exit(1);
    }

    let input_path = PathBuf::from(&args[1]);
    let output_dir = PathBuf::from(&args[2]);

    let start = Instant::now();

    let mut child = Command::new("gzip")
        .arg("-dc")
        .arg(&input_path)
        .stdout(Stdio::piped())
        .spawn()?;

    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| io::Error::other("failed to capture gzip stdout"))?;
    let reader = BufReader::new(stdout);

    let mut counts: Counts = HashMap::new();
    let mut n_reads: u64 = 0;

    for line in reader.lines() {
        let line = line?;
        if let Some((selection_ids, barcode_ids)) = parse_line(&line) {
            let barcode_counts = counts.entry(selection_ids).or_default();
            *barcode_counts.entry(barcode_ids).or_insert(0) += 1;
            n_reads += 1;
        }
    }

    let status = child.wait()?;
    if !status.success() {
        return Err(io::Error::other(format!("gzip exited with status {status}")));
    }

    write_counts(&output_dir, counts)?;

    eprintln!(
        "processed_reads={n_reads} elapsed_s={:.2}",
        start.elapsed().as_secs_f64()
    );

    Ok(())
}
