use std::{
    fs::File,
    io::{BufRead, BufReader, Error, ErrorKind, Lines, Read},
    path::Path,
};

use zoe::{
    data::{err::ResultWithErrorContext, fasta::FastaSeq},
    unwrap_or_return_some_err,
};

pub(crate) struct TsvReader<R> {
    lines: Lines<BufReader<R>>,
}

impl TsvReader<std::fs::File> {
    /// Creates an iterator over the FASTA data contained in a path, using a
    /// buffered reader.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the path does not exist, if there are insufficient
    /// permissions to read from it, or if it contains no data. The path is
    /// included in the error message.
    pub fn from_path<P>(path: P) -> std::io::Result<TsvReader<File>>
    where
        P: AsRef<Path>, {
        let path = path.as_ref();
        let file = File::open(path).with_path_context("Failed to open path", path)?;
        Ok(Self::from_readable(file).with_path_context("Failed to read data at path", path)?)
    }
}

impl<R: Read> TsvReader<R> {
    fn from_readable(reader: R) -> std::io::Result<Self> {
        Self::from_bufreader(BufReader::new(reader))
    }

    /// Creates an iterator over TSV data from a [`BufReader`].
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    fn from_bufreader(mut reader: BufReader<R>) -> std::io::Result<Self> {
        if reader.fill_buf()?.is_empty() {
            return Err(std::io::Error::new(ErrorKind::InvalidData, "No TSV data was found!"));
        }

        Ok(TsvReader { lines: reader.lines() })
    }

    fn parse_line(line: &str) -> Result<FastaSeq, Error> {
        let mut columns = line.split('\t').map(|part| part.trim_ascii().to_string());

        let name = columns.next().unwrap_or_default();
        if name.is_empty() {
            return Err(std::io::Error::new(
                ErrorKind::InvalidData,
                "Invalid TSV: sequence name (first column) cannot be empty",
            ));
        }
        let Some(second) = columns.next() else {
            return Err(std::io::Error::new(
                ErrorKind::InvalidData,
                "Invalid TSV format: expect 2 tab-separated columns but found 1",
            ));
        };
        let sequence = second.as_bytes().to_vec();
        if sequence.is_empty() {
            return Err(std::io::Error::new(
                ErrorKind::InvalidData,
                "Invalid TSV: sequence (second column) cannot be empty",
            ));
        }
        Ok(FastaSeq { name, sequence })
    }
}

impl<R: Read> Iterator for TsvReader<R> {
    type Item = Result<FastaSeq, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = unwrap_or_return_some_err!(self.lines.next()?);

            if !line.trim().is_empty() {
                return Some(Self::parse_line(&line));
            }
        }
    }
}
